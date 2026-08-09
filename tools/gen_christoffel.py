#!/usr/bin/env python3
"""Generate the christoffel*() bodies of OpenACC_cli/cuda_ray.h.

Emits the geodesic right-hand side  ch[i] = -Gamma^i_jk v^j v^k  for the
Kerr-Newman metric in Boyer-Lindquist coordinates, as a straight-line block of
FP temporaries. Two specializations are generated, plus a dispatcher:

  christoffel_general  a != 0, Q arbitrary (full Kerr-Newman).
  christoffel_static   a == 0 (Reissner-Nordstrom if Q != 0, Schwarzschild if
                        Q is also 0 - one body covers both, see below).
  christoffel          picks between them on hole.a == FP(0). hole.a is the
                        same value for every thread in a given render, so this
                        branch is warp-uniform - free on GPU, not divergent.

Why only two cases, not four: spin `a` is what couples t/phi and breaks
spherical symmetry, so a=0 collapses the metric to diagonal and spherically
symmetric - that alone takes the CSE'd op count from 216 to 52 (measured with
--bench-cases). Charge Q only ever appears in an isolated `rs*r - Q**2`
combination and never couples coordinates, so a Q=0-only ("pure Kerr")
specialization was measured and rejected: 223 ops, i.e. no better than
general. Zeroing Q on top of a=0 (true Schwarzschild) only takes 52 ops to 48,
not worth a third branch and a third body to keep verified.

The recipe for each case (run --bench for the numbers behind it):

  1. Contract with v BEFORE generating code.  Only the 4 contracted sums are
     ever needed, never the 64 individual Gamma^i_jk, so CSE sees a far smaller
     problem.
  2. Invert the metric analytically.  The (t,phi) 2x2 block has
     det = -Delta*sin^2(theta) exactly, so g^{ij} is written in closed form and
     the sin^2 cancels out of g^tt and g^tphi instead of being carried through
     CSE as an expanded polynomial.
  3. factor_terms, NOT simplify.  simplify() burns ~10 minutes and yields
     cos(4*theta)/cos(6*theta) multiple-angle terms that make the result
     *worse*.  cancel/together are catastrophic (~400 temps).
  4. cse(optimizations='basic'), then ccode with Pow -> pown.

Floating-point accuracy of the general case was checked against Herbie
(tools/herbie_christoffel.py, see its docstring for the verdict): already
well-conditioned, worst case ~2 bits of float32 error, which is below what the
dopri54 adaptive step controller's own tolerance can resolve. Herbie's
rewrites were not adopted - they optimize each ch[i] independently and would
cost ~5x the ops by losing the cross-output CSE sharing this file exists for.

Usage:
    python tools/gen_christoffel.py             # print the generated block
    python tools/gen_christoffel.py --bench      # compare pipelines (general case)
    python tools/gen_christoffel.py --verify     # numeric check vs raw form, each case
    python tools/gen_christoffel.py --check      # diff against cuda_ray.h

Requires sympy (in the root pixi environment: `pixi run python tools/...`).
"""
import argparse
import math
import random
import re
import time
from pathlib import Path

import sympy as sp
from sympy import sin, cos, symbols, Rational, cse, ccode

HEADER = Path(__file__).resolve().parent.parent / 'OpenACC_cli' / 'cuda_ray.h'
BEGIN_MARK = '// --- BEGIN GENERATED christoffel: tools/gen_christoffel.py ---'
END_MARK = '// --- END GENERATED christoffel ---'

a, Q, rs = symbols('a Q rs', real=True, positive=True)
r, th = symbols('r theta', real=True, positive=True)
v = symbols('v0 v1 v2 v3', real=True)
COORD = [None, r, th, None]   # only r and theta appear in g

# (function name, a substitution, Q substitution, does the body reference `a`)
CASES = [
    ('christoffel_general', a, Q, True),
    ('christoffel_static', sp.Integer(0), Q, False),
]


def covariant_metric(a_val, Q_val):
    """g_{mu nu}, signature (+,-,-,-), matching notebooks/derivation/."""
    delta = r**2 - rs*r + a_val**2 + Q_val**2
    rho2 = r**2 + a_val**2*cos(th)**2
    s2 = sin(th)**2
    g = sp.zeros(4, 4)
    g[0, 0] = delta/rho2 - a_val**2*s2/rho2
    g[1, 1] = -rho2/delta
    g[2, 2] = -rho2
    g[3, 3] = (delta*a_val**2*sin(th)**4 - a_val**4*s2 - 2*a_val**2*r**2*s2 - r**4*s2)/rho2
    g[0, 3] = g[3, 0] = (a_val**3*s2 + a_val*r**2*s2 - delta*a_val*s2)/rho2
    return g, delta, rho2, s2


def contravariant_metric_analytic(a_val, Q_val, delta, rho2, s2):
    """g^{mu nu} in closed form, using det(t,phi block) = -Delta*sin^2(theta)."""
    gi = sp.zeros(4, 4)
    A = (r**2 + a_val**2)**2 - delta*a_val**2*s2
    gi[0, 0] = A/(rho2*delta)
    gi[0, 3] = gi[3, 0] = a_val*(rs*r - Q_val**2)/(rho2*delta)
    gi[3, 3] = -(delta - a_val**2*s2)/(rho2*delta*s2)
    gi[1, 1] = -delta/rho2
    gi[2, 2] = -1/rho2
    return gi


def contravariant_metric_blockwise(g):
    """Same inverse, but letting sympy carry the expanded block determinant.

    Kept so --bench can show what the closed-form determinant buys.
    """
    gi = sp.zeros(4, 4)
    det2 = g[0, 0]*g[3, 3] - g[0, 3]**2
    gi[0, 0] = g[3, 3]/det2
    gi[3, 3] = g[0, 0]/det2
    gi[0, 3] = gi[3, 0] = -g[0, 3]/det2
    gi[1, 1] = 1/g[1, 1]
    gi[2, 2] = 1/g[2, 2]
    return gi


def geodesic_rhs(g, gi, simp=sp.factor_terms):
    """rhs[i] = -Gamma^i_jk v^j v^k, contracted before any codegen."""
    dg = [[[0]*4 for _ in range(4)] for _ in range(4)]
    for l in range(4):
        for k in range(4):
            for d in (1, 2):
                dg[d][l][k] = sp.diff(g[l, k], COORD[d])
    out = []
    for i in range(4):
        acc = 0
        for j in range(4):
            for k in range(4):
                Gam = 0
                for l in range(4):
                    if gi[i, l] == 0:
                        continue
                    term = dg[j][l][k] + dg[k][l][j] - dg[l][j][k]
                    if term != 0:
                        Gam += gi[i, l]*term
                if Gam != 0:
                    acc -= Rational(1, 2)*Gam*v[j]*v[k]
        out.append(simp(acc))
    return out


def rhs_for(a_val, Q_val, simp=sp.factor_terms):
    g, delta, rho2, s2 = covariant_metric(a_val, Q_val)
    gi = contravariant_metric_analytic(a_val, Q_val, delta, rho2, s2)
    return geodesic_rhs(g, gi, simp)


def tally(repl, red):
    """Op counts by class over the whole straight-line block."""
    add = mul = div = trans = pw = 0
    for e in [e for _, e in repl] + list(red):
        for node in sp.preorder_traversal(e):
            if isinstance(node, sp.Add):
                add += len(node.args) - 1
            elif isinstance(node, sp.Mul):
                mul += len(node.args) - 1
            elif isinstance(node, sp.Pow):
                _, ex = node.args
                if ex.is_Integer and ex < 0:
                    div += 1
                    pw += abs(int(ex)) - 1
                elif ex.is_Integer:
                    pw += int(ex) - 1
                else:
                    trans += 1
            elif isinstance(node, (sp.sin, sp.cos, sp.tan)):
                trans += 1
    return dict(temps=len(repl), add=add, mul=mul + pw, div=div, trans=trans,
                total=add + mul + pw + div + trans)


USER_FUNCS = {'Pow': [(lambda b, e: e.is_integer, 'pown'),
                      (lambda b, e: not e.is_integer, 'pow')]}


def to_c_symbols(rhs):
    """Rename r/theta/v_i to the array slots cuda_ray.h uses.

    Done before CSE, not after: the renaming changes the term ordering CSE
    sees, so measuring the symbolic form would not report the cost of the code
    actually emitted.
    """
    subs = {r: sp.Symbol('x[1]'), th: sp.Symbol('y')}
    for i in range(4):
        subs[v[i]] = sp.Symbol(f'v[{i}]')
    return [e.subs(subs) for e in rhs]


POLE_COMMENT = """\
// The phi equation carries an explicit 1/sin(theta).  This is a coordinate
// artifact of Boyer-Lindquist-type charts at the polar axis (theta=0,pi), not a
// physical curvature singularity: every numerator multiplying it is itself
// proportional to v[3]=dphi/de, which vanishes on the axis along with
// sin(theta), so the true ratio stays finite.  A ray integrated exactly through
// the axis would otherwise hit a 0/0-like blow-up.  Flooring |sin(theta)|
// regularizes the chart without perturbing any ray not already on the axis (a
// set of measure zero)."""


def _guard_sin_pole(lines, sin_sym, indent):
    """Floor |sin(theta)| in the one place it is divided by.

    Only the division is guarded; sin(theta) as a factor elsewhere is left
    exact, so the floor cannot perturb anything off the axis.
    """
    div = re.compile(rf'/{sin_sym}\b')
    hits = [i for i, l in enumerate(lines) if div.search(l)]
    if not hits:
        return lines
    safe = f'{sin_sym}_pole_safe'
    lines = [div.sub(f'/{safe}', l) if i in hits else l
             for i, l in enumerate(lines)]
    guard = [''] + [indent + c for c in POLE_COMMENT.split('\n')] + [
        f'{indent}FP const {safe} = (fabs({sin_sym}) > FP(1e-6)) '
        f'? {sin_sym} : copysign(FP(1e-6), {sin_sym});', '']
    return lines[:hits[0]] + guard + lines[hits[0]:]


def emit_function(name, a_val, Q_val, uses_a, indent='    '):
    """Render one full christoffel_*() template function definition."""
    rhs = rhs_for(a_val, Q_val)
    repl, red = cse(to_c_symbols(rhs), optimizations='basic', order='none')

    lines = [
        'template <class FP>',
        f'inline void {name}(kerr_black_hole<FP> const& hole, FP const* '
        'const __restrict__ x, FP const* const __restrict__ v, FP* const '
        '__restrict__ ch)',
        '{',
    ]
    if uses_a:
        lines.append(f'{indent}FP a = hole.a;')
    lines.append(f'{indent}FP Q = hole.Q;')
    lines.append(f'{indent}FP rs = hole.rs;')
    lines.append('')
    lines.append(f'{indent}FP y = x[2];')
    lines.append('')
    lines += [f'{indent}FP {s} = {ccode(e, user_functions=USER_FUNCS)};'
              for s, e in repl]
    lines.append('')
    lines += [f'{indent}ch[{i}] = {ccode(e, user_functions=USER_FUNCS)};'
              for i, e in enumerate(red)]

    sin_syms = [s for s, e in repl if e == sp.sin(sp.Symbol('y'))]
    if sin_syms:
        lines = _guard_sin_pole(lines, sin_syms[0], indent)
    lines.append('}')
    return '\n'.join(lines), tally(repl, red)


DISPATCHER = """\
template <class FP>
inline void christoffel(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ ch)
{
    // a=0 collapses Kerr-Newman to the spherically symmetric, non-frame-
    // dragging case (Reissner-Nordstrom, or Schwarzschild when Q is also 0):
    // no t/phi coupling and no theta-dependence in rho2, so christoffel_static
    // is ~4x fewer ops than christoffel_general (52 vs 216). A Q=0-only
    // specialization was measured and rejected - Q never couples coordinates,
    // so dropping it alone saves nothing. hole.a is the same value for every
    // thread in a render, so this branch is warp-uniform: no GPU divergence.
    if (hole.a == FP(0))
    {
        christoffel_static(hole, x, v, ch);
    }
    else
    {
        christoffel_general(hole, x, v, ch);
    }
}"""


def emit_all():
    parts = []
    total = dict(temps=0, add=0, mul=0, div=0, trans=0, total=0)
    for name, a_val, Q_val, uses_a in CASES:
        body, t = emit_function(name, a_val, Q_val, uses_a)
        parts.append(body)
        for k in total:
            total[k] += t[k]
        print(f'// {name}: temps={t["temps"]} add={t["add"]} mul={t["mul"]} '
              f'div={t["div"]} trig={t["trans"]} total={t["total"]}')
    parts.append(DISPATCHER)
    return '\n\n\n'.join(parts), total


def bench(include_slow=False):
    pipelines = [('raw (none)', lambda e: e),
                 ('factor_terms', sp.factor_terms),
                 ('cancel', sp.cancel),
                 ('simplify  [SLOW]', sp.simplify)]
    g, delta, rho2, s2 = covariant_metric(a, Q)
    inverses = [('analytic', contravariant_metric_analytic(a, Q, delta, rho2, s2)),
                ('blockwise', contravariant_metric_blockwise(g))]
    print(f'{"inverse":10s} {"simplification":18s} {"opt":6s} '
          f'{"temps":>6s} {"add":>5s} {"mul":>5s} {"div":>4s} {"trig":>5s} '
          f'{"total":>6s} {"gen_s":>7s}')
    for iname, gi in inverses:
        for pname, simp in pipelines:
            if 'SLOW' in pname and not include_slow:
                continue
            t0 = time.time()
            rhs = to_c_symbols(geodesic_rhs(g, gi, simp))
            dt = time.time() - t0
            for opt in ('basic', None):
                repl, red = cse(rhs, optimizations=opt, order='none')
                t = tally(repl, red)
                print(f'{iname:10s} {pname:18s} {str(opt):6s} '
                      f'{t["temps"]:6d} {t["add"]:5d} {t["mul"]:5d} '
                      f'{t["div"]:4d} {t["trans"]:5d} {t["total"]:6d} '
                      f'{dt:7.1f}', flush=True)


def bench_cases():
    """Op-count comparison across the a/Q specializations - the numbers
    behind picking christoffel_general + christoffel_static and rejecting a
    standalone Q=0 ("pure Kerr") case."""
    trial_cases = [
        ('general (a,Q free)', a, Q),
        ('Q=0 only ("Kerr")', a, sp.Integer(0)),
        ('a=0 ("static": Reissner-Nordstrom/Schwarzschild)', sp.Integer(0), Q),
        ('a=0, Q=0 (Schwarzschild)', sp.Integer(0), sp.Integer(0)),
    ]
    for name, a_val, Q_val in trial_cases:
        rhs = rhs_for(a_val, Q_val)
        repl, red = cse(to_c_symbols(rhs), optimizations='basic', order='none')
        t = tally(repl, red)
        print(f'{name:52s} temps={t["temps"]:3d} total={t["total"]:4d} '
              f'(add={t["add"]} mul={t["mul"]} div={t["div"]} trig={t["trans"]})')


# ---------------------------------------------------------------------------
# Verification: each generated function, and the block currently in
# cuda_ray.h, are checked against an untouched symbolic contraction.
# ---------------------------------------------------------------------------

def _generated_region():
    src = HEADER.read_text()
    start = src.index(BEGIN_MARK) + len(BEGIN_MARK)
    end = src.index(END_MARK, start)
    return src[start:end]


def _extract_function(text, name):
    start = text.index(f'inline void {name}(')
    start = text.index('\n{', start) + 2
    depth = 1
    i = start
    while depth:
        if text[i] == '{':
            depth += 1
        elif text[i] == '}':
            depth -= 1
        i += 1
    return text[start:i - 1]


def _transliterate(body):
    """Turn one generated function body into a callable python function."""
    lines = []
    for line in body.splitlines():
        line = line.strip()
        if not line or line.startswith('//'):
            continue
        line = re.sub(r'\bFP const\b', '', line)
        line = re.sub(r'\bFP\b', '', line).replace('FP(', '(')
        line = re.sub(r'pown\(', 'pw(', line).rstrip(';')
        line = line.replace('1.0/2.0', '0.5')
        tern = re.match(r'^(.*?)=\s*(.*?)\s*\?\s*(.*?)\s*:\s*(.*)$', line)
        if tern:
            line = (f'{tern.group(1)}= ({tern.group(3)}) '
                    f'if ({tern.group(2)}) else ({tern.group(4)})')
        lines.append('    ' + line.strip())
    src = ['def _f(a, Q, rs, x, v):',
           '    hole = type("H", (), dict(a=a, Q=Q, rs=rs))',
           '    ch = [0.0]*4', *lines, '    return ch']
    ns = {'pw': lambda b, e: b**e, 'sin': math.sin, 'cos': math.cos,
          'fabs': abs, 'copysign': math.copysign}
    exec('\n'.join(src), ns)
    return ns['_f']


def verify(n=4000):
    random.seed(0)
    ok = True
    try:
        region = _generated_region()
    except ValueError:
        print('(could not find generated region in cuda_ray.h)')
        region = None

    for name, a_val, Q_val, uses_a in CASES:
        ref = rhs_for(a_val, Q_val, simp=lambda e: e)   # no simplification
        gen = rhs_for(a_val, Q_val)
        args = (a, Q, rs, r, th, *v)
        f_ref = sp.lambdify(args, ref, 'math')
        f_gen = sp.lambdify(args, gen, 'math')

        f_hdr = None
        if region is not None:
            try:
                f_hdr = _transliterate(_extract_function(region, name))
            except Exception as exc:                    # noqa: BLE001
                print(f'{name}: (could not parse cuda_ray.h: {exc})')

        a_lo, a_hi = (0.0, 0.0) if a_val == 0 else (0.0, 0.95)
        Q_lo, Q_hi = (0.0, 0.0) if Q_val == 0 else (0.0, 0.3)
        worst_gen = worst_hdr = 0.0
        for _ in range(n):
            p = (random.uniform(a_lo, a_hi), random.uniform(Q_lo, Q_hi), 1.0,
                 random.uniform(1.6, 60.0), random.uniform(0.05, math.pi - 0.05),
                 *[random.uniform(-1, 1) for _ in range(4)])
            A = f_ref(*p)
            sc = max(1.0, max(abs(z) for z in A))
            B = f_gen(*p)
            worst_gen = max(worst_gen, max(abs(u - w) for u, w in zip(A, B))/sc)
            if f_hdr:
                C = f_hdr(p[0], p[1], p[2], [0.0, p[3], p[4], 0.0], list(p[5:]))
                worst_hdr = max(worst_hdr, max(abs(u - w) for u, w in zip(A, C))/sc)
        print(f'{name}: generated vs raw = {worst_gen:.3e}'
              + (f'   cuda_ray.h vs raw = {worst_hdr:.3e}' if f_hdr else ''))
        if worst_gen > 1e-9 or (f_hdr and worst_hdr > 1e-9):
            ok = False
    return ok


def check():
    """Report whether cuda_ray.h is still exactly what this script emits.

    Exit status is 0 when they match, 1 when they have drifted, so this can be
    wired into CI.
    """
    body, t = emit_all()
    try:
        live = _generated_region()
    except ValueError:
        print(f'DRIFT: could not find {BEGIN_MARK!r} / {END_MARK!r} '
              'markers in cuda_ray.h.')
        return 1

    def norm(text):
        return [l.strip() for l in text.splitlines()
                if l.strip() and not l.strip().startswith('//')]

    want, got = norm(body), norm(live)
    print(f'cuda_ray.h : {sum(l.startswith("FP x") for l in got)} temporaries')
    print(f'generated  : {t["temps"]} temporaries, {t["total"]} ops '
          f'({t["div"]} div, {t["trans"]} trig)')
    if want == got:
        print('\nin sync: cuda_ray.h matches this generator exactly.')
        return 0
    print('\nDRIFT: cuda_ray.h no longer matches this generator.')
    for i, (w, gt) in enumerate(zip(want, got)):
        if w != gt:
            print(f'  first difference at statement {i}:\n'
                  f'    header:    {gt}\n    generated: {w}')
            break
    else:
        print(f'  length differs: header {len(got)}, generated {len(want)}')
    return 1


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--bench', action='store_true',
                    help='compare simplification pipelines (general case)')
    ap.add_argument('--bench-cases', action='store_true',
                    help='compare op counts across a/Q specializations')
    ap.add_argument('--verify', action='store_true',
                    help='numeric check against the raw symbolic form')
    ap.add_argument('--check', action='store_true',
                    help='compare against the block in cuda_ray.h')
    ap.add_argument('--all', action='store_true',
                    help='include simplify() in --bench (~10 minutes)')
    args = ap.parse_args()
    if args.bench:
        bench(include_slow=args.all)
    elif args.bench_cases:
        bench_cases()
    elif args.verify:
        raise SystemExit(0 if verify() else 1)
    elif args.check:
        raise SystemExit(check())
    else:
        body, t = emit_all()
        print(body)
        print(f'\n// total: temps={t["temps"]} add={t["add"]} mul={t["mul"]} '
              f'div={t["div"]} trig={t["trans"]} total={t["total"]}')
