#!/usr/bin/env python3
"""Generate the canonical (Hamiltonian) geodesic right-hand side of
OpenACC_cli/hamiltonian.h.

The renderer integrates the geodesic equation in second-order form,
x'' = -Gamma^i_jk v^j v^k (tools/gen_christoffel.py).  That form is what an
explicit Runge-Kutta method wants, but it is *not* a Hamiltonian system in the
(x, v) variables: the symplectic form pulled back to the tangent bundle is
g_{mu nu}(x) dx^mu ^ dv^nu, which is not the canonical one, so a "symplectic"
integrator applied to (x, v) is merely symmetric, not symplectic.

The cotangent variables fix that.  With

    H(x, p) = 1/2 g^{mu nu}(x) p_mu p_nu

the geodesic flow is the canonical Hamiltonian flow

    dx^mu / dlambda =  dH/dp_mu = g^{mu nu} p_nu
    dp_mu / dlambda = -dH/dx^mu = -1/2 (d_mu g^{alpha beta}) p_alpha p_beta

on the canonical (x, p), which is what OpenACC_cli/symplectic.h integrates.
Two structural properties come for free and are worth the change of variables
on their own:

  * H is the null constraint.  g^{mu nu} p_mu p_nu = 0 defines a photon, and H
    is conserved exactly by any symplectic method up to a bounded O(h^order)
    oscillation - it does not drift secularly the way it does under Dormand-
    Prince.
  * t and phi are cyclic, so d_0 g = d_3 g = 0 identically and the generated
    code sets dp[0] = dp[3] = 0.  Photon energy and axial angular momentum are
    then conserved to the last bit, for free, by *any* integrator run in these
    variables - not to within a tolerance, exactly.

Emitted, for each of the same two specializations gen_christoffel.py uses
(general a != 0, and static a == 0), plus dispatchers on hole.a:

  hamiltonian_rhs   dx^mu and dp_mu above, contracted with p before codegen.
  raise_momentum    v^mu = g^{mu nu} p_nu    (the dx half on its own)
  lower_velocity    p_mu = g_{mu nu} v^nu    (entry point from the renderer's
                                              (x, v) state)

The metric itself is imported from gen_christoffel rather than restated, so
both generators are derived from one definition of Kerr-Newman.  Recipe and
rationale are otherwise that file's: contract with the momentum before codegen,
analytic inverse, factor_terms (not simplify), cse(optimizations='basic'),
ccode with Pow -> pown.

Usage:
    python tools/gen_hamiltonian.py            # print the generated block
    python tools/gen_hamiltonian.py --write    # splice it into hamiltonian.h
    python tools/gen_hamiltonian.py --verify   # numeric check, incl. against
                                               # the christoffel generator
    python tools/gen_hamiltonian.py --check    # diff against hamiltonian.h

Requires sympy (root pixi environment: `pixi run python tools/...`).
"""
import argparse
import math
import random
import re
from pathlib import Path

import sympy as sp
from sympy import Rational, cse, ccode

from gen_christoffel import (covariant_metric, contravariant_metric_analytic,
                             rhs_for, tally, USER_FUNCS, a, Q, rs, r, th)

HEADER = Path(__file__).resolve().parent.parent / 'OpenACC_cli' / 'hamiltonian.h'
BEGIN_MARK = '// --- BEGIN GENERATED hamiltonian: tools/gen_hamiltonian.py ---'
END_MARK = '// --- END GENERATED hamiltonian ---'

p = sp.symbols('p0 p1 p2 p3', real=True)
v = sp.symbols('v0 v1 v2 v3', real=True)
COORD = [None, r, th, None]     # t and phi are cyclic: only r, theta appear

# sin/cos(theta) as atoms, plus the reciprocal of sin as an *independent*
# symbol so that no 1/sin(theta) can survive into the emitted code without
# passing through the axis guard.  See _pole_safe.
SY, CY, ISY = sp.symbols('sy cy isy', real=True)

CASES = [
    ('general', a, Q, True),
    ('static', sp.Integer(0), Q, False),
]


def hamiltonian_rhs(gi, simp=sp.factor_terms):
    """(dx^mu, dp_mu) for H = 1/2 g^{mu nu} p_mu p_nu, contracted with p."""
    dx = [simp(sum(gi[i, j]*p[j] for j in range(4) if gi[i, j] != 0))
          for i in range(4)]
    dp = []
    for d in range(4):
        if COORD[d] is None:                       # cyclic coordinate
            dp.append(sp.Integer(0))
            continue
        acc = 0
        for i in range(4):
            for j in range(4):
                if gi[i, j] == 0:
                    continue
                acc -= Rational(1, 2)*sp.diff(gi[i, j], COORD[d])*p[i]*p[j]
        dp.append(simp(acc))
    return dx, dp


def raise_rhs(gi, simp=sp.factor_terms):
    """v^mu = g^{mu nu} p_nu."""
    return [simp(sum(gi[i, j]*p[j] for j in range(4) if gi[i, j] != 0))
            for i in range(4)]


def lower_rhs(g, simp=sp.factor_terms):
    """p_mu = g_{mu nu} v^nu."""
    return [simp(sum(g[i, j]*v[j] for j in range(4) if g[i, j] != 0))
            for i in range(4)]


def _pole_safe(exprs):
    """Move every negative power of sin(theta) onto the guarded reciprocal.

    Boyer-Lindquist g^{phi phi} carries 1/sin^2(theta), and d_theta of it
    carries 1/sin^3, so the polar axis is a coordinate singularity of the chart
    exactly as it is for the christoffel block.  It is not a physical one:
    every such term is multiplied by p_phi, which vanishes on the axis along
    with sin(theta).  Rewriting 1/sin^n as isy^n lets the emitted code floor
    |sin| in the one place it is inverted, and nowhere else, so no ray off the
    axis is perturbed.
    """
    out = []
    for e in exprs:
        e = sp.powsimp(e, deep=True)
        e = e.replace(lambda z: (z.is_Pow and z.base == SY
                                 and z.exp.is_Integer and z.exp < 0),
                      lambda z: ISY**(-z.exp))
        # Nothing may reach the printer still dividing by sin(theta): a future
        # sympy that groups the denominator differently must fail loudly here
        # rather than quietly emit a pole.
        for node in sp.preorder_traversal(e):
            if node.is_Pow and node.exp.is_number and node.exp < 0:
                if SY in node.base.free_symbols:
                    raise AssertionError(f'unguarded 1/sin(theta) in {node}')
        out.append(e)
    return out


def to_c_symbols(exprs, arrays):
    """Rename symbols to the array slots the emitted C uses.

    Before CSE, not after, for the reason gen_christoffel gives: the renaming
    changes the ordering CSE sees, so measuring the un-renamed form would not
    report the cost of the code actually emitted.
    """
    subs = {r: sp.Symbol('x[1]'), sp.sin(th): SY, sp.cos(th): CY}
    for name, syms in arrays.items():
        for i, s in enumerate(syms):
            subs[s] = sp.Symbol(f'{name}[{i}]')
    return _pole_safe([e.subs(subs) for e in exprs])


def _no_negative_pown(code):
    """pown(z, -n) -> a division, as the christoffel block already emits.

    ccode routes every integer power through pown(), including the negative
    ones, and pown() reaches them through pown_gen()'s binary-exponentiation
    loop.  That folds at -O2 for a literal exponent, but it is a shape the
    generated code has never handed the GPU compiler before, and there is no
    reason to: 1/pown(z, n) is the same arithmetic written plainly.
    """
    out, i = code, 0
    while True:
        i = out.find('pown(', i)
        if i < 0:
            return out
        depth, j = 0, i + 4
        while True:
            depth += (out[j] == '(') - (out[j] == ')')
            if depth == 0:
                break
            j += 1
        inner = out[i + 5:j]
        base, _, exp = inner.rpartition(', ')
        if exp.startswith('-') and exp[1:].isdigit():
            n = int(exp[1:])
            rep = (f'(FP(1)/({base}))' if n == 1
                   else f'(FP(1)/pown({base}, {n}))')
            out = out[:i] + rep + out[j + 1:]
        else:
            i = j
    return out


POLE_COMMENT = """\
// g^{phi phi} and its theta-derivative carry 1/sin^2 and 1/sin^3(theta).  Like
// the 1/sin in christoffel(), this is a coordinate artifact of the Boyer-
// Lindquist chart at the polar axis, not a curvature singularity: every term
// carrying it is proportional to p[3] = p_phi, which vanishes on the axis
// along with sin(theta), so the true ratio stays finite.  Flooring |sin| in
// the reciprocal - and only there, sin as a factor is left exact - regularizes
// the chart without perturbing any ray not already on the axis."""


def emit_function(kind, name, exprs, in_arrays, out_arrays, uses_a,
                  needs_pole, indent='    '):
    """Render one template function whose body is a CSE'd straight-line block.

    out_arrays is a list of (array name, [expressions]) written in order; the
    expression list passed in is their concatenation, so that CSE sees all
    outputs of one function at once and shares across them.
    """
    repl, red = cse(exprs, optimizations='basic', order='none')
    needs_pole = needs_pole and (ISY in set().union(
        *[e.free_symbols for _, e in repl], *[e.free_symbols for e in red]))
    args = ', '.join(f'FP const* const __restrict__ {n}' for n in in_arrays)
    outs = ', '.join(f'FP* const __restrict__ {n}' for n, _ in out_arrays)
    lines = [
        'template <class FP>',
        f'inline void {kind}_{name}(kerr_black_hole<FP> const& hole, {args}, '
        f'{outs})',
        '{',
    ]
    if uses_a:
        lines.append(f'{indent}FP const a = hole.a;')
    lines.append(f'{indent}FP const Q = hole.Q;')
    lines.append(f'{indent}FP const rs = hole.rs;')
    lines.append('')
    # Only the trig actually referenced: the a=0 specializations drop the
    # theta dependence of rho^2 entirely, and an unused local is a -Wall
    # warning in a header every translation unit includes.
    used = set().union(*[e.free_symbols for _, e in repl],
                       *[e.free_symbols for e in red])
    if SY in used or needs_pole:
        lines.append(f'{indent}FP const sy = sin(x[2]);')
    if CY in used:
        lines.append(f'{indent}FP const cy = cos(x[2]);')
    if needs_pole:
        lines += ['']
        lines += [indent + c for c in POLE_COMMENT.split('\n')]
        lines.append(f'{indent}FP const sy_pole_safe = (fabs(sy) > FP(1e-6)) '
                     f'? sy : copysign(FP(1e-6), sy);')
        lines.append(f'{indent}FP const isy = FP(1) / sy_pole_safe;')
    lines.append('')
    def c(e):
        return _no_negative_pown(ccode(e, user_functions=USER_FUNCS))

    lines += [f'{indent}FP const {s} = {c(e)};' for s, e in repl]
    lines.append('')
    k = 0
    for arr, items in out_arrays:
        for i in range(len(items)):
            lines.append(f'{indent}{arr}[{i}] = {c(red[k])};')
            k += 1
    lines.append('}')
    return '\n'.join(lines), tally(repl, red)


DISPATCH_COMMENT = """\
    // a=0 collapses Kerr-Newman to the spherically symmetric, non-frame-
    // dragging case, which is a much shorter block - the same split, and the
    // same warp-uniform branch, as christoffel() in cuda_ray.h."""


def dispatcher(kind, in_arrays, out_arrays):
    args = ', '.join(f'FP const* const __restrict__ {n}' for n in in_arrays)
    outs = ', '.join(f'FP* const __restrict__ {n}' for n in out_arrays)
    call = ', '.join(['hole', *in_arrays, *out_arrays])
    return '\n'.join([
        'template <class FP>',
        f'inline void {kind}(kerr_black_hole<FP> const& hole, {args}, {outs})',
        '{',
        DISPATCH_COMMENT,
        '    if (hole.a == FP(0))',
        '    {',
        f'        {kind}_static({call});',
        '    }',
        '    else',
        '    {',
        f'        {kind}_general({call});',
        '    }',
        '}',
    ])


# (kind, inputs, output array names, builder -> list of (array, exprs))
def _build(kind, a_val, Q_val):
    g, delta, rho2, s2 = covariant_metric(a_val, Q_val)
    gi = contravariant_metric_analytic(a_val, Q_val, delta, rho2, s2)
    if kind == 'hamiltonian_rhs':
        dx, dp = hamiltonian_rhs(gi)
        return ['x', 'p'], [('dx', dx), ('dp', dp)], {'p': p}
    if kind == 'raise_momentum':
        return ['x', 'p'], [('v', raise_rhs(gi))], {'p': p}
    if kind == 'lower_velocity':
        return ['x', 'v'], [('p', lower_rhs(g))], {'v': v}
    raise ValueError(kind)


KINDS = ['hamiltonian_rhs', 'raise_momentum', 'lower_velocity']


def emit_all(quiet=False):
    parts = []
    total = dict(temps=0, add=0, mul=0, div=0, trans=0, total=0)
    for kind in KINDS:
        for name, a_val, Q_val, uses_a in CASES:
            in_arrays, out_arrays, arrays = _build(kind, a_val, Q_val)
            flat = to_c_symbols([e for _, items in out_arrays for e in items],
                                arrays)
            needs_pole = any(ISY in e.free_symbols for e in flat)
            body, t = emit_function(kind, name, flat, in_arrays, out_arrays,
                                    uses_a, needs_pole)
            parts.append(body)
            for k in total:
                total[k] += t[k]
            if not quiet:
                print(f'// {kind}_{name}: temps={t["temps"]} add={t["add"]} '
                      f'mul={t["mul"]} div={t["div"]} trig={t["trans"]} '
                      f'total={t["total"]}')
        in_arrays, out_arrays, _ = _build(kind, a, Q)
        parts.append(dispatcher(kind, in_arrays, [n for n, _ in out_arrays]))
    return '\n\n\n'.join(parts), total


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def _generated_region():
    src = HEADER.read_text()
    start = src.index(BEGIN_MARK) + len(BEGIN_MARK)
    end = src.index(END_MARK, start)
    return src[start:end]


def _extract_function(text, name):
    start = text.index(f'inline void {name}(')
    start = text.index('\n{', start) + 2
    depth, i = 1, start
    while depth:
        depth += (text[i] == '{') - (text[i] == '}')
        i += 1
    return text[start:i - 1]


def _transliterate(body, in_arrays, out_arrays):
    """Turn one emitted function body into a callable python function.

    The point is to test the text that is actually in hamiltonian.h, not the
    symbolic expressions it was generated from.
    """
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
    src = [f'def _f(a, Q, rs, {", ".join(in_arrays)}):',
           '    hole = type("H", (), dict(a=a, Q=Q, rs=rs))',
           *[f'    {o} = [0.0]*4' for o in out_arrays],
           *lines,
           f'    return {", ".join(out_arrays)}']
    ns = {'pw': lambda b, e: b**e, 'sin': math.sin, 'cos': math.cos,
          'fabs': abs, 'copysign': math.copysign}
    exec('\n'.join(src), ns)
    return ns['_f']


def _sample(a_zero, Q_zero, rng):
    """A random state, and a momentum drawn to be roughly null."""
    a_val = 0.0 if a_zero else rng.uniform(0.0, 0.95)
    Q_val = 0.0 if Q_zero else rng.uniform(0.0, 0.3)
    return (a_val, Q_val, 1.0, rng.uniform(1.6, 60.0),
            rng.uniform(0.05, math.pi - 0.05),
            *[rng.uniform(-1, 1) for _ in range(4)])


def verify(n=3000):
    """Three independent checks, per specialization.

    1. the emitted expressions against the raw, unsimplified symbolic form;
    2. the block currently in hamiltonian.h against the same;
    3. the Hamiltonian system against the *christoffel* generator, i.e. that
       d/dlambda (g_{mu nu} v^nu) computed from x'' = -Gamma(v, v) agrees with
       -1/2 d_mu g^{alpha beta} p_alpha p_beta.  This is the check that the
       change of variables is right, and it ties the two generators together.
    """
    rng = random.Random(0)
    ok = True
    try:
        region = _generated_region()
    except (ValueError, FileNotFoundError):
        print('(no generated region in hamiltonian.h yet)')
        region = None

    args = (a, Q, rs, r, th)
    for name, a_val, Q_val, uses_a in CASES:
        g, delta, rho2, s2 = covariant_metric(a_val, Q_val)
        gi = contravariant_metric_analytic(a_val, Q_val, delta, rho2, s2)

        raw_dx, raw_dp = hamiltonian_rhs(gi, simp=lambda e: e)
        gen_dx, gen_dp = hamiltonian_rhs(gi)
        f_raw = sp.lambdify((*args, *p), raw_dx + raw_dp, 'math')
        f_gen = sp.lambdify((*args, *p), gen_dx + gen_dp, 'math')
        f_low = sp.lambdify((*args, *v), lower_rhs(g), 'math')

        # The same dp, but routed through the second-order form instead:
        #   p_mu' = g_{mu nu} a^nu + d_alpha g_{mu nu} v^alpha v^nu
        # with a = -Gamma(v, v) straight out of gen_christoffel.
        acc = rhs_for(a_val, Q_val, simp=lambda e: e)
        cross = []
        for mu in range(4):
            e = sum(g[mu, nu]*acc[nu] for nu in range(4))
            for alpha in (1, 2):
                for nu in range(4):
                    if g[mu, nu] != 0:
                        e += sp.diff(g[mu, nu], COORD[alpha])*v[alpha]*v[nu]
            cross.append(e)
        f_cross = sp.lambdify((*args, *v), cross, 'math')
        f_raise = sp.lambdify((*args, *p), raise_rhs(gi, simp=lambda e: e),
                              'math')

        hdr = {}
        if region is not None:
            spec = {'hamiltonian_rhs': (['x', 'p'], ['dx', 'dp']),
                    'raise_momentum': (['x', 'p'], ['v']),
                    'lower_velocity': (['x', 'v'], ['p'])}
            for kind, (ins, outs) in spec.items():
                try:
                    hdr[kind] = _transliterate(
                        _extract_function(region, f'{kind}_{name}'), ins, outs)
                except Exception as exc:                    # noqa: BLE001
                    print(f'{kind}_{name}: (could not parse header: {exc})')

        worst = dict(generated=0.0, header=0.0, christoffel=0.0, roundtrip=0.0)
        for _ in range(n):
            s = _sample(a_val == 0, Q_val == 0, rng)
            av, Qv, rsv, rv, thv = s[:5]
            vel = list(s[5:])
            mom = f_low(av, Qv, rsv, rv, thv, *vel)
            base = (av, Qv, rsv, rv, thv)

            A = f_raw(*base, *mom)
            sc = max(1.0, max(abs(z) for z in A))
            B = f_gen(*base, *mom)
            worst['generated'] = max(worst['generated'],
                                     max(abs(u - w) for u, w in zip(A, B))/sc)

            C = f_cross(*base, *vel)
            worst['christoffel'] = max(
                worst['christoffel'],
                max(abs(A[4 + i] - C[i]) for i in range(4))/sc)

            back = f_raise(*base, *mom)
            worst['roundtrip'] = max(worst['roundtrip'],
                                     max(abs(u - w) for u, w in zip(vel, back))
                                     / max(1.0, max(abs(z) for z in vel)))

            if hdr:
                x_arr = [0.0, rv, thv, 0.0]
                hp = hdr['lower_velocity'](av, Qv, rsv, x_arr, vel)
                hdx, hdp = hdr['hamiltonian_rhs'](av, Qv, rsv, x_arr, hp)
                hv = hdr['raise_momentum'](av, Qv, rsv, x_arr, hp)
                d = max(max(abs(u - w) for u, w in zip(A[:4], hdx)),
                        max(abs(u - w) for u, w in zip(A[4:], hdp)),
                        max(abs(u - w) for u, w in zip(A[:4], hv)),
                        max(abs(u - w) for u, w in zip(mom, hp)))
                worst['header'] = max(worst['header'], d/sc)

        print(f'{name:8s} generated vs raw = {worst["generated"]:.3e}   '
              f'vs christoffel = {worst["christoffel"]:.3e}   '
              f'raise(lower(v)) = {worst["roundtrip"]:.3e}'
              + (f'   header vs raw = {worst["header"]:.3e}' if hdr else ''))
        if max(worst.values()) > 1e-9:
            ok = False
    return ok


def _norm(text):
    return [l.strip() for l in text.splitlines()
            if l.strip() and not l.strip().startswith('//')]


def check():
    body, t = emit_all(quiet=True)
    try:
        live = _generated_region()
    except (ValueError, FileNotFoundError):
        print(f'DRIFT: no {BEGIN_MARK!r} / {END_MARK!r} markers in '
              f'{HEADER.name}.')
        return 1
    print(f'generated: {t["temps"]} temporaries, {t["total"]} ops '
          f'({t["div"]} div, {t["trans"]} trig)')
    if _norm(body) == _norm(live):
        print(f'in sync: {HEADER.name} matches this generator exactly.')
        return 0
    print(f'\nDRIFT: {HEADER.name} no longer matches this generator.')
    for i, (w, gt) in enumerate(zip(_norm(body), _norm(live))):
        if w != gt:
            print(f'  first difference at statement {i}:\n'
                  f'    header:    {gt}\n    generated: {w}')
            break
    else:
        print(f'  length differs: header {len(_norm(live))}, '
              f'generated {len(_norm(body))}')
    return 1


def write():
    body, _ = emit_all(quiet=True)
    src = HEADER.read_text()
    start = src.index(BEGIN_MARK) + len(BEGIN_MARK)
    end = src.index(END_MARK, start)
    HEADER.write_text(src[:start] + '\n\n\n' + body + '\n\n\n' + src[end:])
    print(f'wrote generated block into {HEADER}')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--write', action='store_true',
                    help='splice the block into hamiltonian.h in place')
    ap.add_argument('--verify', action='store_true',
                    help='numeric check against the raw symbolic form')
    ap.add_argument('--check', action='store_true',
                    help='compare against the block in hamiltonian.h')
    args = ap.parse_args()
    if args.verify:
        raise SystemExit(0 if verify() else 1)
    elif args.check:
        raise SystemExit(check())
    elif args.write:
        write()
    else:
        body, t = emit_all()
        print(body)
        print(f'\n// total: temps={t["temps"]} add={t["add"]} mul={t["mul"]} '
              f'div={t["div"]} trig={t["trans"]} total={t["total"]}')
