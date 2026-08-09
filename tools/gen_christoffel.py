#!/usr/bin/env python3
"""Generate the CSE-optimized C body of christoffel() in OpenACC_cli/cuda_ray.h.

Emits the geodesic right-hand side  ch[i] = -Gamma^i_jk v^j v^k  for the
Kerr-Newman metric in Boyer-Lindquist coordinates, as a straight-line block of
FP temporaries.  The block in cuda_ray.h was hand-pasted from a notebook and had
no committed generator; this script reproduces it, and documents which knobs
actually matter.

The recipe (run --bench for the numbers behind it):

  1. Contract with v BEFORE generating code.  Only the 4 contracted sums are
     ever needed, never the 64 individual Gamma^i_jk, so CSE sees a far smaller
     problem.
  2. Invert the metric analytically.  The (t,phi) 2x2 block has
     det = -Delta*sin^2(theta) exactly, so g^{ij} is written in closed form and
     the sin^2 cancels out of g^tt and g^tphi instead of being carried through
     CSE as an expanded polynomial.
  3. factor_terms, NOT simplify.  This is the whole difference between the live
     block and the one commented out below it in cuda_ray.h: simplify() burns
     ~10 minutes and yields cos(4*theta)/cos(6*theta) multiple-angle terms that
     make the result *worse*.  cancel/together are catastrophic (~400 temps).
  4. cse(optimizations='basic'), then ccode with Pow -> pown.

Usage:
    python tools/gen_christoffel.py             # print the C body
    python tools/gen_christoffel.py --bench     # compare pipelines
    python tools/gen_christoffel.py --verify    # numeric check vs raw form
    python tools/gen_christoffel.py --check     # diff against cuda_ray.h

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

a, Q, rs = symbols('a Q rs', real=True, positive=True)
r, th = symbols('r theta', real=True, positive=True)
v = symbols('v0 v1 v2 v3', real=True)

delta = r**2 - rs*r + a**2 + Q**2
rho2 = r**2 + a**2*cos(th)**2
s2 = sin(th)**2

# Covariant metric, signature (+,-,-,-), matching notebooks/derivation/.
g = sp.zeros(4, 4)
g[0, 0] = delta/rho2 - a**2*s2/rho2
g[1, 1] = -rho2/delta
g[2, 2] = -rho2
g[3, 3] = (delta*a**2*sin(th)**4 - a**4*s2 - 2*a**2*r**2*s2 - r**4*s2)/rho2
g[0, 3] = g[3, 0] = (a**3*s2 + a*r**2*s2 - delta*a*s2)/rho2

COORD = [None, r, th, None]   # only r and theta appear in g


def inverse_analytic():
    """g^{ij} in closed form, using det(t,phi block) = -Delta*sin^2(theta)."""
    gi = sp.zeros(4, 4)
    A = (r**2 + a**2)**2 - delta*a**2*s2
    gi[0, 0] = A/(rho2*delta)
    gi[0, 3] = gi[3, 0] = a*(rs*r - Q**2)/(rho2*delta)
    gi[3, 3] = -(delta - a**2*s2)/(rho2*delta*s2)
    gi[1, 1] = -delta/rho2
    gi[2, 2] = -1/rho2
    return gi


def inverse_blockwise():
    """Same inverse, but letting sympy carry the expanded block determinant.

    This is what the live cuda_ray.h block was generated from; kept so --bench
    can show what the closed-form determinant buys.
    """
    gi = sp.zeros(4, 4)
    det2 = g[0, 0]*g[3, 3] - g[0, 3]**2
    gi[0, 0] = g[3, 3]/det2
    gi[3, 3] = g[0, 0]/det2
    gi[0, 3] = gi[3, 0] = -g[0, 3]/det2
    gi[1, 1] = 1/g[1, 1]
    gi[2, 2] = 1/g[2, 2]
    return gi


def geodesic_rhs(gi, simp=sp.factor_terms):
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


def emit_c(rhs, indent='        '):
    """Render the contracted RHS as C, in the naming cuda_ray.h expects."""
    repl, red = cse(to_c_symbols(rhs), optimizations='basic', order='none')
    lines = [f'{indent}FP {s} = {ccode(e, user_functions=USER_FUNCS)};'
             for s, e in repl]
    lines.append('')
    lines += [f'{indent}ch[{i}] = {ccode(e, user_functions=USER_FUNCS)};'
              for i, e in enumerate(red)]
    return '\n'.join(lines), tally(repl, red)


def bench(include_slow=False):
    pipelines = [('raw (none)', lambda e: e),
                 ('factor_terms', sp.factor_terms),
                 ('cancel', sp.cancel),
                 ('simplify  [SLOW]', sp.simplify)]
    inverses = [('analytic', inverse_analytic()),
                ('blockwise', inverse_blockwise())]
    print(f'{"inverse":10s} {"simplification":18s} {"opt":6s} '
          f'{"temps":>6s} {"add":>5s} {"mul":>5s} {"div":>4s} {"trig":>5s} '
          f'{"total":>6s} {"gen_s":>7s}')
    for iname, gi in inverses:
        for pname, simp in pipelines:
            if 'SLOW' in pname and not include_slow:
                continue
            t0 = time.time()
            rhs = to_c_symbols(geodesic_rhs(gi, simp))
            dt = time.time() - t0
            for opt in ('basic', None):
                repl, red = cse(rhs, optimizations=opt, order='none')
                t = tally(repl, red)
                print(f'{iname:10s} {pname:18s} {str(opt):6s} '
                      f'{t["temps"]:6d} {t["add"]:5d} {t["mul"]:5d} '
                      f'{t["div"]:4d} {t["trans"]:5d} {t["total"]:6d} '
                      f'{dt:7.1f}', flush=True)


# ---------------------------------------------------------------------------
# Verification: the generated form, and the block currently in cuda_ray.h, are
# both checked against an untouched symbolic contraction.
# ---------------------------------------------------------------------------

def _header_body():
    src = HEADER.read_text()
    start = src.index('FP* const __restrict__ ch)\n{\n    FP a = hole.a;')
    end = src.index('/*FP a = hole.a;', start)      # the superseded block
    return src[src.index('{', start) + 1:end]


def _header_fn():
    """Transliterate the C block in cuda_ray.h into a callable python function."""
    body = []
    for line in _header_body().splitlines():
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
        body.append('    ' + line.strip())
    src = ['def header_ch(a, Q, rs, x, v):',
           '    hole = type("H", (), dict(a=a, Q=Q, rs=rs))',
           '    ch = [0.0]*4', *body, '    return ch']
    ns = {'pw': lambda b, e: b**e, 'sin': math.sin, 'cos': math.cos,
          'fabs': abs, 'copysign': math.copysign}
    exec('\n'.join(src), ns)
    return ns['header_ch']


def verify(n=4000):
    ref = geodesic_rhs(inverse_blockwise(), lambda e: e)   # no simplification
    new = geodesic_rhs(inverse_analytic(), sp.factor_terms)
    args = (a, Q, rs, r, th, *v)
    f_ref = sp.lambdify(args, ref, 'math')
    f_new = sp.lambdify(args, new, 'math')
    try:
        f_hdr = _header_fn()
    except Exception as exc:                                # noqa: BLE001
        print(f'(could not parse cuda_ray.h: {exc})')
        f_hdr = None

    random.seed(0)
    worst_new = worst_hdr = 0.0
    for _ in range(n):
        p = (random.uniform(0, 0.95), random.uniform(0, 0.3), 1.0,
             random.uniform(1.6, 60.0), random.uniform(0.05, math.pi - 0.05),
             *[random.uniform(-1, 1) for _ in range(4)])
        A = f_ref(*p)
        sc = max(1.0, max(abs(z) for z in A))
        B = f_new(*p)
        worst_new = max(worst_new, max(abs(u - w) for u, w in zip(A, B))/sc)
        if f_hdr:
            C = f_hdr(p[0], p[1], p[2], [0.0, p[3], p[4], 0.0], list(p[5:]))
            worst_hdr = max(worst_hdr, max(abs(u - w) for u, w in zip(A, C))/sc)
    print(f'generated  vs raw symbolic reference: {worst_new:.3e}')
    if f_hdr:
        print(f'cuda_ray.h vs raw symbolic reference: {worst_hdr:.3e}')
    return worst_new


def check():
    """Report whether cuda_ray.h still matches what this script would emit."""
    body, t = emit_c(geodesic_rhs(inverse_analytic()))
    live = _header_body()
    live_temps = len(re.findall(r'^\s*FP x\d+ =', live, re.M))
    print(f'cuda_ray.h : {live_temps} temporaries')
    print(f'generated  : {t["temps"]} temporaries, {t["total"]} ops '
          f'({t["div"]} div, {t["trans"]} trig)')
    print('\nThe two are not expected to be textually identical: the live block '
          'came from\nthe blockwise inverse, this script defaults to the '
          'analytic one.  Run --verify\nto confirm both still agree '
          'numerically.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--bench', action='store_true',
                    help='compare simplification pipelines')
    ap.add_argument('--verify', action='store_true',
                    help='numeric check against the raw symbolic form')
    ap.add_argument('--check', action='store_true',
                    help='compare against the block in cuda_ray.h')
    ap.add_argument('--all', action='store_true',
                    help='include simplify() in --bench (~10 minutes)')
    args = ap.parse_args()
    if args.bench:
        bench(include_slow=args.all)
    elif args.verify:
        verify()
    elif args.check:
        check()
    else:
        body, t = emit_c(geodesic_rhs(inverse_analytic()))
        print(body)
        print(f'\n// temps={t["temps"]} add={t["add"]} mul={t["mul"]} '
              f'div={t["div"]} trig={t["trans"]} total={t["total"]}')
