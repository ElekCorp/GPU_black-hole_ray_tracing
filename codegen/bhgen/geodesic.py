"""Metric -> Christoffel symbols -> geodesic acceleration, plus an independent
re-derivation of the same thing for cross-checking.

The two derivations share no intermediate quantity beyond the metric and its
first derivatives:

* :func:`acceleration_from_christoffel` builds Gamma^mu_{ab} from the standard
  formula, which needs the explicit inverse metric, and contracts it with v.
* :func:`acceleration_from_lagrangian` never forms Gamma or g^{-1} at all.  It
  writes down the Euler-Lagrange equations of L = g_{mu nu} x'^mu x'^nu and
  solves the resulting linear system for x'' by LU decomposition.

Proving those two symbolically equal (see bhgen.verify.crossderiv) checks the
index bookkeeping and the contraction - by far the most error-prone part - along
a genuinely different route, so it is not a circular self-test.
"""

from __future__ import annotations

import sympy as sp

from .spec import DIM, MetricSpec


def velocity_symbols() -> list[sp.Symbol]:
    """The symbols v0..v3 standing for the ray's coordinate 4-velocity."""
    return list(sp.symbols("v0 v1 v2 v3", real=True))


def metric_derivatives(spec: MetricSpec) -> list[list[list[sp.Expr]]]:
    """dg[i][j][k] = d(g_ij) / d(x^k)."""
    return [[[sp.diff(spec.g[i, j], spec.coords[k]) for k in range(DIM)]
             for j in range(DIM)] for i in range(DIM)]


def christoffel(spec: MetricSpec, simplify: bool = False) -> list[list[list[sp.Expr]]]:
    """Gamma^mu_{alpha beta} = 1/2 g^{mu k} (d_a g_{k b} + d_b g_{k a} - d_k g_{a b})."""
    ginv = spec.g.inv().applyfunc(sp.cancel)
    dg = metric_derivatives(spec)
    out = [[[sp.S.Zero] * DIM for _ in range(DIM)] for _ in range(DIM)]
    for mu in range(DIM):
        for al in range(DIM):
            for be in range(al, DIM):
                e = sp.S.Zero
                for k in range(DIM):
                    if ginv[mu, k] == 0:
                        continue
                    e += ginv[mu, k] * (dg[k][al][be] + dg[k][be][al] - dg[al][be][k])
                e = sp.Rational(1, 2) * e
                if simplify:
                    e = sp.simplify(e)
                out[mu][al][be] = out[mu][be][al] = e
    return out


def acceleration_from_christoffel(spec: MetricSpec, v: list[sp.Symbol],
                                  simplify: bool = True) -> list[sp.Expr]:
    """A^mu = -Gamma^mu_{alpha beta} v^alpha v^beta, the geodesic RHS.

    The contraction is done before simplification: simplifying the 40 independent
    Christoffel components separately is far slower than simplifying the four
    scalars they collapse into, and produces a worse common-subexpression graph.
    """
    ginv = spec.g.inv().applyfunc(sp.cancel)
    dg = metric_derivatives(spec)
    acc = []
    for mu in range(DIM):
        s = sp.S.Zero
        for al in range(DIM):
            for be in range(DIM):
                gam = sp.S.Zero
                for k in range(DIM):
                    if ginv[mu, k] == 0:
                        continue
                    gam += ginv[mu, k] * (dg[k][al][be] + dg[k][be][al] - dg[al][be][k])
                s -= sp.Rational(1, 2) * gam * v[al] * v[be]
        acc.append(sp.simplify(s) if simplify else s)
    return acc


def acceleration_from_lagrangian(spec: MetricSpec, v: list[sp.Symbol],
                                 simplify: bool = True) -> list[sp.Expr]:
    """The same acceleration, from the Euler-Lagrange equations instead.

    For L = g_{mu nu} x'^mu x'^nu the Euler-Lagrange equations read

        2 g_{a nu} x''^nu + 2 (d_b g_{a nu}) x'^b x'^nu - (d_a g_{mu nu}) x'^mu x'^nu = 0

    which is a linear system M x'' = b with M = 2g.  Solving it by LU touches
    neither the closed-form inverse metric nor the Christoffel formula.
    """
    dg = metric_derivatives(spec)
    rhs = sp.zeros(DIM, 1)
    for al in range(DIM):
        e = sp.S.Zero
        for mu in range(DIM):
            for nu in range(DIM):
                e += dg[mu][nu][al] * v[mu] * v[nu]          # + d_a g_{mu nu} v v
                e -= 2 * dg[al][nu][mu] * v[mu] * v[nu]      # - 2 d_b g_{a nu} v^b v^nu
        rhs[al] = e
    acc = (2 * spec.g).LUsolve(rhs)
    out = [acc[i] for i in range(DIM)]
    return [sp.simplify(e) for e in out] if simplify else out
