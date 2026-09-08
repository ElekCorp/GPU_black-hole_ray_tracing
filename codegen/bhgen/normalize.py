"""Expression normalisation applied before both code emission and verification.

Both stages must see the *same* expression or the proofs are about a different
program than the one that ships, so this lives in one place and is called by
both.
"""

from __future__ import annotations

import sympy as sp

_RATIOS = {
    sp.tan: lambda u: sp.sin(u) / sp.cos(u),
    sp.cot: lambda u: sp.cos(u) / sp.sin(u),
    sp.sec: lambda u: 1 / sp.cos(u),
    sp.csc: lambda u: 1 / sp.sin(u),
}


def rewrite_trig_ratios(expr):
    """Replace tan/cot/sec/csc by explicit sin and cos quotients.

    Two independent reasons, one numerical and one for the proofs:

    * ``tan`` overflows at theta = pi/2 while ``cot`` = cos/sin is perfectly
      finite there.  Computing cot as 1/tan therefore routes a well-behaved
      quantity through an intermediate infinity, which costs precision near the
      equatorial plane - exactly where the accretion disk is.
    * A stand-in variable for ``tan`` has an unbounded range on any interval
      containing pi/2, so Gappa can prove nothing about an expression containing
      it.  As sin/cos the same quantity is a bounded quotient.
    """
    expr = sp.sympify(expr)
    for fn, repl in _RATIOS.items():
        expr = expr.replace(fn, repl)
    return expr


def normalize(expr):
    return rewrite_trig_ratios(expr)
