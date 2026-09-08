"""The three symbolic proof obligations.

Each returns a residual that is *identically zero* when the obligation holds.
"Identically zero" here means SymPy reduced it to the literal 0, which is a
proof modulo trusting SymPy's simplifier - the same standing assumption every
computer-algebra derivation in this field rests on.  It is strictly stronger
than sampling: no finite set of test points can establish an identity.
"""

from __future__ import annotations

import sympy as sp

from ..spec import DIM, MetricSpec
from .. import geodesic, tetrad
from .algebra import all_zero


def check_cross_derivation(spec: MetricSpec, acc=None, v=None):
    """Christoffel-formula RHS minus Euler-Lagrange RHS.

    The two derivations share only g and dg.  One inverts the metric explicitly
    and contracts the Christoffel formula; the other forms no Christoffel symbol
    at all and LU-solves the Euler-Lagrange system.  Agreement is therefore
    evidence about the index bookkeeping and the contraction, not a tautology.
    """
    v = v or geodesic.velocity_symbols()
    acc = acc if acc is not None else geodesic.acceleration_from_christoffel(spec, v)
    # Deliberately do NOT simplify the Lagrangian side on its own: simplifying two
    # large equivalent expressions separately costs far more than simplifying
    # their difference, which is designed to collapse to 0.
    lag = geodesic.acceleration_from_lagrangian(spec, v, simplify=False)
    # Left unreduced on purpose: algebra.all_zero decides these exactly, and far
    # faster than simplify would on a Kerr-sized expression.
    return [acc[i] - lag[i] for i in range(DIM)]


def check_cse_equivalence(exprs, reps, reduced):
    """Prove that the emitted CSE'd program computes exactly `exprs`.

    Common-subexpression elimination is a program transformation applied between
    "the mathematics we proved" and "the code we ship", so it needs its own
    obligation.  Back-substituting every temporary in reverse order reconstructs
    what the generated straight-line code actually evaluates; simplifying the
    difference against the original expression to 0 proves the transformation was
    meaning-preserving.
    """
    out = []
    for original, red in zip(exprs, reduced):
        e = red
        for sym, val in reversed(reps):
            e = e.subs(sym, val)
        out.append(e - original)
    return out


def check_tetrad(spec: MetricSpec, frame):
    """e^T g e - eta: zero exactly iff the frame is orthonormal for g.

    This is the obligation whose violation caused the wrong camera frame in the
    shipped tracer.  Discharging it symbolically covers every camera position and
    every parameter value at once.
    """
    return tetrad.orthonormality_residual(spec.g, frame, spec.eta)


def is_zero(residual):
    """Decide a residual list or matrix; returns (verdict, method, first_nonzero)."""
    return all_zero(residual)
