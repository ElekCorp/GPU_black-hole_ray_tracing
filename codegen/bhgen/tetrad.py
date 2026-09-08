"""Orthonormal frames for an arbitrary metric.

The camera sits at one point, so the frame it carries is a single 4x4 matrix.
The pipeline therefore builds it on the HOST, once per frame, and uploads the 16
numbers - the device never touches metric-specific frame algebra at all.  That is
deliberate: the bug this replaces was a mistranscribed closed-form Carter tetrad,
and a numerically constructed frame that is checked against e^T g e = eta cannot
be mistranscribed.

:func:`gram_schmidt` works symbolically or numerically depending on what you feed
it, so the same routine generates the reference frame for the proofs and runs at
runtime in C++ (see the emitted `tetrad` helper).
"""

from __future__ import annotations

import sympy as sp

from .spec import DIM, MetricSpec


def zamo(g: sp.Matrix, signature: int) -> list[sp.Expr]:
    """Zero-angular-momentum observer: the unit normal to the x^0 = const slices.

    u_mu ~ -delta^0_mu, i.e. u^mu ~ g^{mu 0}.  Unlike a static observer (u ~ d_t)
    this stays timelike inside an ergosphere, which makes it the right generic
    default for a rotating spacetime.
    """
    ginv = g.inv()
    norm = sp.sqrt(signature * ginv[0, 0])
    return [ginv[mu, 0] / norm for mu in range(DIM)]


def gram_schmidt(g, seed, signature, order=None, simplify=None):
    """Lorentzian Gram-Schmidt: an orthonormal frame with `seed` as its time leg.

    Returns e[b][mu], the frame vectors as coordinate components, satisfying

        g_{mu nu} e[a]^mu e[b]^nu = eta_{ab},  eta = diag(s, -s, -s, -s).

    The spatial legs are grown from the coordinate directions d_1, d_2, d_3 in
    `order` (coordinate order by default, which for a (t, r, theta, phi) chart
    reproduces the familiar radial/polar/azimuthal triad).  Each candidate has the
    already-built legs projected out with the Lorentzian projector

        w = c - sum_b [ g(c, e_b) / g(e_b, e_b) ] e_b

    and is then normalised.  No assumption is made about the chart: `order` only
    has to name three directions that stay linearly independent of the seed.
    """
    simplify = simplify if simplify is not None else (lambda e: e)
    order = list(order) if order is not None else [1, 2, 3]

    def dot(u, w):
        return sum(g[i, j] * u[i] * w[j] for i in range(DIM) for j in range(DIM))

    e0 = list(seed)
    n0 = simplify(dot(e0, e0))
    e0 = [simplify(c / sp.sqrt(signature * n0)) for c in e0]

    frame = [e0]
    for idx in order:
        c = [sp.S.One if i == idx else sp.S.Zero for i in range(DIM)]
        w = list(c)
        for eb in frame:
            coef = simplify(dot(c, eb) / dot(eb, eb))
            w = [simplify(w[i] - coef * eb[i]) for i in range(DIM)]
        nw = simplify(dot(w, w))
        frame.append([simplify(c_ / sp.sqrt(-signature * nw)) for c_ in w])
    return frame


def orthonormality_residual(g: sp.Matrix, frame, eta: sp.Matrix) -> sp.Matrix:
    """e^T g e - eta.  Zero exactly iff `frame` is an orthonormal frame for g."""
    e = sp.Matrix(DIM, DIM, lambda b, mu: frame[b][mu])
    return sp.simplify(e * g * e.T - eta)


def spec_frame(spec: MetricSpec, order=None, simplify=sp.simplify):
    """The frame the pipeline will use for `spec`: its own observer, else ZAMO."""
    seed = list(spec.observer) if spec.observer is not None else zamo(spec.g, spec.signature)
    return gram_schmidt(spec.g, seed, spec.signature, order=order, simplify=simplify)
