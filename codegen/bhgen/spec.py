"""Chart-agnostic description of a spacetime, and of the scene drawn in it.

Everything the generated C++ needs is declared here as SymPy expressions in the
coordinates and parameters.  Nothing in this module - or anywhere downstream -
assumes the chart is (t, r, theta, phi); a metric written in Cartesian
Kerr-Schild coordinates is described exactly the same way.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Sequence

import sympy as sp

DIM = 4


class EventKind(Enum):
    """How a scalar event field turns into "the ray stops here"."""

    #: Fires as soon as f(x) < 0 anywhere along the step.  Use for regions the
    #: ray falls into and never leaves: inside a horizon, outside the sky sphere.
    REGION_NEG = "region_neg"

    #: Fires when f(x) changes sign from + to - between the step endpoints.
    CROSS_POS_TO_NEG = "cross_pos_to_neg"

    #: Fires when f(x) changes sign from - to + between the step endpoints.
    CROSS_NEG_TO_POS = "cross_neg_to_pos"


@dataclass(frozen=True)
class Event:
    """A termination condition, expressed as sign data of scalar fields.

    Sign conditions on scalars are the chart-independent way to say "the ray hit
    something": they carry over unchanged to Cartesian, Kerr-Schild or horizon-
    penetrating coordinates, whereas "r < r_horizon" or "theta crossed pi/2" do
    not.

    Attributes
    ----------
    name
        Identifier; becomes an enum constant in the generated header.
    kind
        See :class:`EventKind`.
    field
        The scalar f(x, params) whose sign is tested.
    guards
        Extra scalars that must all be > 0 for the event to fire.  This is how an
        accretion disk becomes "crossed the equatorial surface WHILE inside the
        annulus" without baking an annulus into the event field itself.
    """

    name: str
    kind: EventKind
    field: sp.Expr
    guards: Sequence[sp.Expr] = ()


@dataclass(frozen=True)
class Param:
    """A runtime scalar the generated code reads out of the parameter vector."""

    symbol: sp.Symbol
    default: float
    doc: str = ""

    @property
    def name(self) -> str:
        return str(self.symbol)


@dataclass
class MetricSpec:
    """A spacetime plus the scene rendered in it.

    Attributes
    ----------
    name
        Used for file and namespace names; must be a valid C identifier.
    coords
        The four coordinate symbols, in the order the ray tracer stores them.
    params
        Every runtime scalar the generated code may reference: metric parameters
        (mass, spin, charge) and scene parameters (sky radius, disk radii) alike.
        Keeping them in one vector is what lets `events` mention both.
    g
        The 4x4 metric, symmetric, in terms of `coords` and `params`.
    signature
        +1 for the (+,-,-,-) convention, -1 for (-,+,+,+).  The generated code
        and the tetrad builder both need to know which one you used; the geodesic
        equation itself is insensitive to it.
    observer
        Optional 4-velocity of the camera as a coordinate vector.  Determines the
        camera's rest frame, hence aberration and which way "at rest" points.  If
        omitted the pipeline falls back to the zero-angular-momentum observer
        derived numerically from `g`, which exists wherever g^{00} is timelike -
        including inside an ergosphere, where a static observer does not.
    events
        Termination conditions, in priority order.
    chart_doc
        Free text describing the chart, reproduced in the generated header.
    """

    name: str
    coords: Sequence[sp.Symbol]
    params: Sequence[Param]
    g: sp.Matrix
    signature: int = 1
    observer: Sequence[sp.Expr] | None = None
    events: Sequence[Event] = ()
    chart_doc: str = ""

    def __post_init__(self) -> None:
        if len(self.coords) != DIM:
            raise ValueError(f"{self.name}: expected {DIM} coordinates, got {len(self.coords)}")
        if self.g.shape != (DIM, DIM):
            raise ValueError(f"{self.name}: metric must be {DIM}x{DIM}, got {self.g.shape}")
        if self.signature not in (1, -1):
            raise ValueError(f"{self.name}: signature must be +1 or -1")
        asym = sp.simplify(self.g - self.g.T)
        if not all(e == 0 for e in asym):
            raise ValueError(f"{self.name}: metric is not symmetric")
        if self.observer is not None and len(self.observer) != DIM:
            raise ValueError(f"{self.name}: observer must have {DIM} components")
        allowed = set(self.coords) | {p.symbol for p in self.params}
        for label, exprs in (("metric", list(self.g)),
                             ("observer", list(self.observer or ())),
                             ("events", [e.field for e in self.events]
                                        + [gexpr for e in self.events for gexpr in e.guards])):
            free = set().union(*(sp.sympify(e).free_symbols for e in exprs)) if exprs else set()
            stray = free - allowed
            if stray:
                raise ValueError(f"{self.name}: {label} references undeclared symbols {stray}")

    @property
    def eta(self) -> sp.Matrix:
        """Minkowski metric in the spec's own signature convention."""
        s = self.signature
        return sp.diag(s, -s, -s, -s)

    @property
    def param_symbols(self) -> list[sp.Symbol]:
        return [p.symbol for p in self.params]
