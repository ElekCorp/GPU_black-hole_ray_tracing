"""Flat spacetime in Cartesian coordinates (t, x, y, z).

Not a black hole - it is here because it is the cheapest possible check that
nothing in the pipeline assumes a (t, r, theta, phi) chart.  Every scalar the
scene needs is written in terms of x, y, z, the geodesic RHS must come out
identically zero, and rays must travel in straight lines.
"""

import sympy as sp

from bhgen.spec import Event, EventKind, MetricSpec, Param

t, x, y, z = sp.symbols("x0 x1 x2 x3", real=True)
R_sky, r_in, r_out = sp.symbols("R_sky r_in r_out", real=True)

g = sp.diag(1, -1, -1, -1)
rad = sp.sqrt(x**2 + y**2 + z**2)

SPEC = MetricSpec(
    name="minkowski_cartesian",
    coords=[t, x, y, z],
    params=[
        Param(R_sky, 100.0, "sky sphere radius"),
        Param(r_in, 3.0, "disk inner radius"),
        Param(r_out, 12.0, "disk outer radius"),
    ],
    g=g,
    signature=1,
    observer=[1, 0, 0, 0],
    events=[
        Event("SKY", EventKind.REGION_NEG, R_sky - rad),
        # "The equatorial plane" is z = 0 here, not theta = pi/2 - which is the
        # whole point: the event is a sign condition on a scalar, so it does not
        # care what the coordinates mean.
        Event("DISK_TOP", EventKind.CROSS_NEG_TO_POS, z,
              guards=(sp.sqrt(x**2 + y**2) - r_in, r_out - sp.sqrt(x**2 + y**2))),
        Event("DISK_BOTTOM", EventKind.CROSS_POS_TO_NEG, z,
              guards=(sp.sqrt(x**2 + y**2) - r_in, r_out - sp.sqrt(x**2 + y**2))),
    ],
    chart_doc="Cartesian (t, x, y, z); signature (+,-,-,-).",
)


# ---------------------------------------------------------------- verification
VERIFY_BOX = {t: (-1e3, 1e3), x: (-50.0, 50.0), y: (-50.0, 50.0), z: (-50.0, 50.0)}
VERIFY_CONSTRAINTS = []
ROUNDOFF_BOX = {
    sp.Symbol("x[1]"): (-50.0, 50.0),
    sp.Symbol("x[2]"): (-50.0, 50.0),
    sp.Symbol("x[3]"): (-50.0, 50.0),
    sp.Symbol("v[0]"): (-4.0, 4.0),
    sp.Symbol("v[1]"): (-4.0, 4.0),
    sp.Symbol("v[2]"): (-4.0, 4.0),
    sp.Symbol("v[3]"): (-4.0, 4.0),
}
