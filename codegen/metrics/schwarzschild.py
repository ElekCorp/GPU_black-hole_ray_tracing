"""Schwarzschild in Schwarzschild coordinates - the a = Q = 0 sanity case."""

import sympy as sp

from bhgen.spec import Event, EventKind, MetricSpec, Param

t, r, th, ph = sp.symbols("x0 x1 x2 x3", real=True)
rs, eps_h, R_sky, r_in, r_out = sp.symbols("rs eps_h R_sky r_in r_out", real=True)

f = 1 - rs / r
g = sp.diag(f, -1 / f, -r**2, -r**2 * sp.sin(th)**2)

SPEC = MetricSpec(
    name="schwarzschild",
    coords=[t, r, th, ph],
    params=[
        Param(rs, 0.05, "Schwarzschild radius, rs = 2M"),
        Param(eps_h, 1e-3, "stop this far outside the horizon"),
        Param(R_sky, 1.01, "sky sphere radius"),
        Param(r_in, 0.1, "disk inner radius"),
        Param(r_out, 0.5, "disk outer radius"),
    ],
    g=g,
    signature=1,
    observer=[1, 0, 0, 0],          # static observer; fine, there is no ergosphere
    events=[
        Event("HORIZON", EventKind.REGION_NEG, r - rs - eps_h),
        Event("SKY", EventKind.REGION_NEG, R_sky - r),
        Event("DISK_TOP", EventKind.CROSS_NEG_TO_POS, th - sp.pi / 2,
              guards=(r - r_in, r_out - r)),
        Event("DISK_BOTTOM", EventKind.CROSS_POS_TO_NEG, th - sp.pi / 2,
              guards=(r - r_in, r_out - r)),
    ],
    chart_doc="Schwarzschild (t, r, theta, phi); signature (+,-,-,-).",
)


# ---------------------------------------------------------------- verification
# The region the tracer actually visits, and the parameter ranges the CLI allows.
# Coordinates first: the branch-and-bound splits the widest variable, so listing
# the coordinates before the parameters keeps the search focused where the
# geometry varies.
VERIFY_BOX = {
    t:  (-1e3, 1e3),
    r:  (0.02, 50.0),
    th: (0.02, sp.pi.evalf() - 0.02),
    ph: (-1e3, 1e3),
    rs: (0.01, 1.0),
}
# Outside the horizon: everything the ray tracer does is conditioned on this.
VERIFY_CONSTRAINTS = [r - rs - sp.Rational(1, 1000)]

# Gappa needs a bounded box in every variable, velocities included.
ROUNDOFF_BOX = {
    sp.Symbol("x[1]"): (0.06, 50.0),
    sp.Symbol("x[2]"): (0.02, 3.12),
    sp.Symbol("p[0]"): (0.01, 0.05),
    sp.Symbol("v[0]"): (-4.0, 4.0),
    sp.Symbol("v[1]"): (-4.0, 4.0),
    sp.Symbol("v[2]"): (-4.0, 4.0),
    sp.Symbol("v[3]"): (-4.0, 4.0),
}
