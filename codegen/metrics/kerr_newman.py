"""Kerr-Newman in Boyer-Lindquist coordinates.

Written exactly as the repo's own SageManifolds and sympy.diffgeom notebooks
write it, in the (+,-,-,-) convention, so the generated code can be compared
against the block currently pasted into cuda_ray.h.
"""

import sympy as sp

from bhgen.spec import Event, EventKind, MetricSpec, Param

t, r, th, ph = sp.symbols("x0 x1 x2 x3", real=True)
a, Q, rs = sp.symbols("a Q rs", real=True)
eps_h, R_sky, r_in, r_out = sp.symbols("eps_h R_sky r_in r_out", real=True)

m = rs / 2
delta = r**2 - 2 * m * r + a**2 + Q**2
rho2 = r**2 + a**2 * sp.cos(th)**2

g = sp.zeros(4, 4)
g[0, 0] = -a**2 * sp.sin(th)**2 / rho2 + delta / rho2
g[1, 1] = -rho2 / delta
g[2, 2] = -rho2
g[3, 3] = (delta * a**2 * sp.sin(th)**4 - a**4 * sp.sin(th)**2
           - 2 * a**2 * r**2 * sp.sin(th)**2 - r**4 * sp.sin(th)**2) / rho2
g[0, 3] = g[3, 0] = (a**3 * sp.sin(th)**2 + a * r**2 * sp.sin(th)**2
                     - delta * a * sp.sin(th)**2) / rho2

SPEC = MetricSpec(
    name="kerr_newman",
    coords=[t, r, th, ph],
    params=[
        Param(a, 0.0, "spin, a = J/M"),
        Param(Q, 0.0, "electric charge"),
        Param(rs, 0.05, "Schwarzschild radius, rs = 2M"),
        Param(eps_h, 1e-3, "stop this far outside the horizon, in units of Delta"),
        Param(R_sky, 1.01, "sky sphere radius; beyond it the ray has escaped"),
        Param(r_in, 0.1, "accretion disk inner radius"),
        Param(r_out, 0.5, "accretion disk outer radius"),
    ],
    g=g,
    signature=1,
    # Carter's observer: the frame the tracer has always used, so images are
    # unchanged.  Only the direction matters - gram_schmidt normalises it.
    observer=[r**2 + a**2, 0, 0, a],
    events=[
        # Delta <= 0 IS the horizon, so no separate root-finding and no special
        # case for a naked singularity: when rs^2 < 4(a^2+Q^2), Delta is positive
        # everywhere and this event simply never fires.
        Event("HORIZON", EventKind.REGION_NEG, delta - eps_h),
        Event("SKY", EventKind.REGION_NEG, R_sky - r),
        Event("DISK_TOP", EventKind.CROSS_NEG_TO_POS, th - sp.pi / 2,
              guards=(r - r_in, r_out - r)),
        Event("DISK_BOTTOM", EventKind.CROSS_POS_TO_NEG, th - sp.pi / 2,
              guards=(r - r_in, r_out - r)),
    ],
    chart_doc="Boyer-Lindquist (t, r, theta, phi); signature (+,-,-,-).",
)


# ---------------------------------------------------------------- verification
VERIFY_BOX = {
    t:  (-1e3, 1e3),
    r:  (0.02, 50.0),
    th: (0.02, sp.pi.evalf() - 0.02),
    ph: (-1e3, 1e3),
    a:  (-0.6, 0.6),
    Q:  (0.0, 0.4),
    rs: (0.01, 1.0),
}
# Delta > 0 is exactly "outside the outer horizon".
VERIFY_CONSTRAINTS = [delta - sp.Rational(1, 1000)]

ROUNDOFF_BOX = {
    sp.Symbol("x[1]"): (0.06, 50.0),
    sp.Symbol("x[2]"): (0.02, 3.12),
    sp.Symbol("p[0]"): (-0.4, 0.4),
    sp.Symbol("p[1]"): (0.0, 0.2),
    sp.Symbol("p[2]"): (0.01, 0.05),
    sp.Symbol("v[0]"): (-4.0, 4.0),
    sp.Symbol("v[1]"): (-4.0, 4.0),
    sp.Symbol("v[2]"): (-4.0, 4.0),
    sp.Symbol("v[3]"): (-4.0, 4.0),
}
