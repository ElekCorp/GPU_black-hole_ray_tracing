// Checks the canonical-form geodesic right-hand side (hamiltonian.h) and the
// Gauss-Legendre integrators built on it (symplectic.h) against the
// second-order form the renderer ships (cuda_ray.h).
//
// The claims under test, in order:
//
//   1. raising and lowering are inverses, so the change of variables the
//      symplectic methods need is not itself a source of error;
//   2. the Hamiltonian right-hand side is the same geodesic as christoffel():
//      dx/dlambda is v, and dp/dlambda agrees with what differentiating
//      p = g v along the second-order equation gives;
//   3. each Gauss-Legendre member converges at its advertised order, 2s.  This
//      is what catches a mistyped tableau coefficient - nothing else in the
//      code would notice;
//   4. energy and axial angular momentum are conserved to the last bit, not to
//      a tolerance, because dp[0] and dp[3] are identically zero;
//   5. the null constraint does not drift secularly under a symplectic step,
//      while it does under Dormand-Prince, on the same ray at matched cost;
//   6. the two integrators agree on where a lensed ray actually goes.
//
// Build and run:  make -C .. test

#include <cmath>
#include <cstdint>
#include <cstdio>

#include "../black_hole.h"
#include "../cuda_ray.h"
#include "../hamiltonian.h"
#include "../symplectic.h"

namespace {

int failures = 0;

void check(bool ok, char const* what, double measured, double tolerance)
{
    printf("  %-4s  %-52s  %.3e (tolerance %.1e)\n", ok ? "ok" : "FAIL", what, measured, tolerance);
    if (!ok) ++failures;
}

kerr_black_hole<double> make_hole(double a, double Q, double rs = 2.0)
{
    double const x[4] = {0.0, 60.0, M_PI / 2, 0.0};
    double const omega[3] = {0.0, 0.0, 0.0};
    return kerr_black_hole<double>(64, 64, x, omega, a, Q, rs,
                                   1e-10, 0.05, 1.0, 1.0, 200.0, 6.0, 20.0);
}

// A photon at (r, theta) aimed by (impact parameter b, out-of-plane tilt),
// with v normalised so that g_{mu nu} v^mu v^nu = 0 to roundoff.  Solving the
// null condition for v^t is what makes it exactly null rather than nearly so.
void make_photon(kerr_black_hole<double> const& hole, double r, double theta,
                 double vr, double vtheta, double vphi, double* x, double* v)
{
    x[0] = 0.0; x[1] = r; x[2] = theta; x[3] = 0.0;
    double const s2 = std::sin(theta) * std::sin(theta);
    double const sigma = r * r + hole.a * hole.a * std::cos(theta) * std::cos(theta);
    double const delta = r * r - hole.rs * r + hole.a * hole.a + hole.Q * hole.Q;
    // Signature (+,-,-,-), matching the generators.
    double const gtt = delta / sigma - hole.a * hole.a * s2 / sigma;
    double const grr = -sigma / delta;
    double const gthth = -sigma;
    double const gphph = (delta * hole.a * hole.a * s2 * s2 - hole.a * hole.a * hole.a * hole.a * s2
                          - 2 * hole.a * hole.a * r * r * s2 - r * r * r * r * s2) / sigma;
    double const gtph = (hole.a * hole.a * hole.a * s2 + hole.a * r * r * s2
                         - delta * hole.a * s2) / sigma;

    // gtt vt^2 + 2 gtph vt vphi + (spatial part) = 0
    double const c = grr * vr * vr + gthth * vtheta * vtheta + gphph * vphi * vphi;
    double const bq = 2 * gtph * vphi;
    double const disc = bq * bq - 4 * gtt * c;
    double const vt = (-bq + std::sqrt(disc)) / (2 * gtt);
    x[0] = 0.0;
    v[0] = vt; v[1] = vr; v[2] = vtheta; v[3] = vphi;
}

double max_abs_diff(double const* a, double const* b, int n)
{
    double worst = 0.0;
    for (int i = 0; i < n; ++i) worst = std::fmax(worst, std::fabs(a[i] - b[i]));
    return worst;
}


// ---------------------------------------------------------------------------
// 1 + 2: the change of variables, and that it is the same geodesic.
// ---------------------------------------------------------------------------
void test_variable_change(kerr_black_hole<double> const& hole, char const* label)
{
    double worst_roundtrip = 0.0, worst_dx = 0.0, worst_dp = 0.0;
    unsigned seed = 12345u;
    auto rnd = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return double((seed >> 8) % 100000u) / 100000.0;
    };

    for (int trial = 0; trial < 500; ++trial)
    {
        double const r = 2.2 + 40.0 * rnd();
        double const theta = 0.15 + (M_PI - 0.30) * rnd();
        double x[D], v[D], p[D], back[D];
        make_photon(hole, r, theta, 2 * rnd() - 1, 0.4 * rnd() - 0.2, 0.2 * rnd() - 0.1, x, v);

        lower_velocity(hole, x, v, p);
        raise_momentum(hole, x, p, back);
        double scale = 0.0;
        for (int n = 0; n < D; ++n) scale = std::fmax(scale, std::fabs(v[n]));
        worst_roundtrip = std::fmax(worst_roundtrip, max_abs_diff(v, back, D) / scale);

        double dx[D], dp[D], acc[D];
        hamiltonian_rhs(hole, x, p, dx, dp);
        christoffel(hole, x, v, acc);
        worst_dx = std::fmax(worst_dx, max_abs_diff(dx, v, D) / scale);

        // p_mu' = g_{mu nu} a^nu + (d/dlambda g_{mu nu}) v^nu, evaluated by
        // differencing lower_velocity along the trajectory rather than by
        // restating d_alpha g here: it tests the shipped function.
        double const eps = 1e-6;
        double dp_fd[D];
        for (int n = 0; n < D; ++n) dp_fd[n] = 0.0;
        double x_plus[D], v_plus[D], p_plus[D], x_minus[D], v_minus[D], p_minus[D];
        for (int n = 0; n < D; ++n)
        {
            x_plus[n] = x[n] + eps * v[n];
            v_plus[n] = v[n] + eps * acc[n];
            x_minus[n] = x[n] - eps * v[n];
            v_minus[n] = v[n] - eps * acc[n];
        }
        lower_velocity(hole, x_plus, v_plus, p_plus);
        lower_velocity(hole, x_minus, v_minus, p_minus);
        double p_scale = 0.0;
        for (int n = 0; n < D; ++n) p_scale = std::fmax(p_scale, std::fabs(p[n]));
        for (int n = 0; n < D; ++n)
        {
            dp_fd[n] = (p_plus[n] - p_minus[n]) / (2 * eps);
            worst_dp = std::fmax(worst_dp, std::fabs(dp_fd[n] - dp[n]) / p_scale);
        }
    }
    char what[96];
    snprintf(what, sizeof what, "%s: raise(lower(v)) == v", label);
    check(worst_roundtrip < 1e-12, what, worst_roundtrip, 1e-12);
    snprintf(what, sizeof what, "%s: dx/dlambda == v", label);
    check(worst_dx < 1e-12, what, worst_dx, 1e-12);
    snprintf(what, sizeof what, "%s: dp/dlambda == d(g v)/dlambda", label);
    check(worst_dp < 1e-7, what, worst_dp, 1e-7);   // central difference, eps=1e-6
}


// ---------------------------------------------------------------------------
// 3: observed order of convergence.
// ---------------------------------------------------------------------------
// `base_steps` has to be chosen per method: a 6th-order method over the same
// interval as a 2nd-order one is already at double roundoff, where the error
// ratio measures nothing.  The assertion below refuses to pass on a ratio
// taken between two errors that small.
template <int S>
void test_order(kerr_black_hole<double> const& hole, char const* label, int base_steps)
{
    double x0[D], v0[D], p0[D];
    make_photon(hole, 12.0, M_PI / 2 - 0.3, -0.9, 0.05, 0.028, x0, v0);
    lower_velocity(hole, x0, v0, p0);

    double const total = 4.0;
    // Reference: the same method at a step small enough that its own error is
    // far below the coarsest one being measured, so the ratio is not polluted.
    double ref_x[D], ref_p[D];
    for (int n = 0; n < D; ++n) { ref_x[n] = x0[n]; ref_p[n] = p0[n]; }
    {
        int const steps = 20480;
        int evals = 0;
        for (int k = 0; k < steps; ++k)
            gauss_legendre_step<double, S>(hole, ref_x, ref_p, total / steps, 1e-15, 60, evals);
    }

    double err[2];
    int const coarse[2] = {base_steps, 2 * base_steps};
    for (int which = 0; which < 2; ++which)
    {
        double x[D], p[D];
        for (int n = 0; n < D; ++n) { x[n] = x0[n]; p[n] = p0[n]; }
        int evals = 0;
        for (int k = 0; k < coarse[which]; ++k)
            gauss_legendre_step<double, S>(hole, x, p, total / coarse[which], 1e-15, 60, evals);
        err[which] = std::fmax(max_abs_diff(x, ref_x, D), max_abs_diff(p, ref_p, D));
    }

    double const observed = std::log2(err[0] / err[1]);
    double const expected = double(gauss_legendre<double, S>::order);
    char what[96];
    snprintf(what, sizeof what, "%s: observed order %.2f (expected %.0f)", label, observed, expected);
    bool const measurable = err[1] > 1e-13;
    if (!measurable)
        printf("        (finer error %.2e is at roundoff - ratio is meaningless)\n", err[1]);
    check(measurable && std::fabs(observed - expected) < 0.35, what,
          std::fabs(observed - expected), 0.35);
}


// ---------------------------------------------------------------------------
// 4 + 5: what the two integrators conserve along a long ray.
// ---------------------------------------------------------------------------
void test_conservation(kerr_black_hole<double> const& hole)
{
    double x0[D], v0[D], p0[D];
    // Grazes the photon sphere and winds several times before escaping - the
    // case a photon-ring render is made of, and the hardest one here.
    make_photon(hole, 30.0, M_PI / 2, -1.0, 0.004, 0.0177, x0, v0);
    lower_velocity(hole, x0, v0, p0);

    double const total = 300.0;
    double const h = 0.02;
    int const steps = int(total / h);

    // Symplectic, fixed step.
    double x[D], p[D];
    for (int n = 0; n < D; ++n) { x[n] = x0[n]; p[n] = p0[n]; }
    double worst_H = 0.0, worst_E = 0.0, worst_L = 0.0, worst_C = 0.0;
    double const C0 = carter_constant(hole, x0, p0);
    int evals = 0;
    for (int k = 0; k < steps; ++k)
    {
        gauss_legendre_step<double, 2>(hole, x, p, h, 1e-14, 60, evals);
        worst_H = std::fmax(worst_H, std::fabs(hamiltonian_value(hole, x, p)));
        worst_E = std::fmax(worst_E, std::fabs(p[0] - p0[0]));
        worst_L = std::fmax(worst_L, std::fabs(p[3] - p0[3]));
        worst_C = std::fmax(worst_C, std::fabs(carter_constant(hole, x, p) - C0));
    }
    check(worst_E == 0.0, "GL4: energy conserved bit-exactly", worst_E, 0.0);
    check(worst_L == 0.0, "GL4: angular momentum conserved bit-exactly", worst_L, 0.0);
    check(worst_H < 1e-12, "GL4: |H| (null constraint) stays at roundoff", worst_H, 1e-12);
    check(worst_C < 1e-9, "GL4: Carter constant drift", worst_C, 1e-9);

    // Dormand-Prince, at a tolerance chosen so it uses *more* right-hand side
    // evaluations than the symplectic run above, so this is not a comparison
    // at unequal cost.
    kerr_black_hole<double> loose = hole;
    loose.errormax = 1e-11;
    loose.de0 = 0.5;
    double xd[D], vd[D];
    for (int n = 0; n < D; ++n) { xd[n] = x0[n]; vd[n] = v0[n]; }
    double de = 0.01, lambda = 0.0;
    double deriv_x[D], deriv_v[D];
    bool have_deriv = false;
    double dopri_worst_H = 0.0;
    int dopri_steps = 0;
    while (lambda < total)
    {
        double const before = de;
        step(loose, xd, vd, de, deriv_x, deriv_v, have_deriv);
        lambda += before;
        ++dopri_steps;
        double pd[D];
        lower_velocity(loose, xd, vd, pd);
        dopri_worst_H = std::fmax(dopri_worst_H, std::fabs(hamiltonian_value(loose, xd, pd)));
    }
    printf("        (GL4 %d rhs evaluations over %d steps; dopri54 ~%d evaluations over %d steps)\n",
           evals, steps, 7 * dopri_steps, dopri_steps);
    check(worst_H < dopri_worst_H, "GL4 holds the null constraint tighter than dopri54",
          worst_H / dopri_worst_H, 1.0);
}


// ---------------------------------------------------------------------------
// 6: the two integrators describe the same ray.
// ---------------------------------------------------------------------------
void test_agreement(kerr_black_hole<double> const& hole)
{
    double x0[D], v0[D];
    make_photon(hole, 40.0, M_PI / 2 - 0.2, -1.0, 0.01, 0.012, x0, v0);

    double const total = 60.0;
    double xs[D], ps[D];
    for (int n = 0; n < D; ++n) xs[n] = x0[n];
    lower_velocity(hole, x0, v0, ps);
    int evals = 0;
    int const steps = 24000;
    for (int k = 0; k < steps; ++k)
        gauss_legendre_step<double, 3>(hole, xs, ps, total / steps, 1e-15, 60, evals);

    kerr_black_hole<double> tight = hole;
    tight.errormax = 1e-13;
    tight.de0 = 0.05;
    double xd[D], vd[D], deriv_x[D], deriv_v[D];
    for (int n = 0; n < D; ++n) { xd[n] = x0[n]; vd[n] = v0[n]; }
    double de = 1e-3, lambda = 0.0;
    bool have_deriv = false;
    while (lambda < total)
    {
        double h = de;
        if (lambda + h > total) h = total - lambda;   // land exactly on `total`
        de = h;
        step(tight, xd, vd, de, deriv_x, deriv_v, have_deriv);
        lambda += h;
    }

    double vs[D];
    raise_momentum(hole, xs, ps, vs);
    double worst = 0.0;
    for (int n = 0; n < D; ++n)
        worst = std::fmax(worst, std::fabs(xs[n] - xd[n]) / (1.0 + std::fabs(xd[n])));
    check(worst < 1e-8, "GL6 and dopri54 land on the same point", worst, 1e-8);
}

}   // namespace


int main()
{
    printf("== hamiltonian.h / symplectic.h\n");

    kerr_black_hole<double> const kerr = make_hole(0.9, 0.0);
    kerr_black_hole<double> const charged = make_hole(0.6, 0.25);
    kerr_black_hole<double> const schwarzschild = make_hole(0.0, 0.0);

    test_variable_change(kerr, "kerr a=0.9");
    test_variable_change(charged, "kerr-newman a=0.6 Q=0.25");
    test_variable_change(schwarzschild, "schwarzschild");

    test_order<1>(kerr, "implicit midpoint", 40);
    test_order<2>(kerr, "gauss-legendre 4", 40);
    test_order<3>(kerr, "gauss-legendre 6", 8);

    test_conservation(kerr);
    test_agreement(kerr);

    printf("%s (%d failures)\n", failures ? "FAILED" : "passed", failures);
    return failures ? 1 : 0;
}
