// Work-precision comparison of the geodesic integrators.
//
//   dopri54  the shipped explicit Dormand-Prince 5(4), adaptive, on the
//            second-order form (x, v) - cuda_ray.h
//   rk4      classical explicit Runge-Kutta 4, fixed step, on the canonical
//            form (x, p).  A control: same variables as the symplectic
//            methods, same explicitness as dopri54, no symplectic structure.
//            It is here so that any difference between the other two columns
//            can be attributed to symplecticity rather than to the change of
//            variables.
//   gl2/4/6  implicit Gauss-Legendre, fixed step, symplectic - symplectic.h
//   gl4-tt   the same, under the Sundman time transformation with kappa = 2,
//            which is how a symplectic method is allowed to vary its step
//
// Reported per configuration: right-hand side evaluations (the currency both
// families are billed in), wall time, error against a converged reference, and
// the drift of the constants of motion.
//
// Build and run:  make -C .. bench

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>

#include "../black_hole.h"
#include "../cuda_ray.h"
#include "../hamiltonian.h"
#include "../symplectic.h"

namespace {

constexpr double PI = 3.14159265358979323846;

struct State
{
    double x[D];
    double p[D];
};

struct Outcome
{
    State end;
    long evals = 0;
    long steps = 0;
    double seconds = 0.0;
    double max_H = 0.0;         // drift of the null constraint / mass shell
    double max_C = 0.0;         // drift of Carter's constant
    // Running max of |dH| snapshotted at each quarter of the run.  A bounded
    // oscillation saturates across these four numbers; a secular drift keeps
    // climbing.  That difference, not the size of any one of them, is what
    // separates a symplectic method from an accurate one.
    double drift[4] = {0.0, 0.0, 0.0, 0.0};
    double r_min = 1e30;
    double phi_total = 0.0;
    bool ok = true;
};

struct Scenario
{
    char const* name;
    char const* what;
    kerr_black_hole<double> hole;
    State start;
    double lambda_end;
    double max_step;            // dopri54's de0 cap: it must not bind
    double H0;
    double C0;
};


kerr_black_hole<double> make_hole(double a, double Q)
{
    double const x[4] = {0.0, 50.0, PI / 2, 0.0};
    double const omega[3] = {0.0, 0.0, 0.0};
    return kerr_black_hole<double>(64, 64, x, omega, a, Q, 2.0,
                                   1e-10, 1.0, 1.0, 1.0, 400.0, 6.0, 20.0);
}

// Covariant metric, signature (+,-,-,-) as the generators use.  Only needed to
// build initial conditions that satisfy g v v = norm exactly; everything after
// that goes through the generated code.
void metric(kerr_black_hole<double> const& hole, double r, double theta, double g[4][4])
{
    double const s2 = std::sin(theta) * std::sin(theta);
    double const a = hole.a;
    double const sigma = r * r + a * a * std::cos(theta) * std::cos(theta);
    double const delta = r * r - hole.rs * r + a * a + hole.Q * hole.Q;
    std::memset(g, 0, sizeof(double) * 16);
    g[0][0] = delta / sigma - a * a * s2 / sigma;
    g[1][1] = -sigma / delta;
    g[2][2] = -sigma;
    g[3][3] = (delta * a * a * s2 * s2 - a * a * a * a * s2 - 2 * a * a * r * r * s2
               - r * r * r * r * s2) / sigma;
    g[0][3] = g[3][0] = (a * a * a * s2 + a * r * r * s2 - delta * a * s2) / sigma;
}

// A geodesic with g_{mu nu} v^mu v^nu = norm (0 for light, 1 for a unit-mass
// particle), by solving the quadratic for v^t.
State make_state(kerr_black_hole<double> const& hole, double r, double theta,
                 double vr, double vtheta, double vphi, double norm)
{
    double g[4][4];
    metric(hole, r, theta, g);
    double const c = g[1][1] * vr * vr + g[2][2] * vtheta * vtheta + g[3][3] * vphi * vphi - norm;
    double const b = 2 * g[0][3] * vphi;
    double const vt = (-b + std::sqrt(b * b - 4 * g[0][0] * c)) / (2 * g[0][0]);

    State s;
    s.x[0] = 0.0; s.x[1] = r; s.x[2] = theta; s.x[3] = 0.0;
    double const v[D] = {vt, vr, vtheta, vphi};
    lower_velocity(hole, s.x, v, s.p);
    return s;
}


// ---------------------------------------------------------------------------
// The methods.  Each integrates a scenario from lambda = 0 to lambda_end and
// reports the same measurements.
// ---------------------------------------------------------------------------

void observe(Scenario const& sc, State const& s, Outcome& out, double lambda)
{
    out.max_H = std::fmax(out.max_H, std::fabs(hamiltonian_value(sc.hole, s.x, s.p) - sc.H0));
    out.max_C = std::fmax(out.max_C, std::fabs(carter_constant(sc.hole, s.x, s.p) - sc.C0));
    out.r_min = std::fmin(out.r_min, s.x[1]);
    // Latched on first crossing of each quarter, not overwritten afterwards.
    for (int k = 0; k < 4; ++k)
        if (out.drift[k] == 0.0 && lambda >= sc.lambda_end * (k + 1) / 4.0)
            out.drift[k] = out.max_H;
}

// The shipped adaptive explicit method.  It carries (x, v), so the constants
// are measured by lowering back to p at each step - a diagnostic, not part of
// the integration.
Outcome run_dopri54(Scenario const& sc, double tolerance)
{
    kerr_black_hole<double> hole = sc.hole;
    hole.errormax = tolerance;
    hole.de0 = sc.max_step;

    State s = sc.start;
    double v[D];
    raise_momentum(hole, s.x, s.p, v);

    Outcome out;
    double de = 1e-3, lambda = 0.0;
    double deriv_x[D], deriv_v[D];
    bool have_deriv = false;

    auto const t0 = std::chrono::steady_clock::now();
    while (lambda < sc.lambda_end && out.steps < 20000000)
    {
        if (lambda + de > sc.lambda_end) de = sc.lambda_end - lambda;
        // The controller may shrink the step it was handed, so the affine
        // parameter has to come from what it reports it took, not from what it
        // was asked for.  Comparing two integrators at different lambda would
        // otherwise show a difference that is not an error in either.
        double taken = 0.0;
        step(hole, s.x, v, de, deriv_x, deriv_v, have_deriv, taken);
        lambda += taken;
        ++out.steps;
        lower_velocity(hole, s.x, v, s.p);
        observe(sc, s, out, lambda);
    }
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    // 7 christoffel evaluations per accepted step (FSAL), plus rejected
    // attempts, which are not counted here - so this is a lower bound on its
    // cost, i.e. the comparison is being generous to it.
    out.evals = 7 * out.steps;
    out.end = s;
    return out;
}

// Explicit RK4 on the canonical system: the non-symplectic control.
Outcome run_rk4(Scenario const& sc, double h)
{
    State s = sc.start;
    Outcome out;
    long const steps = long(sc.lambda_end / h + 0.5);

    auto const t0 = std::chrono::steady_clock::now();
    for (long k = 0; k < steps; ++k)
    {
        double kx[4][D], kp[4][D], tx[D], tp[D];
        double const node[4] = {0.0, 0.5, 0.5, 1.0};
        for (int stage = 0; stage < 4; ++stage)
        {
            for (int n = 0; n < D; ++n)
            {
                tx[n] = s.x[n] + (stage ? h * node[stage] * kx[stage - 1][n] : 0.0);
                tp[n] = s.p[n] + (stage ? h * node[stage] * kp[stage - 1][n] : 0.0);
            }
            hamiltonian_rhs(sc.hole, tx, tp, kx[stage], kp[stage]);
        }
        for (int n = 0; n < D; ++n)
        {
            s.x[n] += h / 6 * (kx[0][n] + 2 * kx[1][n] + 2 * kx[2][n] + kx[3][n]);
            s.p[n] += h / 6 * (kp[0][n] + 2 * kp[1][n] + 2 * kp[2][n] + kp[3][n]);
        }
        ++out.steps;
        out.evals += 4;
        observe(sc, s, out, (k + 1) * h);
    }
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    out.end = s;
    return out;
}

// Implicit symplectic Gauss-Legendre, optionally time-transformed.  With
// kappa != 0 the step is constant in the fictitious time tau, so the loop runs
// until the affine parameter it has accumulated reaches lambda_end; the last
// step is trimmed in tau to land on it.
template <int S>
Outcome run_gl(Scenario const& sc, double h, int kappa)
{
    State s = sc.start;
    Outcome out;
    double const tol = 1e-14;
    int const max_iter = 50;
    // Only the combination h / r_ref^kappa matters, so r_ref is a choice of
    // units for the fictitious time.  Anchoring it at the Schwarzschild photon
    // sphere makes dtau read as "the step this would be near the photon
    // sphere", which is the region that sets the accuracy.
    double const r_ref = 3.0;
    double lambda = 0.0;

    auto const t0 = std::chrono::steady_clock::now();
    while (lambda < sc.lambda_end && out.steps < 20000000)
    {
        double step_tau = h;
        double const remaining = sc.lambda_end - lambda;
        // Landing exactly on lambda_end matters: the methods are being
        // compared at one affine parameter, so a sloppy last step would show
        // up as an error in the integrator.  Without a time transformation the
        // last step is simply trimmed.  With one, the step is in tau and the
        // affine parameter it buys is only known after the fact, so the trim
        // is corrected from what the previous attempt actually advanced.
        bool const trimming = (kappa != 0)
            ? (remaining < h * std::pow(s.x[1] / r_ref, kappa))
            : (remaining < h);
        if (trimming && kappa == 0) step_tau = remaining;
        else if (trimming) step_tau = remaining / std::pow(s.x[1] / r_ref, kappa);

        State const before = s;
        for (int attempt = 0; attempt < 4; ++attempt)
        {
            double advanced = 0.0;
            int evals = 0;
            s = before;
            bool const ok = gauss_legendre_step_safe<double, S>(
                sc.hole, s.x, s.p, step_tau, tol, max_iter, kappa, r_ref, sc.H0, advanced, evals);
            out.evals += evals;
            out.ok = out.ok && ok;
            if (!ok) { lambda = sc.lambda_end; break; }
            if (trimming && kappa != 0 && advanced > 0.0
                && std::fabs(advanced - remaining) > 1e-13 * sc.lambda_end && attempt < 3)
            {
                step_tau *= remaining / advanced;
                continue;
            }
            lambda += advanced;
            break;
        }
        ++out.steps;
        observe(sc, s, out, lambda);
    }
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    out.end = s;
    return out;
}


// ---------------------------------------------------------------------------
// Reference solution and error measure.
// ---------------------------------------------------------------------------

double position_error(Scenario const& sc, State const& s, State const& ref)
{
    double worst = 0.0;
    for (int n = 0; n < D; ++n)
        worst = std::fmax(worst, std::fabs(s.x[n] - ref.x[n]) / (1.0 + std::fabs(ref.x[n])));
    (void)sc;
    return worst;
}

// Two independent high-accuracy solutions - 6th-order implicit symplectic and
// 4th-order explicit non-symplectic, at steps small enough that both are far
// converged - so the reference comes with a measured uncertainty rather than
// an assumed one.  dopri54 is deliberately not used for this: it floors its
// own tolerance at 1e-8 (see cuda_ray.h), so it cannot resolve a reference.
State reference_for(Scenario const& sc, double& uncertainty)
{
    Outcome const a = run_gl<3>(sc, sc.lambda_end / 100000.0, 0);
    Outcome const b = run_rk4(sc, sc.lambda_end / 400000.0);
    uncertainty = position_error(sc, b.end, a.end);
    return a.end;
}


void header(char const* title)
{
    printf("\n%s\n", title);
    printf("%-10s %10s %10s %9s %10s %11s %11s\n",
           "method", "control", "rhs evals", "ms", "pos err", "max |dH|", "max |dC|");
    printf("%s\n", "-------------------------------------------------------------------------------");
}

void report(char const* name, char const* control, Outcome const& out, Scenario const& sc,
            State const& ref)
{
    printf("%-10s %10s %10ld %9.3f %10.2e %11.2e %11.2e%s\n",
           name, control, out.evals, out.seconds * 1e3,
           position_error(sc, out.end, ref), out.max_H, out.max_C,
           out.ok ? "" : "   (stage iteration did not converge)");
}

// How |H| grows along the run, which is the whole point of the exercise and is
// invisible in a single end-of-run number.  A symplectic method at constant
// step conserves the Hamiltonian of a nearby modified system, so its |dH|
// oscillates within a bound set by h and then stops growing; a non-symplectic
// method of any order accumulates.  The last column is the tell: ~1 means
// bounded, ~4 means it is still tracking the elapsed parameter.
void drift_row(char const* name, char const* control, Outcome const& out)
{
    printf("%-10s %10s %10ld", name, control, out.evals);
    for (int k = 0; k < 4; ++k) printf(" %10.2e", out.drift[k]);
    // Measured over the second half of the run.  Taking it from the first
    // quarter instead would be meaningless for the two ray cases: a ray picks
    // up essentially all of its error in the one strong-field pass, so the
    // early quarters compare "before the encounter" against "after it" rather
    // than measuring any drift.
    double const growth = out.drift[1] > 0.0 ? out.drift[3] / out.drift[1] : 0.0;
    printf("   %6.2fx  %s\n", growth, growth < 1.3 ? "bounded" : "growing");
}

// Repeat a run until it has taken at least `budget` seconds, and report the
// best time seen.  A single lensed ray is well under a millisecond.
template <class Fn>
Outcome timed(Fn fn, double budget = 0.05)
{
    Outcome best = fn();
    double spent = best.seconds;
    while (spent < budget)
    {
        Outcome const again = fn();
        spent += again.seconds;
        if (again.seconds < best.seconds) best.seconds = again.seconds;
    }
    return best;
}


void run_scenario(Scenario& sc, double const* dopri_tolerances, int n_tolerances,
                  double const* steps, int n_steps, double const* taus, int n_taus,
                  bool include_rk4)
{
    double uncertainty = 0.0;
    State const ref = reference_for(sc, uncertainty);

    Outcome const probe = run_gl<3>(sc, sc.lambda_end / 20000.0, 0);
    printf("\n== %s\n   %s\n", sc.name, sc.what);
    printf("   a=%.2f Q=%.2f, lambda 0..%.0f, r_min %.3f, total dphi %.2f rad\n",
           sc.hole.a, sc.hole.Q, sc.lambda_end, probe.r_min,
           std::fabs(probe.end.x[3] - sc.start.x[3]));
    printf("   reference: gauss-legendre 6 at h=%.2e, agreeing with dopri54 "
           "@1e-14 to %.1e\n", sc.lambda_end / 400000.0, uncertainty);

    header("   explicit, adaptive, second-order form (what the renderer ships)");
    for (int i = 0; i < n_tolerances; ++i)
    {
        char label[16];
        snprintf(label, sizeof label, "%.0e", dopri_tolerances[i]);
        double const tol = dopri_tolerances[i];
        report("dopri54", label, timed([&] { return run_dopri54(sc, tol); }), sc, ref);
    }

    if (include_rk4)
    {
        header("   explicit, fixed step, canonical form (non-symplectic control)");
        for (int i = 0; i < n_steps; ++i)
        {
            char label[16];
            snprintf(label, sizeof label, "h=%.3g", steps[i]);
            double const h = steps[i];
            report("rk4", label, timed([&] { return run_rk4(sc, h); }), sc, ref);
        }
    }

    header("   implicit, fixed step, canonical form (symplectic)");
    for (int i = 0; i < n_steps; ++i)
    {
        char label[16];
        snprintf(label, sizeof label, "h=%.3g", steps[i]);
        double const h = steps[i];
        report("gl2", label, timed([&] { return run_gl<1>(sc, h, 0); }), sc, ref);
    }
    for (int i = 0; i < n_steps; ++i)
    {
        char label[16];
        snprintf(label, sizeof label, "h=%.3g", steps[i]);
        double const h = steps[i];
        report("gl4", label, timed([&] { return run_gl<2>(sc, h, 0); }), sc, ref);
    }
    for (int i = 0; i < n_steps; ++i)
    {
        char label[16];
        snprintf(label, sizeof label, "h=%.3g", steps[i]);
        double const h = steps[i];
        report("gl6", label, timed([&] { return run_gl<3>(sc, h, 0); }), sc, ref);
    }

    header("   implicit, symplectic, Sundman time transformation (r_ref = 3)");
    for (int kappa : {1, 2})
    {
        for (int i = 0; i < n_taus; ++i)
        {
            char name[16], label[16];
            snprintf(name, sizeof name, "gl4-tt%d", kappa);
            // dtau is the step near r = r_ref; a given dtau therefore buys a
            // far-field step that grows as r^kappa, so the two kappas are
            // swept over ranges that put them at comparable total cost.
            double const h = taus[i] / (kappa == 2 ? 10.0 : 1.0);
            snprintf(label, sizeof label, "dtau=%.4g", h);
            int const k = kappa;
            report(name, label, timed([&] { return run_gl<2>(sc, h, k); }), sc, ref);
        }
    }

    printf("\n   growth of the |H| error: running max at 1/4, 1/2, 3/4, 4/4 of the run\n");
    printf("%-10s %10s %10s %10s %10s %10s %10s %9s\n",
           "method", "control", "rhs evals", "1/4", "1/2", "3/4", "4/4", "2nd half");
    printf("%s\n", "-------------------------------------------------------------------------------------------");
    {
        char label[16];
        snprintf(label, sizeof label, "%.0e", dopri_tolerances[n_tolerances - 1]);
        drift_row("dopri54", label, run_dopri54(sc, dopri_tolerances[n_tolerances - 1]));
        snprintf(label, sizeof label, "h=%.3g", steps[1]);
        drift_row("rk4", label, run_rk4(sc, steps[1]));
        drift_row("gl2", label, run_gl<1>(sc, steps[1], 0));
        drift_row("gl4", label, run_gl<2>(sc, steps[1], 0));
        drift_row("gl6", label, run_gl<3>(sc, steps[1], 0));
        snprintf(label, sizeof label, "dtau=%.3g", taus[1]);
        drift_row("gl4-tt1", label, run_gl<2>(sc, taus[1], 1));
    }
}

}   // namespace


int main(int argc, char** argv)
{
    bool const quick = (argc > 1 && std::strcmp(argv[1], "--quick") == 0);

    // 1. A strongly lensed ray: the ray tracer's actual workload.  Starts at
    //    the camera radius, grazes the photon sphere, escapes.
    Scenario lensed = {
        "strongly lensed ray (Kerr a=0.9)",
        "camera at r=50 -> grazes the photon sphere -> escapes; one pixel of a render",
        make_hole(0.9, 0.0), State(), 140.0, 4.0, 0.0, 0.0,
    };
    lensed.start = make_state(lensed.hole, 50.0, PI / 2 - 0.35, -1.0, 0.0, 0.00214, 0.0);
    lensed.H0 = hamiltonian_value(lensed.hole, lensed.start.x, lensed.start.p);
    lensed.C0 = carter_constant(lensed.hole, lensed.start.x, lensed.start.p);

    // 2. The same, but tuned to the critical impact parameter to within 1e-4:
    //    the photon-ring case the zoom movies are made of.  It winds 4 times
    //    around the prograde photon sphere at r = 1.557 - deep in the strong
    //    field, where Delta is small - before escaping, and the winding makes
    //    it exponentially sensitive to integration error.
    Scenario ring = {
        "photon-ring ray (Kerr a=0.9)",
        "critical impact parameter to 1e-4: 4 windings at r_min = 1.56, then escapes",
        make_hole(0.9, 0.0), State(), 115.0, 4.0, 0.0, 0.0,
    };
    ring.start = make_state(ring.hole, 50.0, PI / 2 - 0.02, -1.0, 0.0, 0.0011545143, 0.0);
    ring.H0 = hamiltonian_value(ring.hole, ring.start.x, ring.start.p);
    ring.C0 = carter_constant(ring.hole, ring.start.x, ring.start.p);

    // 3. A bound, precessing timelike orbit: not a ray, but the regime where
    //    the long-time behaviour of the two families separates most clearly,
    //    and the same integrators drive it.
    Scenario orbit = {
        "bound timelike orbit (Kerr a=0.9)",
        "eccentric massive-particle orbit, tens of revolutions - long-time drift test",
        make_hole(0.9, 0.0), State(), 20000.0, 100.0, 0.5, 0.0,
    };
    orbit.start = make_state(orbit.hole, 20.0, PI / 2 - 0.3, 0.0, 0.0, 0.0102, 1.0);
    orbit.H0 = hamiltonian_value(orbit.hole, orbit.start.x, orbit.start.p);
    orbit.C0 = carter_constant(orbit.hole, orbit.start.x, orbit.start.p);

    // dopri54 floors its own tolerance at 1e-8, so asking for less than that
    // changes nothing - the sweep stops where the shipped controller does.
    double const tolerances[] = {1e-4, 1e-5, 1e-6, 1e-7, 1e-8};
    double const ray_steps[] = {0.4, 0.2, 0.1, 0.05};
    double const ray_taus[] = {0.2, 0.1, 0.05, 0.025};
    double const orbit_steps[] = {2.0, 1.0, 0.5};
    double const orbit_taus[] = {0.4, 0.2, 0.1};

    printf("Geodesic integrator comparison.  Kerr-Newman in Boyer-Lindquist, "
           "rs=2 (M=1), double precision.\n");

    int const n = quick ? 2 : 4;
    run_scenario(lensed, tolerances, 5, ray_steps, n, ray_taus, n, true);
    run_scenario(ring, tolerances, 5, ray_steps, n, ray_taus, n, true);
    if (!quick)
        run_scenario(orbit, tolerances, 5, orbit_steps, 3, orbit_taus, 3, true);

    return 0;
}
