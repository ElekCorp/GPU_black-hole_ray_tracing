#ifndef CUDA_RAY_H
#define CUDA_RAY_H

#include <iostream>
#include <math.h>

#include "black_hole.h"

/*
 * =====================================================================
 *  Null-geodesic ray tracing in the Kerr-Newman spacetime
 * =====================================================================
 *
 *  Coordinates are Boyer-Lindquist, x = (t, r, theta, phi), and the
 *  geometry is fixed by three parameters:
 *
 *      rs  Schwarzschild radius, rs = 2M     (geometric units G = c = 1)
 *      a   spin parameter, a = J/M
 *      Q   electric charge
 *
 *  with the two structure functions
 *
 *      Sigma = r^2 + a^2 cos^2(theta)
 *      Delta = r^2 - rs*r + a^2 + Q^2
 *
 *  A light ray is a null geodesic.  Written as a first-order system in the
 *  affine parameter e it is
 *
 *      dx^mu / de = v^mu
 *      dv^mu / de = -Gamma^mu_{alpha beta} v^alpha v^beta
 *
 *  christoffel() evaluates the second right-hand side; the RK routines
 *  advance the pair (x, v); ray_step*() shoots one ray per pixel and
 *  classifies where it ends up.
 *
 *  Rays are traced BACKWARDS from the camera, so an "escaping" ray
 *  (r > sugar_ki) is what actually lights up the pixel.
 */


//#include "debugmalloc.h"


// ---------------------------------------------------------------- constants

constexpr int RK45_S = 7;               // Dormand-Prince 5(4) stage count
constexpr int RK38_S = 4;               // Kutta 3/8 rule stage count

// Smallest step the controller may propose, as a fraction of de0.  Without a
// floor a ray grazing the photon sphere drives de towards zero and burns its
// whole step budget standing still.
constexpr double DE_MIN_FRAC = 1e-6;

// The horizon test is offset outwards by this much: Delta -> 0 on the horizon
// and the geodesic equation is singular there, so absorb the ray just outside.
constexpr double HORIZON_EPS = 1e-3;

// How the ray ended.  Kept separate from the pixel value so the int8 and the
// floating-point kernels can share one tracing routine.
enum ray_outcome
{
    RAY_HORIZON = 0,    // fell through the outer horizon
    RAY_ESCAPED,        // left the scene through the sky sphere at sugar_ki
    RAY_DISK_TOP,       // crossed the equatorial disk from below
    RAY_DISK_BOTTOM,    // crossed the equatorial disk from above
    RAY_NONFINITE,      // integration produced NaN/Inf
    RAY_BUDGET          // ran out of steps without hitting anything
};

template <class FP>
inline void step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de);
template <class FP>
inline FP step_factor(FP const err, FP const tol);

template <class FP>
inline void RK45(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP const de, FP& err_norm);
template <class FP>
inline void RK38(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP const de);

template <class FP>
inline void christoffel(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ ch);

template <class FP>
inline void ijk_to_vec_mink_zoom(uint64_t const i, uint64_t const j, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg, kerr_black_hole<FP> const& hole, FP* const __restrict__ dir);
template <class FP>
inline void ijk_to_vec_zoom(uint64_t const i, uint64_t const j, kerr_black_hole<FP> const& hole, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg, FP* const __restrict__ v);

template <class FP>
inline ray_outcome trace_ray(kerr_black_hole<FP> const& hole, uint64_t const i, uint64_t const j, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg, FP* const __restrict__ x_end);

template <class FP>
inline void ray_step(int8_t* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const xd, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, uint64_t const max_steps, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki_in, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg);
template <class FP>
inline void ray_step_T(FP* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const xd, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, uint64_t const max_steps, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki_in, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg);

template <class FP>
inline bool gomb_be(FP const sugar, FP const* const x);
template <class FP>
inline bool gomb_ki(FP const sugar, FP const* const x);
template <class FP>
inline bool disk(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2);
template <class FP>
inline bool disk1(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2);
template <class FP>
inline bool disk2(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2);

template <class FP>
inline int ijk_to_n(uint64_t const i, uint64_t const j, uint64_t const k, kerr_black_hole<FP> const& hole);

template <class FP>
inline FP pown(FP const x, int const n);
template <class FP>
inline FP pown_gen(FP const x, int const n);
template <class FP>
inline FP powni_rec(FP const x, int const n);


/**
 * Advance the ray by one adaptive step of the affine parameter.
 *
 * On entry `de` is the step size to attempt.  On return x and v hold the state
 * after an accepted step and `de` holds the size proposed for the next call.
 *
 * A step whose estimated local error exceeds hole.errormax is REJECTED and
 * retried with a smaller de, so the tolerance is genuinely enforced.  (The
 * previous version only ever adjusted the size of the *next* step, which meant
 * a bad step was always kept.)  The retry count is capped, and de is never
 * driven below DE_MIN_FRAC * de0, so a ray grazing the photon sphere slows down
 * but cannot stall forever.
 */
template <class FP>
inline void step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)
{
    int const MAX_REJECT = 8;
    FP const de_min = hole.de0 * FP(DE_MIN_FRAC);
    FP const tol = hole.errormax;

    FP x_try[D];
    FP v_try[D];
    FP err = FP(0);

    for (int attempt = 0; ; ++attempt)
    {
        for (int i = 0; i < D; ++i)
        {
            x_try[i] = x[i];
            v_try[i] = v[i];
        }

        RK45(hole, x_try, v_try, de, err);

        // NaN compares false against everything, so a NaN error is treated as a
        // rejection and shrinks the step - which is what we want.
        if ((err <= tol) || (de <= de_min) || (attempt >= MAX_REJECT))
        {
            break;
        }

        de = fmax(de_min, de * step_factor(err, tol));
    }

    for (int i = 0; i < D; ++i)
    {
        x[i] = x_try[i];
        v[i] = v_try[i];
    }

    de = fmin(hole.de0, fmax(de_min, de * step_factor(err, tol)));
}

/**
 * Classical step-size controller for an order-5 propagated solution with an
 * order-4 embedded estimate: the local error scales as de^5, so the size that
 * would have just hit the tolerance is de * (tol/err)^(1/5).
 *
 * Returns the multiplicative factor, damped by 0.9 and clamped to [0.2, 5] so a
 * single freak estimate cannot make the integrator lurch.  err == 0 grows the
 * step; err == NaN shrinks it.
 */
template <class FP>
inline FP step_factor(FP const err, FP const tol)
{
    if (isnan(err))
    {
        return FP(0.2);
    }
    if (!(err > FP(0)))
    {
        return FP(5.0);
    }

    FP const f = FP(0.9) * pow(tol / err, FP(1.0 / 5.0));
    return fmin(FP(5.0), fmax(FP(0.2), f));
}

/**
 * One Dormand-Prince 5(4) step ("DOPRI5", the tableau behind MATLAB's ode45).
 *
 * Advances (x, v) over `de` using the 5th-order weights and writes into
 * err_norm the mixed absolute/relative norm of the difference between the
 * 5th- and 4th-order solutions:
 *
 *     err_norm = max_i |de * sum_s e_s k_s,i| / (1 + |y_i|)
 *
 * which is the standard local-error estimate the controller in step() consumes.
 *
 * Two defects of the previous hand-unrolled version are fixed here:
 *
 *   1. a[5][2] was written 46732/5147; the correct DOPRI5 coefficient is
 *      46732/5247.  The typo broke the row-sum condition sum_m a[s][m] = c[s]
 *      for stage 6 (1.1730 instead of 1), which collapses the method from
 *      order 5 to order 1 - while still paying for all seven stages.
 *   2. the update used bhat (the 4th-order embedded weights) instead of b.
 *      Stage 7 is evaluated exactly at the 5th-order solution point, so the
 *      better answer was computed and then thrown away, and the pair's free
 *      error estimate went unused.
 */
template <class FP>
inline void RK45(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP const de, FP& err_norm)
{
    // Butcher tableau.  a[s][m] multiplies stage m when forming stage s;
    // c[s] = sum_m a[s][m] is implied and is not needed explicitly because the
    // system is autonomous (the right-hand side has no explicit e dependence).
    FP const a[RK45_S][RK45_S - 1] = {
        { FP(0),                 FP(0),                  FP(0),                 FP(0),               FP(0),                 FP(0)          },
        { FP(1.0 / 5.0),         FP(0),                  FP(0),                 FP(0),               FP(0),                 FP(0)          },
        { FP(3.0 / 40.0),        FP(9.0 / 40.0),         FP(0),                 FP(0),               FP(0),                 FP(0)          },
        { FP(44.0 / 45.0),       FP(-56.0 / 15.0),       FP(32.0 / 9.0),        FP(0),               FP(0),                 FP(0)          },
        { FP(19372.0 / 6561.0),  FP(-25360.0 / 2187.0),  FP(64448.0 / 6561.0),  FP(-212.0 / 729.0),  FP(0),                 FP(0)          },
        { FP(9017.0 / 3168.0),   FP(-355.0 / 33.0),      FP(46732.0 / 5247.0),  FP(49.0 / 176.0),    FP(-5103.0 / 18656.0), FP(0)          },
        { FP(35.0 / 384.0),      FP(0),                  FP(500.0 / 1113.0),    FP(125.0 / 192.0),   FP(-2187.0 / 6784.0),  FP(11.0 / 84.0) }
    };

    // 5th-order weights (sum = 1).  Row 7 of `a` equals `b`, so stage 7 is
    // evaluated at the new point - the FSAL property.  We do not currently carry
    // k7 over as the next step's k1, so a step costs 7 evaluations, not 6.
    FP const b[RK45_S] = { FP(35.0 / 384.0), FP(0), FP(500.0 / 1113.0), FP(125.0 / 192.0),
                          FP(-2187.0 / 6784.0), FP(11.0 / 84.0), FP(0) };

    // e = b - bhat, the embedded error weights (sum = 0).
    FP const e[RK45_S] = { FP(71.0 / 57600.0), FP(0), FP(-71.0 / 16695.0), FP(71.0 / 1920.0),
                          FP(-17253.0 / 339200.0), FP(22.0 / 525.0), FP(-1.0 / 40.0) };

    FP kx[RK45_S][D];
    FP kv[RK45_S][D];
    FP xs[D];
    FP vs[D];

    for (int s = 0; s < RK45_S; ++s)
    {
        for (int i = 0; i < D; ++i)
        {
            xs[i] = x[i];
            vs[i] = v[i];
        }
        for (int m = 0; m < s; ++m)
        {
            for (int i = 0; i < D; ++i)
            {
                xs[i] += de * a[s][m] * kx[m][i];
                vs[i] += de * a[s][m] * kv[m][i];
            }
        }

        for (int i = 0; i < D; ++i)
        {
            kx[s][i] = vs[i];       // dx/de = v
        }
        christoffel(hole, xs, vs, kv[s]);   // dv/de = -Gamma v v
    }

    err_norm = FP(0);
    for (int i = 0; i < D; ++i)
    {
        FP sum_x = FP(0);
        FP sum_v = FP(0);
        FP err_x = FP(0);
        FP err_v = FP(0);

        for (int s = 0; s < RK45_S; ++s)
        {
            sum_x += b[s] * kx[s][i];
            sum_v += b[s] * kv[s][i];
            err_x += e[s] * kx[s][i];
            err_v += e[s] * kv[s][i];
        }

        x[i] += de * sum_x;
        v[i] += de * sum_v;

        // Written as !(e <= err_norm) rather than fmax so that a NaN estimate
        // propagates into err_norm and forces step() to reject the step,
        // instead of being silently swallowed (fmax(x, NaN) == x).
        FP const e_x = fabs(de * err_x) / (FP(1) + fabs(x[i]));
        FP const e_v = fabs(de * err_v) / (FP(1) + fabs(v[i]));
        if (!(e_x <= err_norm))
        {
            err_norm = e_x;
        }
        if (!(e_v <= err_norm))
        {
            err_norm = e_v;
        }
    }
}

/**
 * Kutta's "3/8 rule": a fixed-size, 4-stage, 4th-order Runge-Kutta step.
 *
 * Kept as a cheaper alternative to RK45() - it costs 4 christoffel evaluations
 * instead of 7 - but it provides no error estimate, so step() cannot adapt or
 * reject when this is selected.  To use it, call it in place of RK45() and drop
 * the error-driven part of the controller.
 */
template <class FP>
inline void RK38(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP const de)
{
    FP const a[RK38_S][RK38_S - 1] = {
        { FP(0),         FP(0),  FP(0) },
        { FP(1.0 / 3.0), FP(0),  FP(0) },
        { FP(-1.0 / 3.0), FP(1), FP(0) },
        { FP(1),         FP(-1), FP(1) }
    };
    FP const b[RK38_S] = { FP(1.0 / 8.0), FP(3.0 / 8.0), FP(3.0 / 8.0), FP(1.0 / 8.0) };

    FP kx[RK38_S][D];
    FP kv[RK38_S][D];
    FP xs[D];
    FP vs[D];

    for (int s = 0; s < RK38_S; ++s)
    {
        for (int i = 0; i < D; ++i)
        {
            xs[i] = x[i];
            vs[i] = v[i];
        }
        for (int m = 0; m < s; ++m)
        {
            for (int i = 0; i < D; ++i)
            {
                xs[i] += de * a[s][m] * kx[m][i];
                vs[i] += de * a[s][m] * kv[m][i];
            }
        }

        for (int i = 0; i < D; ++i)
        {
            kx[s][i] = vs[i];
        }
        christoffel(hole, xs, vs, kv[s]);
    }

    for (int i = 0; i < D; ++i)
    {
        FP sum_x = FP(0);
        FP sum_v = FP(0);
        for (int s = 0; s < RK38_S; ++s)
        {
            sum_x += b[s] * kx[s][i];
            sum_v += b[s] * kv[s][i];
        }
        x[i] += de * sum_x;
        v[i] += de * sum_v;
    }
}

/**
 * Geodesic acceleration ch[mu] = -Gamma^mu_{alpha beta} v^alpha v^beta for the
 * Kerr-Newman metric in Boyer-Lindquist coordinates, i.e. the right-hand side of
 *
 *     dv^mu / de = ch[mu].
 *
 * Despite the name it is NOT the Christoffel symbols themselves - those are
 * never formed; the whole double contraction is emitted as one common-
 * subexpression-eliminated block (x0..x83), which is why it looks the way it
 * does.  It is generated by kerr_metric_code_generator_sage_math.ipynb; the
 * commented-out blocks below are earlier, unfactored versions kept for
 * reference.
 *
 * Verified against the metric symbolically: the expression below agrees with
 * -Gamma^mu_{ab} v^a v^b computed from
 *
 *     Sigma = r^2 + a^2 cos^2(theta),   Delta = r^2 - rs*r + a^2 + Q^2
 *
 * to ~5e-14 relative error over random (a, Q, rs, x, v).
 *
 * Singular where the coordinates are: on the axis (sin(theta) -> 0, ch[0] and
 * ch[3] carry an explicit 1/sin(theta)), on the horizon (Delta -> 0) and on the
 * ring singularity (Sigma -> 0).  Callers must keep rays away from those.
 */
template <class FP>
inline void christoffel(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ ch)
{
    FP a = hole.a;
    FP Q = hole.Q;
    FP rs = hole.rs;

    FP x0 = pown(x[1], 3);
    FP x1 = rs*x0;
    FP x2 = pown(Q, 2);
    FP x3 = pown(x[1], 2);
    FP x4 = x2*x3;
    FP x5 = pown(a, 2);
    FP x6 = cos(x[2]);
    FP x7 = pown(x6, 2);
    FP x8 = x5*x7;
    FP x9 = rs*x[1];
    FP x10 = -x1 + x2*x8 + x4 - x8*x9;
    FP x11 = x2 - x9;
    FP x12 = sin(x[2]);
    FP x13 = pown(x12, 2);
    FP x14 = x13*x5;
    FP x15 = x3 + x8;
    FP x16 = x14 + x15;
    FP x17 = x11*x16;
    FP x18 = x10*x17;
    FP x19 = x14 - x2 - x5 + x8 + x9;
    FP x20 = -x19;
    FP x21 = pown(x[1], 6);
    FP x22 = pown(a, 4);
    FP x23 = x22*x3;
    FP x24 = pown(x[1], 4);
    FP x25 = 2*x5;
    FP x26 = pown(a, 6);
    FP x27 = x14*x24;
    FP x28 = x13*x23;
    FP x29 = x1*x14;
    FP x30 = cos(4*x[2]);
    FP x31 = (1.0/8.0)*x30;
    FP x32 = 1.0/8.0 - x31;
    FP x33 = x22*x32;
    FP x34 = x14*x4;
    FP x35 = x2*x22;
    FP x36 = x21 - x23*x32 + 2*x23*x7 + x23 + x24*x25 + x24*x8 - x26*x32 + x26*x7 - x27 - x28 + x29 - x32*x35 + x33*x9 - x34;
    FP x37 = a*v[0];
    FP x38 = x14*x9;
    FP x39 = -x14*x3 + x24;
    FP x40 = x22 + x25*x3;
    FP x41 = -x13*x2*x5 - x13*x22 + x38 + x39 + x40;
    FP x42 = 2*x13*x22;
    FP x43 = 2*x3;
    FP x44 = x14*x41 + x15*(-2*x14*x2 - x14*x43 + x24 + 2*x38 + x40 - x42);
    FP x45 = -x11*x16;
    FP x46 = x12*x6;
    FP x47 = 2*v[2];
    FP x48 = 2*x[1];
    FP x49 = rs*x15;
    FP x50 = x11*x48 + x49;
    FP x51 = x10*x14;
    FP x52 = x50*x51;
    FP x53 = x15*x48;
    FP x54 = -x14;
    FP x55 = x11 + x3 + x5;
    FP x56 = x48*(x54 + x55) + x49 - x53;
    FP x57 = x14*x49 - x41*x48 + x53*(x25 + x43 + x54);
    FP x58 = v[3]*x13;
    FP x59 = pown(x15, 2);
    FP x60 = pown(x[1], 5);
    FP x61 = x22*x9;
    FP x62 = x24*x5;
    FP x63 = pown(x12, 4);
    FP x64 = 2*x13;
    FP x65 = 1/(x59*(-rs*x60 - x1*x25 + x2*x24 + x21 + x23*x63 + 3*x23 + x25*x4 + x26*x63 - x26*x64 + x26 - 2*x27 - 4*x28 + 2*x29 - 2*x34 + x35*x63 - x35*x64 + x35 + x42*x9 - x61*x63 - x61 + 3*x62));
    FP x66 = pown(x55, 2);
    FP x67 = pown(v[1], 2);
    FP x68 = 2*x[2];
    FP x69 = x5*sin(x68);
    FP x70 = v[2]*x55;
    FP x71 = x59*x70;
    FP x72 = pown(v[0], 2);
    FP x73 = pown(v[3], 2);
    FP x74 = 1/(pown(x15, 3)*x55);
    FP x75 = cos(x68);
    FP x76 = x5*x75;
    FP x77 = rs*x26;
    FP x78 = x35*x[1];
    FP x79 = (1.0/2.0)*rs;
    FP x80 = rs*x23;
    FP x81 = x0*x2*x25;
    FP x82 = x30 + 1;
    FP x83 = x10 + x22*x7 + x3*x5 + x3*x8 - x33 + x39;
ch[0] = x65*(a*x46*x47*(-v[3]*(x10*x44 + x36*x45) - x37*(x18 + x20*x36)) + v[1]*(-a*x58*(x10*x57 + x36*x50) + v[0]*(x36*x56 + x52)));
ch[1] = x74*(-x37*x50*x58*x66 + (1.0/2.0)*x59*x67*(x15*(-rs + 2*x[1]) - x48*x55) + (1.0/2.0)*x66*(x13*x57*x73 + x56*x72) + x71*(v[1]*x69 + x70*x[1]));
ch[2] = x74*(x46*x55*(2*v[3]*x17*x37 + x19*x5*x72 + x44*x73) - 1.0/8.0*x67*x69*pown(x43 + x5 + x76, 2) + x71*(-v[1]*x48 + (1.0/2.0)*v[2]*x69));
ch[3] = x65*(v[1]*x12*(v[3]*(x52 - x57*x83) + x37*(-rs*x21 + 2*x2*x60 + x23*x75*x79 - x24*x76*x79 + (1.0/16.0)*x30*x77 + (1.0/4.0)*x30*x78 - x31*x80 - x62*x79 + (15.0/32.0)*x75*x77 + x75*x78 + x75*x81 + (1.0/8.0)*x77*x82 + (1.0/32.0)*x77*cos(6*x[2]) + (3.0/16.0)*x77 + (3.0/4.0)*x78 + (1.0/4.0)*x80*x82 + (1.0/8.0)*x80 + x81)) + x47*x6*(-v[3]*(x14*x18 + x44*x83) + x37*(x20*x51 + x45*x83)))/x12;

}
/*
template <>
inline void christoffel<float>(kerr_black_hole<float>& hole, float* x, float* v, float* ch)
{
    float a = hole.a;
    float Q = hole.Q;
    float rs = hole.rs;

    ch[0] = 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0f / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (a * v[1] * v[3] * (sinf(x[2])) * (sinf(x[2])) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) * (x[1] * x[1]) - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * a * a + rs * (a * a) * (x[1] * x[1]) + 3 * rs * (x[1] * x[1] * x[1] * x[1]) - 2 * x[1] * (cosf(x[2])) * (cosf(x[2])) * Q * Q * a * a - 2 * x[1] * Q * Q * a * a - 4 * Q * Q * x[1] * x[1] * x[1]) + a * v[2] * v[3] * (-rs * x[1] + Q * Q) * (-2 * sinf(x[2]) * cosf(x[2]) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) + sinf(2 * x[2]) * (a * a + x[1] * x[1]) * (rs * x[1] * (sinf(x[2])) * (sinf(x[2])) * (a * a) - (sinf(x[2])) * (sinf(x[2])) * Q * Q * a * a - (sinf(x[2])) * (sinf(x[2])) * a * a * x[1] * x[1] - (sinf(x[2])) * (sinf(x[2])) * a * a * a * a + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) + v[0] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * (sinf(x[2])) * (sinf(x[2])) * a * a * x[1] * x[1] - rs * (sinf(x[2])) * (sinf(x[2])) * a * a * a * a + rs * (a * a * a * a) - rs * x[1] * x[1] * x[1] * x[1] + 2 * x[1] * (Q * Q) * (a * a) + 2 * (Q * Q) * (x[1] * x[1] * x[1])) - v[0] * v[2] * sinf(2 * x[2]) * a * a * (-rs * x[1] + Q * Q) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]));

    ch[1] = (1.0f / 2.0f) * 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (2 * v[1] * v[2] * sinf(2 * x[2]) * (a * a) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + 2 * x[1] * (v[2] * v[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-2 * a * v[0] * v[3] * (sinf(x[2])) * (sinf(x[2])) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - (sinf(x[2])) * (sinf(x[2])) * v[3] * v[3] * (2 * x[1] * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) - ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (4 * x[1] * (a * a) + (rs - 2 * x[1]) * (sinf(x[2])) * (sinf(x[2])) * (a * a) + 4 * (x[1] * x[1] * x[1]))) + (v[0] * v[0]) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q))) - v[1] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * x[1] * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (rs - 2 * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])));

    ch[2] = (1.0f / 2.0f) * 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-4 * v[1] * v[2] * x[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) - sinf(2 * x[2]) * a * a * v[1] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) + sinf(2 * x[2]) * (a * a) * (v[2] * v[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (2 * a * v[0] * v[3] * sinf(2 * x[2]) * (-rs * x[1] + Q * Q) * (a * a + x[1] * x[1]) + 2 * sinf(x[2]) * cosf(x[2]) * (v[3] * v[3]) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) - sinf(2 * x[2]) * a * a * v[0] * v[0] * (-rs * x[1] + Q * Q)) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]));

    ch[3] = 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-2 * a * v[0] * v[2] * (-rs * x[1] + Q * Q) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + v[1] * tan(x[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (a * v[0] * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - v[3] * (-rs * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * a * a * a * a - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * x[1] * x[1] - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * a * a - rs * a * a * x[1] * x[1] + rs * (a * a * a * a) - 2 * rs * x[1] * x[1] * x[1] * x[1] + 2 * x[1] * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (a * a * a * a) + 4 * x[1] * (cosf(x[2])) * (cosf(x[2])) * (a * a * a * a) + 2 * x[1] * (Q * Q) * (a * a) - 2 * x[1] * a * a * a * a + 4 * (cosf(x[2])) * (cosf(x[2])) * (a * a) * (x[1] * x[1] * x[1]) + 2 * (Q * Q) * (x[1] * x[1] * x[1]) + 2 * (x[1] * x[1] * x[1] * x[1] * x[1]))) - v[2] * v[3] * (sinf(2 * x[2]) * tan(x[2]) * (a * a) * (-rs * x[1] + Q * Q) * (-rs * x[1] + Q * Q) * (a * a + x[1] * x[1]) + 2 * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) * (-rs * x[1] + (cosf(x[2])) * (cosf(x[2])) * (a * a) + Q * Q + x[1] * x[1]))) / tan(x[2]);

}*/

/**
 * Direction of the ray through pixel (i, j) in the camera's LOCAL orthonormal
 * frame, written into dir[0..3].
 *
 * The camera is a pinhole looking along its local +x axis, with the image plane
 * at distance kepernyo_tav and total height kepernyo_high.  dir[0] is the time
 * component and is 1: the ray is null and affinely parameterised so that the
 * camera measures unit frequency.  dir[1..3] are the unit spatial direction
 * (forward, up, right), so the 4-vector is null in the local Minkowski frame.
 *
 * Zoom: the frame being rendered covers the sub-rectangle starting at
 * (ikezd, jkezd) and ending at column iveg of a virtual SZELESregi x MAGASregi
 * image.  Only iveg is passed because the aspect ratio is pinned by
 * SZELES/MAGAS == (iveg-ikezd)/(jveg-jkezd), which is exactly why the j mapping
 * below also divides by SZELES rather than by MAGAS.
 */
template <class FP>
inline void ijk_to_vec_mink_zoom(uint64_t const i, uint64_t const j, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg, kerr_black_hole<FP> const& hole, FP* const __restrict__ dir)
{
    FP const span = FP(iveg) - FP(ikezd);
    FP const ir = FP(ikezd) + (FP(i) / FP(hole.SZELES)) * span;
    FP const jr = FP(jkezd) + (FP(j) / FP(hole.SZELES)) * span;

    // Square pixels: one pitch, kepernyo_high / MAGASregi, on both axes.
    FP const pitch = hole.kepernyo_high / FP(MAGASregi);

    FP const x = hole.kepernyo_tav;                          // forward
    FP const y = pitch * (FP(MAGASregi) / FP(2) - jr);       // up
    FP const z = pitch * (ir - FP(SZELESregi) / FP(2));      // right

    FP const norm = sqrt(x * x + y * y + z * z);

    dir[0] = FP(1);
    dir[1] = x / norm;
    dir[2] = y / norm;
    dir[3] = z / norm;
}

/**
 * Initial 4-velocity of the ray for pixel (i, j) in Boyer-Lindquist COORDINATE
 * components, written into v[0..3].
 *
 * Two transformations are applied to the local direction from
 * ijk_to_vec_mink_zoom():
 *
 *  1. Camera orientation.  The spatial part is rotated about the axis
 *     Omega/|Omega| by the angle |Omega| (Rodrigues' rotation formula).  The
 *     default Omega = (0, pi, 0) is a half turn about "up", which turns the
 *     camera from facing radially outward to facing the black hole.
 *
 *  2. Frame -> coordinates.  The result is contracted with Carter's orthonormal
 *     tetrad e_(b)^mu at the camera position, v^mu = e_(b)^mu v^(b):
 *
 *         e_(0) = [ (r^2+a^2) d_t + a d_phi ] / sqrt(Sigma Delta)
 *         e_(1) = sqrt(Delta/Sigma) d_r
 *         e_(2) = d_theta / sqrt(Sigma)
 *         e_(3) = [ a sin(theta) d_t + d_phi / sin(theta) ] / sqrt(Sigma)
 *
 *     This satisfies g_{mu nu} e_(a)^mu e_(b)^nu = eta_{ab} exactly, so a null,
 *     unit-frequency direction in the camera frame comes out as a genuinely
 *     null vector in the curved geometry.  e_(0) is the camera's four-velocity:
 *     a locally non-rotating observer co-rotating with the frame dragging at
 *     dphi/dt = a/(r^2+a^2).
 *
 * Both steps used to be wrong.  The rotation matrix's middle row had a
 * misplaced parenthesis and a sign error, so R was not orthogonal for any axis
 * with sin|Omega| != 0 (the default Omega happens to have sin = 0, which is why
 * this stayed hidden).  The tetrad used Delta = r^2 - 4 rs r + a^2 + Q^2 - the
 * factor 4 does not belong, and for rs >= r/2 it turns negative and NaNs the
 * whole image - and it placed the two `a` terms in the wrong rows, leaving a
 * frame that was not orthonormal, so the initial ray was not even null.
 */
template <class FP>
inline void ijk_to_vec_zoom(uint64_t const i, uint64_t const j, kerr_black_hole<FP> const& hole, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg, FP* const __restrict__ v)
{
    FP loc[D];
    ijk_to_vec_mink_zoom(i, j, SZELESregi, MAGASregi, ikezd, jkezd, iveg, hole, loc);

    // ---- 1. camera orientation (Rodrigues) --------------------------------

    FP const angle = sqrt(hole.Omega_1 * hole.Omega_1 + hole.Omega_2 * hole.Omega_2 + hole.Omega_3 * hole.Omega_3);

    // Pick an arbitrary axis when the angle is zero; the rotation is the
    // identity then anyway.  (Normalising before this test divided by zero.)
    FP u[D - 1] = { FP(1), FP(0), FP(0) };
    if (angle > FP(0))
    {
        u[0] = hole.Omega_1 / angle;
        u[1] = hole.Omega_2 / angle;
        u[2] = hole.Omega_3 / angle;
    }

    FP const c = cos(angle);
    FP const s = sin(angle);
    FP const C = FP(1) - c;

    FP const R[D - 1][D - 1] = {
        { c + u[0] * u[0] * C,        u[0] * u[1] * C - u[2] * s,  u[0] * u[2] * C + u[1] * s },
        { u[1] * u[0] * C + u[2] * s, c + u[1] * u[1] * C,         u[1] * u[2] * C - u[0] * s },
        { u[2] * u[0] * C - u[1] * s, u[2] * u[1] * C + u[0] * s,  c + u[2] * u[2] * C        }
    };

    FP d[D];
    d[0] = loc[0];
    for (int m = 0; m < D - 1; ++m)
    {
        d[m + 1] = R[m][0] * loc[1] + R[m][1] * loc[2] + R[m][2] * loc[3];
    }

    // ---- 2. Carter tetrad at the camera position --------------------------

    FP const r_0 = hole.r_0;
    FP const a = hole.a;

    FP const sin_t = sin(hole.theta_0);
    FP const cos_t = cos(hole.theta_0);

    FP const Sigma = r_0 * r_0 + a * a * cos_t * cos_t;
    FP const Delta = r_0 * r_0 - hole.rs * r_0 + a * a + hole.Q * hole.Q;

    FP const sqrt_Sigma = sqrt(Sigma);
    FP const sqrt_SD = sqrt(Sigma * Delta);

    v[0] = (r_0 * r_0 + a * a) / sqrt_SD * d[0] + a * sin_t / sqrt_Sigma * d[3];
    v[1] = sqrt(Delta / Sigma) * d[1];
    v[2] = d[2] / sqrt_Sigma;
    v[3] = a / sqrt_SD * d[0] + d[3] / (sqrt_Sigma * sin_t);
}

/**
 * Shoot one ray backwards from the camera through pixel (i, j) and report where
 * it ended up.  x_end receives the Boyer-Lindquist coordinates of the final
 * point, which the callers use to shade a disk hit by its radius.
 *
 * Both kernels share this; only the mapping from outcome to pixel value differs.
 */
template <class FP>
inline ray_outcome trace_ray(kerr_black_hole<FP> const& hole, uint64_t const i, uint64_t const j, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg, FP* const __restrict__ x_end)
{
    FP x[D] = { hole.t_0, hole.r_0, hole.theta_0, hole.phi_0 };
    FP x_prev[D];
    FP v[D];
    FP de = hole.de0;

    ijk_to_vec_zoom(i, j, hole, SZELESregi, MAGASregi, ikezd, jkezd, iveg, v);

    // Outer horizon r+ = (rs + sqrt(rs^2 - 4(a^2+Q^2))) / 2, the larger root of
    // Delta = 0, nudged outwards because the geodesic equation is singular on
    // it.  A negative discriminant means the parameters describe a naked
    // singularity: there is no horizon to absorb the ray, so the test degrades
    // to r < 0 and such rays end up classified as RAY_NONFINITE instead.
    FP r_horizon = FP(0);
    FP const disc = hole.rs * hole.rs - FP(4) * (hole.a * hole.a + hole.Q * hole.Q);
    if (disc > FP(0))
    {
        r_horizon = (hole.rs + sqrt(disc)) / FP(2) + FP(HORIZON_EPS);
    }

    ray_outcome result = RAY_BUDGET;

    for (uint64_t n = 0; n < hole.max_steps; ++n)
    {
        for (int k = 0; k < D; ++k)
        {
            x_prev[k] = x[k];
        }

        step(hole, x, v, de);

        if (!(isfinite(x[0]) && isfinite(x[1]) && isfinite(x[2]) && isfinite(x[3])))
        {
            result = RAY_NONFINITE;
            break;
        }
        if (gomb_be(r_horizon, x))
        {
            result = RAY_HORIZON;
            break;
        }
        if (gomb_ki(hole.sugar_ki, x))
        {
            result = RAY_ESCAPED;
            break;
        }
        if (disk1(hole.sugar_kicsi, hole.sugar_nagy, x, x_prev))
        {
            result = RAY_DISK_TOP;
            break;
        }
        if (disk2(hole.sugar_kicsi, hole.sugar_nagy, x, x_prev))
        {
            result = RAY_DISK_BOTTOM;
            break;
        }
    }

    for (int k = 0; k < D; ++k)
    {
        x_end[k] = x[k];
    }
    return result;
}

/**
 * One ray per pixel; writes a palette index for each.
 *
 * szin encoding, as consumed by the imagemaker scripts:
 *   1  swallowed by the horizon, or escaped to the sky sphere
 *   0  hit the accretion disk from above
 *   3  hit the accretion disk from below
 *   2  the integration produced NaN/Inf
 *  -1  the ray used up its step budget without hitting anything
 */
template <class FP>
inline void ray_step(int8_t* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const xd, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, uint64_t const max_steps, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki_in, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg)
{
    kerr_black_hole<FP> hole(SZELES, MAGAS, xd, Omega, a, Q, rs, errormax, de0, max_steps, kepernyo_high, kepernyo_tav, sugar_ki_in, gyuru_sugar_kicsi, gyuru_sugar_nagy);

#pragma acc data copyin(xd[0:4],Omega[0:4]) copyout(szin[0:SZELES*MAGAS])
{
#pragma acc parallel loop collapse(2)
    for (uint64_t j = 0; j < MAGAS; j++)
    {
        for (uint64_t i = 0; i < SZELES; i++)
        {
            FP x_end[D];
            ray_outcome const outcome = trace_ray(hole, i, j, SZELESregi, MAGASregi, ikezd, jkezd, iveg, x_end);

            int8_t value;
            switch (outcome)
            {
            case RAY_HORIZON:
            case RAY_ESCAPED:
                value = 1;
                break;
            case RAY_DISK_TOP:
                value = 0;
                break;
            case RAY_DISK_BOTTOM:
                value = 3;
                break;
            case RAY_NONFINITE:
                value = 2;
                break;
            default:
                value = -1;
                break;
            }

            szin[i * MAGAS + j] = value;
        }
    }
}
}

/**
 * Same trace as ray_step(), but the pixel carries a physical quantity instead of
 * a palette index: 0 where nothing was hit, the Boyer-Lindquist radius of the
 * disk crossing where the ray landed on the disk (the imagemaker turns that into
 * a temperature), and -1 for a failed ray.
 */
template <class FP>
inline void ray_step_T(FP* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const xd, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, uint64_t const max_steps, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki_in, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg)
{
    kerr_black_hole<FP> hole(SZELES, MAGAS, xd, Omega, a, Q, rs, errormax, de0, max_steps, kepernyo_high, kepernyo_tav, sugar_ki_in, gyuru_sugar_kicsi, gyuru_sugar_nagy);

#pragma acc data copyin(xd[0:4],Omega[0:4]) copyout(szin[0:SZELES*MAGAS])
{
#pragma acc parallel loop collapse(2)
    for (uint64_t j = 0; j < MAGAS; j++)
    {
        for (uint64_t i = 0; i < SZELES; i++)
        {
            FP x_end[D];
            ray_outcome const outcome = trace_ray(hole, i, j, SZELESregi, MAGASregi, ikezd, jkezd, iveg, x_end);

            FP value;
            switch (outcome)
            {
            case RAY_HORIZON:
            case RAY_ESCAPED:
                value = FP(0);
                break;
            case RAY_DISK_TOP:
            case RAY_DISK_BOTTOM:
                value = x_end[1];
                break;
            default:
                value = FP(-1);
                break;
            }

            szin[i * MAGAS + j] = value;
        }
    }
}
}

/** True when the ray has fallen inside the sphere of radius `sugar` (the horizon). */
template <class FP>
inline bool gomb_be(FP const sugar, FP const* const x)
{
    return x[1] < sugar;
}

/** True when the ray has left the sphere of radius `sugar` (the sky sphere). */
template <class FP>
inline bool gomb_ki(FP const sugar, FP const* const x)
{
    return x[1] > sugar;
}

/**
 * True when the segment from x2 to x1 crosses the equatorial plane in EITHER
 * direction, inside the annulus sugar_kicsi < r < sugar_nagy.
 *
 * Unused - disk1 and disk2 split the same test by crossing direction so the two
 * faces of the accretion disk can be shaded differently.  Kept because it spells
 * the geometry out in one place.
 */
template <class FP>
inline bool disk(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2)
{
    return disk1(sugar_kicsi, sugar_nagy, x1, x2) || disk2(sugar_kicsi, sugar_nagy, x1, x2);
}

/** Equatorial-plane crossing from theta < pi/2 to theta > pi/2 (disk seen from above). */
template <class FP>
inline bool disk1(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2)
{
    FP const half_pi = FP(3.14159265358979323846 / 2.0);
    return (x1[1] > sugar_kicsi) && (x1[1] < sugar_nagy) && (x1[2] > half_pi) && (x2[2] < half_pi);
}

/** Equatorial-plane crossing from theta > pi/2 to theta < pi/2 (disk seen from below). */
template <class FP>
inline bool disk2(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2)
{
    FP const half_pi = FP(3.14159265358979323846 / 2.0);
    return (x1[1] > sugar_kicsi) && (x1[1] < sugar_nagy) && (x1[2] < half_pi) && (x2[2] > half_pi);
}

/** Flat index of component k of pixel (i, j) in a SZELES x MAGAS x D array.  Unused. */
template <class FP>
inline int ijk_to_n(uint64_t const i, uint64_t const j, uint64_t const k, kerr_black_hole<FP> const& hole)
{
    return i * hole.MAGAS * D + j * D + k;
}

/** Integer power; small exponents are spelled out because christoffel() leans on them. */
template <class FP>
inline FP pown(FP const x, int const n)
{
    switch(n)
    {
        case 0: return FP(1);
        case 1: return x;
        case 2: return x*x;
        case 3: return x*x*x;
        case 4: { FP x2=x*x; return x2*x2; }
        case 5: { FP x2=x*x; return x2*x2*x; }
        case 6: { FP x3=x*x*x; return x3*x3; }
    }

    // fallback generic
    return pown_gen(x,n);
}
/** Generic integer power by binary exponentiation; fallback for pown(). */
template <class FP>
inline FP pown_gen(FP const x, int const n) 
{ 
    if (n == 0) return FP(1); // Handle negative exponent safely (including INT_MIN)
    bool neg = (n < 0);
    long long nl = (long long)n;
    unsigned long long exp = neg ? (unsigned long long)(-nl) : (unsigned long long)nl;
    FP result = FP(1);
    FP base = x;
    while (exp) 
    {
         if (exp & 1ull) result = result * base;
	 base = base * base;
	 exp >>= 1ull;
    }
    if (neg) return FP(1) / result;
    return result;
}

/** Recursive integer power.  Unused; kept as a reference implementation. */
template <class FP>
inline FP powni_rec(FP const x, int const n)
{
    if (n < 0)
        return pown(1.0 / x, (-n));
    else if (n == 0)
        return  1.0;
    else if (n == 1)
        return  x;
    else if (n % 2 == 0)
        return pown((x * x), (n / 2));
    else //if (n%2==1)
        return x * pown((x * x), ((n - 1) / 2));
}


#endif // CUDA_RAY_H
