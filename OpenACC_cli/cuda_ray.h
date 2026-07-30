#ifndef CUDA_RAY_H
#define CUDA_RAY_H

#include <iostream>
#include <math.h>

#include "black_hole.h"

//#include "debugmalloc.h"


template <class FP>
inline void step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de);
template <class FP>
inline void dopri54_step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de);

template <class FP>
inline FP ijk_to_vec_mink_zoom(uint64_t const i, uint64_t const j, uint64_t const k, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg, kerr_black_hole<FP> const& hole);
template <class FP>
inline FP ijk_to_vec_zoom(uint64_t const i, uint64_t const j, uint64_t const k, kerr_black_hole<FP> const& hole, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg);



template <class FP>
inline void RK38(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de);//4-ed foku legpontosabb
template <class FP>
inline void RK6(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de);

template <class FP>
inline void christoffel(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ ch);

template <class FP>
inline void ray_step(int8_t* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const x, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg);
template <class FP>
inline void ray_step_T(FP* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const x, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg);

// Ratio of photon frequency measured by the camera to that measured by a
// prograde circular emitter in the equatorial Kerr-Newman disk.  The ray's
// covariant energy and angular momentum are conserved, so this includes both
// gravitational and special-relativistic Doppler shifts.
template <class FP>
inline FP disk_redshift(kerr_black_hole<FP> const& hole, FP const* const x, FP const* const v);




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
inline FP pown_rec(FP const x, int const n);


template <class FP>
inline void step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)//adaptiv step size
{
    dopri54_step(hole, x, v, de);
}

// Embedded Dormand-Prince 5(4) controller for the first-order form of the
// geodesic equation: x' = v, v' = Gamma(x, v).  Unlike the earlier
// curvature-magnitude heuristic, the 5th-vs-4th-order difference is a local
// truncation-error estimate.  It also removes the old post-step Christoffel
// evaluation (7 evaluations per accepted attempt instead of 8).
template <class FP>
inline void dopri54_step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)
{
    FP const max_step = hole.de0;
    FP const min_step = fmax(FP(1e-8), max_step * FP(1e-5));
    FP const tolerance = fmax(hole.errormax, FP(1e-8));
    FP h = fmin(max_step, fmax(min_step, de));
    FP kx[7][D], kv[7][D], trial_x[D], trial_v[D], fifth_x[D], fifth_v[D], fourth_x[D], fourth_v[D], acceleration[D];

    for (int attempt = 0; attempt < 10; ++attempt)
    {
        christoffel(hole, x, v, acceleration);
        for (int n = 0; n < D; ++n) { kx[0][n] = v[n]; kv[0][n] = acceleration[n]; }

        for (int n = 0; n < D; ++n) { trial_x[n] = x[n] + h * FP(1.0/5.0) * kx[0][n]; trial_v[n] = v[n] + h * FP(1.0/5.0) * kv[0][n]; }
        christoffel(hole, trial_x, trial_v, acceleration);
        for (int n = 0; n < D; ++n) { kx[1][n] = trial_v[n]; kv[1][n] = acceleration[n]; }

        for (int n = 0; n < D; ++n) { trial_x[n] = x[n] + h * (FP(3.0/40.0) * kx[0][n] + FP(9.0/40.0) * kx[1][n]); trial_v[n] = v[n] + h * (FP(3.0/40.0) * kv[0][n] + FP(9.0/40.0) * kv[1][n]); }
        christoffel(hole, trial_x, trial_v, acceleration);
        for (int n = 0; n < D; ++n) { kx[2][n] = trial_v[n]; kv[2][n] = acceleration[n]; }

        for (int n = 0; n < D; ++n) { trial_x[n] = x[n] + h * (FP(44.0/45.0) * kx[0][n] - FP(56.0/15.0) * kx[1][n] + FP(32.0/9.0) * kx[2][n]); trial_v[n] = v[n] + h * (FP(44.0/45.0) * kv[0][n] - FP(56.0/15.0) * kv[1][n] + FP(32.0/9.0) * kv[2][n]); }
        christoffel(hole, trial_x, trial_v, acceleration);
        for (int n = 0; n < D; ++n) { kx[3][n] = trial_v[n]; kv[3][n] = acceleration[n]; }

        for (int n = 0; n < D; ++n) { trial_x[n] = x[n] + h * (FP(19372.0/6561.0) * kx[0][n] - FP(25360.0/2187.0) * kx[1][n] + FP(64448.0/6561.0) * kx[2][n] - FP(212.0/729.0) * kx[3][n]); trial_v[n] = v[n] + h * (FP(19372.0/6561.0) * kv[0][n] - FP(25360.0/2187.0) * kv[1][n] + FP(64448.0/6561.0) * kv[2][n] - FP(212.0/729.0) * kv[3][n]); }
        christoffel(hole, trial_x, trial_v, acceleration);
        for (int n = 0; n < D; ++n) { kx[4][n] = trial_v[n]; kv[4][n] = acceleration[n]; }

        for (int n = 0; n < D; ++n) { trial_x[n] = x[n] + h * (FP(9017.0/3168.0) * kx[0][n] - FP(355.0/33.0) * kx[1][n] + FP(46732.0/5247.0) * kx[2][n] + FP(49.0/176.0) * kx[3][n] - FP(5103.0/18656.0) * kx[4][n]); trial_v[n] = v[n] + h * (FP(9017.0/3168.0) * kv[0][n] - FP(355.0/33.0) * kv[1][n] + FP(46732.0/5247.0) * kv[2][n] + FP(49.0/176.0) * kv[3][n] - FP(5103.0/18656.0) * kv[4][n]); }
        christoffel(hole, trial_x, trial_v, acceleration);
        for (int n = 0; n < D; ++n) { kx[5][n] = trial_v[n]; kv[5][n] = acceleration[n]; }

        for (int n = 0; n < D; ++n)
        {
            fifth_x[n] = x[n] + h * (FP(35.0/384.0) * kx[0][n] + FP(500.0/1113.0) * kx[2][n] + FP(125.0/192.0) * kx[3][n] - FP(2187.0/6784.0) * kx[4][n] + FP(11.0/84.0) * kx[5][n]);
            fifth_v[n] = v[n] + h * (FP(35.0/384.0) * kv[0][n] + FP(500.0/1113.0) * kv[2][n] + FP(125.0/192.0) * kv[3][n] - FP(2187.0/6784.0) * kv[4][n] + FP(11.0/84.0) * kv[5][n]);
        }
        christoffel(hole, fifth_x, fifth_v, acceleration);
        for (int n = 0; n < D; ++n) { kx[6][n] = fifth_v[n]; kv[6][n] = acceleration[n]; }

        FP error_norm = FP(0);
        for (int n = 0; n < D; ++n)
        {
            fourth_x[n] = x[n] + h * (FP(5179.0/57600.0) * kx[0][n] + FP(7571.0/16695.0) * kx[2][n] + FP(393.0/640.0) * kx[3][n] - FP(92097.0/339200.0) * kx[4][n] + FP(187.0/2100.0) * kx[5][n] + FP(1.0/40.0) * kx[6][n]);
            fourth_v[n] = v[n] + h * (FP(5179.0/57600.0) * kv[0][n] + FP(7571.0/16695.0) * kv[2][n] + FP(393.0/640.0) * kv[3][n] - FP(92097.0/339200.0) * kv[4][n] + FP(187.0/2100.0) * kv[5][n] + FP(1.0/40.0) * kv[6][n]);
            FP const x_scale = tolerance * (FP(1) + fmax(fabs(x[n]), fabs(fifth_x[n])));
            FP const v_scale = tolerance * (FP(1) + fmax(fabs(v[n]), fabs(fifth_v[n])));
            error_norm = fmax(error_norm, fabs(fifth_x[n] - fourth_x[n]) / x_scale);
            error_norm = fmax(error_norm, fabs(fifth_v[n] - fourth_v[n]) / v_scale);
        }

        FP factor = error_norm > FP(0) ? FP(0.9) * pow(error_norm, FP(-0.2)) : FP(5);
        factor = fmin(FP(5), fmax(FP(0.2), factor));
        if (error_norm <= FP(1) || h <= min_step || attempt == 9)
        {
            for (int n = 0; n < D; ++n) { x[n] = fifth_x[n]; v[n] = fifth_v[n]; }
            de = fmin(max_step, fmax(min_step, h * factor));
            return;
        }
        h = fmax(min_step, h * factor);
    }
}



template <class FP>
inline void RK38(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)
{
    FP ch[D];

    christoffel(hole, x, v, ch);

    FP kx1[D];
    FP kv1[D];

    for (int i = 0; i < D; ++i)
    {
        kx1[i] = v[i];
    }
    for (int i = 0; i < D; ++i)
    {
        kv1[i] = ch[i];
    }

    FP x1[D];
    FP v1[D];

    for (int i = 0; i < D; ++i)
    {
        x1[i] = x[i] + kx1[i] * (de / 3.0);
    }
    for (int i = 0; i < D; ++i)
    {
        v1[i] = v[i] + kv1[i] * (de / 3.0);
    }

    christoffel(hole, x1, v1, ch);

    FP kx2[D];
    FP kv2[D];

    for (int i = 0; i < D; ++i)
    {
        kx2[i] = v1[i];
    }
    for (int i = 0; i < D; ++i)
    {
        kv2[i] = ch[i];
    }

    for (int i = 0; i < D; ++i)
    {
        x1[i] = x[i] - kx1[i] * (de / 3.0) + kx2[i] * de;
    }
    for (int i = 0; i < D; ++i)
    {
        v1[i] = v[i] - kv1[i] * (de / 3.0) + kv2[i] * de;
    }

    christoffel(hole, x1, v1, ch);

    FP kx3[D];
    FP kv3[D];

    for (int i = 0; i < D; ++i)
    {
        kx3[i] = v1[i];
    }
    for (int i = 0; i < D; ++i)
    {
        kv3[i] = ch[i];
    }


    for (int i = 0; i < D; ++i)
    {
        x1[i] = x[i] + (kx1[i] - kx2[i] + kx3[i]) * de;
    }
    for (int i = 0; i < D; ++i)
    {
        v1[i] = v[i] + (kv1[i] - kv2[i] + kv3[i]) * de;
    }

    christoffel(hole, x1, v1, ch);

    FP kx4[D];
    FP kv4[D];

    for (int i = 0; i < D; ++i)
    {
        kx4[i] = v1[i];
    }
    for (int i = 0; i < D; ++i)
    {
        kv4[i] = ch[i];
    }

    for (int i = 0; i < D; ++i)
    {
        x[i] = x[i] + (kx1[i] / 8.0 + kx2[i] * (3.0 / 8.0) + kx3[i] * (3.0 / 8.0) + kx4[i] / 8.0) * de;
        v[i] = v[i] + (kv1[i] / 8.0 + kv2[i] * (3.0 / 8.0) + kv3[i] * (3.0 / 8.0) + kv4[i] / 8.0) * de;
    }
}

template <class FP>
inline void RK6(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)
{
    FP ch[D];

    christoffel(hole, x, v, ch);

    FP kx1[D];
    FP kv1[D];

    for (int i = 0; i < D; ++i)
    {
        kx1[i] = v[i];
    }
    for (int i = 0; i < D; ++i)
    {
        kv1[i] = ch[i];
    }

    FP x1[D];
    FP v1[D];

    for (int i = 0; i < D; ++i)
    {
        x1[i] = x[i] + kx1[i] * (de / 5.0);
    }
    for (int i = 0; i < D; ++i)
    {
        v1[i] = v[i] + kv1[i] * (de / 5.0);
    }

    christoffel(hole, x1, v1, ch);

    FP kx2[D];
    FP kv2[D];

    for (int i = 0; i < D; ++i)
    {
        kx2[i] = v1[i];
    }
    for (int i = 0; i < D; ++i)
    {
        kv2[i] = ch[i];
    }

    for (int i = 0; i < D; ++i)
    {
        x1[i] = x[i] + kx1[i] * (de * (3.0 / 40.0)) + kx2[i] * (de * (9.0 / 40.0));
    }
    for (int i = 0; i < D; ++i)
    {
        v1[i] = v[i] + kv1[i] * (de * (3.0 / 40.0)) + kv2[i] * (de * (9.0 / 40.0));
    }

    christoffel(hole, x1, v1, ch);

    FP kx3[D];
    FP kv3[D];

    for (int i = 0; i < D; ++i)
    {
        kx3[i] = v1[i];
    }
    for (int i = 0; i < D; ++i)
    {
        kv3[i] = ch[i];
    }


    for (int i = 0; i < D; ++i)
    {
        x1[i] = x[i] + (kx1[i] * (44.0 / 45.0) - kx2[i] * (56.0 / 15.0) + kx3[i] * (32.0 / 9.0)) * de;
    }
    for (int i = 0; i < D; ++i)
    {
        v1[i] = v[i] + (kv1[i] * (44.0 / 45.0) - kv2[i] * (56.0 / 15.0) + kv3[i] * (32.0 / 9.0)) * de;
    }

    christoffel(hole, x1, v1, ch);

    FP kx4[D];
    FP kv4[D];

    for (int i = 0; i < D; ++i)
    {
        kx4[i] = v1[i];
    }
    for (int i = 0; i < D; ++i)
    {
        kv4[i] = ch[i];
    }


    for (int i = 0; i < D; ++i)
    {
        x1[i] = x[i] + (kx1[i] * (19372.0 / 6561.0) - kx2[i] * (25360.0 / 2187.0) + kx3[i] * (64448.0 / 6561.0) - kx4[i] * (212.0 / 729.0)) * de;
    }
    for (int i = 0; i < D; ++i)
    {
        v1[i] = v[i] + (kv1[i] * (19372.0 / 6561.0) - kv2[i] * (25360.0 / 2187.0) + kv3[i] * (64448.0 / 6561.0) - kv4[i] * (212.0 / 729.0)) * de;
    }

    christoffel(hole, x1, v1, ch);

    FP kx5[D];
    FP kv5[D];

    for (int i = 0; i < D; ++i)
    {
        kx5[i] = v1[i];
    }
    for (int i = 0; i < D; ++i)
    {
        kv5[i] = ch[i];
    }



    for (int i = 0; i < D; ++i)
    {
        x1[i] = x[i] + (kx1[i] * (9017.0 / 3168.0) - kx2[i] * (355.0 / 33.0) + kx3[i] * (46732.0 / 5147.0) + kx4[i] * (49.0 / 176.0) - kx5[i] * (5103.0 / 18656.0)) * de;
    }
    for (int i = 0; i < D; ++i)
    {
        v1[i] = v[i] + (kv1[i] * (9017.0 / 3168.0) - kv2[i] * (355.0 / 33.0) + kv3[i] * (46732.0 / 5147.0) + kv4[i] * (49.0 / 176.0) - kv5[i] * (5103.0 / 18656.0)) * de;
    }

    christoffel(hole, x1, v1, ch);

    FP kx6[D];
    FP kv6[D];

    for (int i = 0; i < D; ++i)
    {
        kx6[i] = v1[i];
    }
    for (int i = 0; i < D; ++i)
    {
        kv6[i] = ch[i];
    }


    for (int i = 0; i < D; ++i)
    {
        x1[i] = x[i] + (kx1[i] * (35.0 / 384.0) + kx3[i] * (500.0 / 1113.0) + kx4[i] * (125.0 / 192.0) - kx5[i] * (2187.0 / 6784.0) + kx6[i] * (11.0 / 84.0)) * de;
    }
    for (int i = 0; i < D; ++i)
    {
        v1[i] = v[i] + (kv1[i] * (35.0 / 384.0) + kv3[i] * (500.0 / 1113.0) + kv4[i] * (125.0 / 192.0) - kv5[i] * (2187.0 / 6784.0) + kv6[i] * (11.0 / 84.0)) * de;
    }

    christoffel(hole, x1, v1, ch);

    FP kx7[D];
    FP kv7[D];

    for (int i = 0; i < D; ++i)
    {
        kx7[i] = v1[i];
    }
    for (int i = 0; i < D; ++i)
    {
        kv7[i] = ch[i];
    }



    for (int i = 0; i < D; ++i)
    {
        x[i] = x[i] + (kx1[i] * (5179.0 / 57600.0) + kx3[i] * (7571.0 / 16695.0) + kx4[i] * (393.0 / 640.0) - kx5[i] * (92097.0 / 339200.0) + kx6[i] * (187.0 / 2100.0) + kx7[i] * (1.0 / 40.0)) * de;
        v[i] = v[i] + (kv1[i] * (5179.0 / 57600.0) + kv3[i] * (7571.0 / 16695.0) + kv4[i] * (393.0 / 640.0) - kv5[i] * (92097.0 / 339200.0) + kv6[i] * (187.0 / 2100.0) + kv7[i] * (1.0 / 40.0)) * de;
    }
}

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

template <class FP>
inline FP ijk_to_vec_mink_zoom(uint64_t i, uint64_t j, uint64_t k, uint64_t SZELESregi, uint64_t MAGASregi, uint64_t ikezd, uint64_t jkezd, uint64_t iveg, kerr_black_hole<FP>& hole)//futtassuk le zoom nélkül és válasszunk ki egy pontot ez lesz az új kép bal felsõ sarka ezek az ikezd,jkezd// a jobb alsó sarok pedig a iveg,jveg de jveg a SZELES MAGAS arányából következik
{
    if (k == 0)
    {
        return 1.0;//idõszerû komponens
    }
    else//térszerû komponensek
    {
        double x = hole.kepernyo_tav;
        double kemernyo_high = hole.kepernyo_high;
        int SZELES = hole.SZELES;
        //int MAGAS = hole.MAGAS;

        FP ir = FP(ikezd) + (FP(i) / FP(SZELES)) * (FP(iveg) - FP(ikezd));
        FP jr = FP(jkezd) + (FP(j) / FP(SZELES)) * (FP(iveg) - FP(ikezd));//igen igy jo csak bele kell gondolni SZELES/MAGAS=(iveg-ikezd)/(jveg-jkezd)

        FP y = (kemernyo_high / MAGASregi) * (FP(MAGASregi) / 2 - FP(jr));
        FP z = (kemernyo_high / MAGASregi) * (FP(ir) - FP(SZELESregi) / 2);

        FP norm = sqrt(x * x + y * y + z * z);


        if (k == 1)
        {
            return x / norm;
        }
        else if (k == 2)
        {
            return y / norm;
        }
        else if (k == 3)
        {
            return z / norm;
        }
        else
        {
            //std::cout << "ijk_to_vec_mink fv.-ben k tul lett indexelve\ni=" << i << "\nj=" << j << "\nk=" << k << "\n";
            return 0;
        }
    }

}

template <class FP>
inline FP ijk_to_vec_zoom(uint64_t i, uint64_t j, uint64_t k, kerr_black_hole<FP>& hole, uint64_t SZELESregi, uint64_t MAGASregi, uint64_t ikezd, uint64_t jkezd, uint64_t iveg)
{
    //elforgatás és át transzformálás a görbült téridõ adott pontjában lakó vektorokra
    //....
    FP x[4] = { ijk_to_vec_mink_zoom(i, j, 0,SZELESregi,MAGASregi,ikezd,jkezd,iveg,hole), ijk_to_vec_mink_zoom(i, j, 1,SZELESregi,MAGASregi,ikezd,jkezd,iveg,hole),ijk_to_vec_mink_zoom(i, j, 2,SZELESregi,MAGASregi,ikezd,jkezd,iveg,hole),ijk_to_vec_mink_zoom(i, j, 3,SZELESregi,MAGASregi,ikezd,jkezd,iveg,hole) };

    FP phi = sqrt(hole.Omega_1 * hole.Omega_1 + hole.Omega_2 * hole.Omega_2 + hole.Omega_3 * hole.Omega_3);//180.0 / 180.0 * 3.14159265;//172.6 / 180.0 * 3.14159265;//
    FP u[D - 1] = { hole.Omega_1 / phi,hole.Omega_2 / phi,hole.Omega_3 / phi };

    if (phi == 0.0)
    {
        u[0] = 1.0;
        u[1] = 0.0;
        u[2] = 0.0;
    }



    //forgatás



    FP x1_tmp = (cos(phi) + u[0] * u[0] * (1 - cos(phi))) * x[1] + (u[0] * u[1] * (1 - cos(phi)) - u[2] * sin(phi)) * x[2] + (u[0] * u[2] * (1 - cos(phi)) + u[1] * sin(phi)) * x[3];
    FP x2_tmp = (u[0] * u[1] * (1 - cos(phi) + u[2] * sin(phi))) * x[1] + (cos(phi) + u[1] * u[1] * (1 - cos(phi))) * x[2] + (u[1] * u[2] * (1 - cos(phi) + u[0] * sin(phi))) * x[3];
    FP x3_tmp = (u[0] * u[2] * (1 - cos(phi)) - u[1] * sin(phi)) * x[1] + (u[1] * u[2] * (1 - cos(phi)) + u[0] * sin(phi)) * x[2] + (cos(phi) + u[2] * u[2] * (1 - cos(phi))) * x[3];

    x[1]=x1_tmp;
    x[2]=x2_tmp;
    x[3]=x3_tmp;

    FP r_0 = hole.r_0;
    FP theta_0 = hole.theta_0;
    FP rs = hole.rs;
    FP a = hole.a;
    FP Q = hole.Q;

    FP delta = r_0 * r_0 - 4 * rs * r_0 + a * a + Q * Q;
    FP rho = sqrt(r_0 * r_0 + a * a * cos(theta_0) * cos(theta_0));

    FP x0_tmp = x[0] * (a * a + r_0 * r_0) * rho / ((a * a * cos(theta_0) * cos(theta_0) + r_0 * r_0) * sqrt(delta)) + x[3] * a * rho / (sqrt(delta) * (a * a * cos(theta_0) * cos(theta_0) + r_0 * r_0));
    x1_tmp = sqrt(delta) / rho * x[1];
    x2_tmp = x[2] / rho;
    x3_tmp = x[0] * a * rho * sin(theta_0) / (a * a * cos(theta_0) * cos(theta_0) + r_0 * r_0) + x[3] * rho / (sin(theta_0) * (r_0 * r_0 + a * a * cos(theta_0) * cos(theta_0)));

    x[1]=x1_tmp;
    x[2]=x2_tmp;
    x[3]=x3_tmp;
    x[0]=x0_tmp;


    return x[k];
}

template <class FP>
inline void ray_step(int8_t* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const xd, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki_in, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg)//kernel
{
    kerr_black_hole<FP> hole(SZELES, MAGAS, xd, Omega, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki_in, gyuru_sugar_kicsi, gyuru_sugar_nagy);

#pragma acc data copyin(xd[0:4],Omega[0:4]) copyout(szin[0:SZELES*MAGAS])
{
#pragma acc parallel loop collapse(2)
    for(uint64_t j=0; j<MAGAS; j++)
    {
//#pragma acc loop
        for(uint64_t i=0; i<SZELES; i++)
        {
            FP x[D] = { hole.t_0,hole.r_0,hole.theta_0,hole.phi_0 };;
            FP v[D];
            FP x_le[D]; //lemaradó hely koordináták
            FP de = de0;

            for (int k = 0; k < D; ++k)
            {
                v[k] = ijk_to_vec_zoom(i, j, k, hole, SZELESregi, MAGASregi, ikezd, jkezd, iveg);
            }




            for (int k = 0; k < D; ++k)
            {
                x_le[k] = x[k];
            }
            step(hole, x, v, de);

            bool fut = true;
            FP sugar_be = 0;

            if ((rs * rs - 4 * (a * a + Q * Q)) > 0.0)
            {
                sugar_be = (rs + sqrt(rs * rs - 4 * (a * a + Q * Q))) / 2 + 0.001;
            }
            else
            {

            }

            //(rs - sqrt(rs * rs - 4 * a * a * cos(x[2]) * cos(x[2]))) / 2 + 0.0001;






            FP sugar_ki = hole.sugar_ki;

            FP sugar_kicsi = hole.sugar_kicsi;
            FP sugar_nagy = hole.sugar_nagy;

            //itt következik az ütközés detektálás

            int idokorlat = 0;

            while (fut)
            {
                if (gomb_be(sugar_be, x))
                {
                    szin[i * MAGAS + j] = 1;//-1;
                    fut = false;
                }
                else if (gomb_ki(sugar_ki, x))
                {
                    szin[i * MAGAS + j] = 1;
                    fut = false;
                }
                else if (disk1(sugar_kicsi, sugar_nagy, x, x_le))
                {
                    szin[i * MAGAS + j] = 0;
                    fut = false;
                }
                else if (disk2(sugar_kicsi, sugar_nagy, x, x_le))
                {
                    szin[i * MAGAS + j] = 3;
                    fut = false;
                }
                else if (isnan(x[0]) || isnan(x[1]) || isnan(x[2]) || isnan(x[3]))
                {
                    //printf("%d\t%d\t%f\tnan\n", i, j, de);
                    szin[i * MAGAS + j] = 2;
                    fut = false;
                }
                else if (isinf(x[0]) || isinf(x[1]) || isinf(x[2]) || isinf(x[3]))
                {
                    //printf("%d\t%d\t%f\tinf\n", i, j, de);
                    szin[i * MAGAS + j] = 2;
                    fut = false;
                }
                else
                {

                }

                ++idokorlat;
                if (idokorlat >= int(1.0 / errormax))//if (idokorlat >= int(1.0 / errormax))
                {
                    //printf("%d\t%d\t%f\tmegunta\n", i, j, de);
                    szin[i * MAGAS + j] = -1;
                    fut = false;
                }


                for (int k = 0; k < D; ++k)
                {
                    x_le[k] = x[k];
                }
                step(hole, x, v, de);




            }



            //printf("%d\t%d\t%f\n", i, j, de);
        }
    }
}
}

template <class FP>
inline void ray_step_T(FP* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const xd, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki_in, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg)//kernel
{
    kerr_black_hole<FP> hole(SZELES, MAGAS, xd, Omega, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki_in, gyuru_sugar_kicsi, gyuru_sugar_nagy);

#pragma acc data copyin(xd[0:4],Omega[0:4]) copyout(szin[0:3*SZELES*MAGAS])
{
#pragma acc parallel loop collapse(2)
    for(uint64_t j=0; j<MAGAS; j++)
    {
//#pragma acc loop
        for(uint64_t i=0; i<SZELES; i++)
        {

            FP x[D] = { hole.t_0,hole.r_0,hole.theta_0,hole.phi_0 };;
            FP v[D];
            FP x_le[D]; //lemaradó hely koordináták
            FP de = de0;

            for (int k = 0; k < D; ++k)
            {
                v[k] = ijk_to_vec_zoom(i, j, k, hole, SZELESregi, MAGASregi, ikezd, jkezd, iveg);
            }




            for (int k = 0; k < D; ++k)
            {
                x_le[k] = x[k];
            }
            step(hole, x, v, de);

            bool fut = true;
            FP sugar_be = 0;

            if ((rs * rs - 4 * (a * a + Q * Q)) > 0.0)
            {
                sugar_be = (rs + sqrt(rs * rs - 4 * (a * a + Q * Q))) / 2 + 0.001;
            }
            else
            {

            }

            //(rs - sqrt(rs * rs - 4 * a * a * cos(x[2]) * cos(x[2]))) / 2 + 0.0001;






            FP sugar_ki = hole.sugar_ki;

            FP sugar_kicsi = hole.sugar_kicsi;
            FP sugar_nagy = hole.sugar_nagy;

            //itt következik az ütközés detektálás

            int idokorlat = 0;

            while (fut)
            {
                //0 fekete, -1 hiba kezeléses piros, egyébként meg egy FP ami reprezental egy szint
                if (gomb_be(sugar_be, x))
                {
                    szin[3 * (i * MAGAS + j)] = 0;//;
                    szin[3 * (i * MAGAS + j) + 1] = 0;
                    szin[3 * (i * MAGAS + j) + 2] = 0;
                    fut = false;
                }
                else if (gomb_ki(sugar_ki, x))
                {
                    // The ray reached the outer integration sphere: it
                    // escaped to the sky, rather than being captured.
                    szin[3 * (i * MAGAS + j)] = -1;
                    szin[3 * (i * MAGAS + j) + 1] = 0;
                    szin[3 * (i * MAGAS + j) + 2] = 0;
                    fut = false;
                }
                else if (disk1(sugar_kicsi, sugar_nagy, x, x_le))
                {
                    szin[3 * (i * MAGAS + j)] = x[1];
                    szin[3 * (i * MAGAS + j) + 1] = disk_redshift(hole, x, v);
                    szin[3 * (i * MAGAS + j) + 2] = x[3];
                    fut = false;
                }
                else if (disk2(sugar_kicsi, sugar_nagy, x, x_le))
                {
                    szin[3 * (i * MAGAS + j)] = x[1];
                    szin[3 * (i * MAGAS + j) + 1] = disk_redshift(hole, x, v);
                    szin[3 * (i * MAGAS + j) + 2] = x[3];
                    fut = false;
                }
                else if (isnan(x[0]) || isnan(x[1]) || isnan(x[2]) || isnan(x[3]))
                {
                    //printf("%d\t%d\t%f\tnan\n", i, j, de);
                    szin[3 * (i * MAGAS + j)] = -1;
                    szin[3 * (i * MAGAS + j) + 1] = 0;
                    szin[3 * (i * MAGAS + j) + 2] = 0;
                    fut = false;
                }
                else if (isinf(x[0]) || isinf(x[1]) || isinf(x[2]) || isinf(x[3]))
                {
                    //printf("%d\t%d\t%f\tinf\n", i, j, de);
                    szin[3 * (i * MAGAS + j)] = -1;
                    szin[3 * (i * MAGAS + j) + 1] = 0;
                    szin[3 * (i * MAGAS + j) + 2] = 0;
                    fut = false;
                }
                else
                {

                }

                ++idokorlat;
                if (idokorlat >= int(1.0 / errormax))//if (idokorlat >= int(1.0 / errormax))
                {
                    //printf("%d\t%d\t%f\tmegunta\n", i, j, de);
                    szin[3 * (i * MAGAS + j)] = -1;
                    szin[3 * (i * MAGAS + j) + 1] = 0;
                    szin[3 * (i * MAGAS + j) + 2] = 0;
                    fut = false;
                }


                for (int k = 0; k < D; ++k)
                {
                    x_le[k] = x[k];
                }
                step(hole, x, v, de);




            }



            //printf("%d\t%d\t%f\n", i, j, de);

        }
    }
}    
}

template <class FP>
inline FP disk_redshift(kerr_black_hole<FP> const& hole, FP const* const x, FP const* const v)
{
    // Boyer-Lindquist Kerr-Newman metric coefficients at the disk hit.
    FP const r = x[1];
    FP const theta = x[2];
    FP const sin_theta = sin(theta);
    FP const sin2 = sin_theta * sin_theta;
    FP const sigma = r * r + hole.a * hole.a * cos(theta) * cos(theta);
    FP const delta = r * r - hole.rs * r + hole.a * hole.a + hole.Q * hole.Q;
    FP const gtt = -(FP(1) - (hole.rs * r - hole.Q * hole.Q) / sigma);
    FP const gtphi = -hole.a * (hole.rs * r - hole.Q * hole.Q) * sin2 / sigma;
    FP const gphiphi = sin2 * ((r * r + hole.a * hole.a) * (r * r + hole.a * hole.a)
        - hole.a * hole.a * delta * sin2) / sigma;

    FP const orbital_term = hole.rs * r / FP(2) - hole.Q * hole.Q;
    if (orbital_term <= FP(0) || -gtt <= FP(0)) return FP(1);
    FP const root = sqrt(orbital_term);
    FP const omega = root / (r * r + hole.a * root); // prograde circular orbit
    FP const emitter_norm = -(gtt + FP(2) * gtphi * omega + gphiphi * omega * omega);
    if (emitter_norm <= FP(0)) return FP(1);

    FP const p_t = gtt * v[0] + gtphi * v[3];
    FP const p_phi = gtphi * v[0] + gphiphi * v[3];

    // The camera is static in these coordinates.  Its frequency measurement
    // is evaluated at the actual camera location, not approximated at infinity.
    FP const r0 = hole.r_0;
    FP const sigma0 = r0 * r0 + hole.a * hole.a * cos(hole.theta_0) * cos(hole.theta_0);
    FP const gtt0 = -(FP(1) - (hole.rs * r0 - hole.Q * hole.Q) / sigma0);
    if (-gtt0 <= FP(0)) return FP(1);

    FP const nu_camera = -p_t / sqrt(-gtt0);
    FP const nu_emitter = -(p_t + omega * p_phi) / sqrt(emitter_norm);
    if (nu_emitter == FP(0)) return FP(1);
    FP const shift = fabs(nu_camera / nu_emitter);
    return fmin(FP(5), fmax(FP(0.05), shift));
}


//ha a gömbön belül van akkor igaz
template <class FP>
inline bool gomb_be(FP const sugar, FP const* const x)
{
    if (x[1] < sugar)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//ha a gömbön kívül van akkor igaz
template <class FP>
inline bool gomb_ki(FP const sugar, FP const* const x)
{
    if (x[1] > sugar)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//ha x1 és x2 között van a disk akkor igaz disk síkja minkowskiban van és zy síkban
template <class FP>
inline bool disk(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2)
{
    if ((x1[1] > sugar_kicsi) && (x1[1] < sugar_nagy))
    {
        if (x1[2] > asin(1.0) && x2[2] < asin(1.0))
        {
            return true;
        }
        else if (x1[2] < asin(1.0) && x2[2] > asin(1.0))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

template <class FP>
inline bool disk1(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2)
{
    if ((x1[1] > sugar_kicsi) && (x1[1] < sugar_nagy))
    {
        if (x1[2] > asin(1.0) && x2[2] < asin(1.0))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

template <class FP>
inline bool disk2(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2)
{
    if ((x1[1] > sugar_kicsi) && (x1[1] < sugar_nagy))
    {
        if (x1[2] < asin(1.0) && x2[2] > asin(1.0))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

template <class FP>
inline int ijk_to_n(uint64_t const i, uint64_t const j, uint64_t const k, kerr_black_hole<FP> const& hole)
{
    return i * hole.MAGAS * D + j * D + k;
}

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
