#ifndef CUDA_RAY_H
#define CUDA_RAY_H

#include <iostream>
#include <math.h>

#include "black_hole.h"

//#include "debugmalloc.h"


// Floats ray_step_T writes per pixel.  Channel 0 keeps the original hit code
// (-1 escaped, 0 captured, r > 0 an accretion-disk hit); the rest record
// quantities the integration already produced, so the shading layer can be
// physical without re-deriving anything or changing a single geodesic:
//
//   0  hit radius / escape code
//   1  disk redshift g (disk hits only)
//   2  disk hit Boyer-Lindquist azimuth phi (disk hits only)
//   3  coordinate time elapsed along the ray, i.e. the light travel time from
//      the emission event to the camera.  The shading layer subtracts it to
//      evaluate the disk at each pixel's own retarded time, so the strongly
//      lensed images (which took a longer path) correctly show the disk as it
//      was earlier than the direct image does.
//   4  escape direction polar angle on the sky (escaped rays only)
//   5  escape direction azimuth on the sky (escaped rays only)
//
// Channels 4/5 let the background sky be sampled as a texture on the celestial
// sphere along the ray's own outgoing direction, so the star field is lensed by
// the same geodesics as everything else instead of being pasted on in screen
// space.
#define RAY_CHANNELS 6


template <class FP>
inline void step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de, FP* const __restrict__ deriv_x, FP* const __restrict__ deriv_v, bool& have_deriv);
template <class FP>
inline void dopri54_step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de, FP* const __restrict__ deriv_x, FP* const __restrict__ deriv_v, bool& have_deriv);

template <class FP>
inline FP ijk_to_vec_mink_zoom(uint64_t const i, uint64_t const j, uint64_t const k, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg, kerr_black_hole<FP> const& hole);
template <class FP>
inline FP ijk_to_vec_zoom(uint64_t const i, uint64_t const j, uint64_t const k, kerr_black_hole<FP> const& hole, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg);



template <class FP>
inline void RK38(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de);//4-ed foku legpontosabb
template <class FP>
inline void RK6(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de);

template <class FP>
inline void christoffel_general(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ ch);
template <class FP>
inline void christoffel_static(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ ch);
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

// Direction the ray is travelling when it leaves the outer sphere, expressed
// as (polar, azimuth) on the celestial sphere.  See RAY_CHANNELS.
template <class FP>
inline void sky_direction(FP const* const x, FP const* const v, FP& sky_theta, FP& sky_phi);

// Estimate of the state a fraction `frac` of the way from the previous step to
// the current one.  Used to recover the moment a boundary was actually reached
// from the two steps that straddle it.
template <class FP>
inline void interpolate_state(FP const frac, FP const* const x_le, FP const* const v_le,
                              FP const* const x, FP const* const v, FP* const x_hit, FP* const v_hit);

// Estimate of the state at the moment the ray actually crossed the disk plane,
// interpolated between the last step that was on one side and the first that
// was on the other.
template <class FP>
inline void disk_crossing(FP const* const x_le, FP const* const v_le, FP const* const x, FP const* const v,
                          FP* const x_hit, FP* const v_hit);

// The same, for the moment the ray reached the outer sphere at radius `sugar`.
template <class FP>
inline void sphere_crossing(FP const sugar, FP const* const x_le, FP const* const v_le,
                            FP const* const x, FP const* const v, FP* const x_hit, FP* const v_hit);




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
inline void step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de, FP* const __restrict__ deriv_x, FP* const __restrict__ deriv_v, bool& have_deriv)//adaptiv step size
{
    dopri54_step(hole, x, v, de, deriv_x, deriv_v, have_deriv);
}

// Embedded Dormand-Prince 5(4) controller for the first-order form of the
// geodesic equation: x' = v, v' = Gamma(x, v).  Unlike the earlier
// curvature-magnitude heuristic, the 5th-vs-4th-order difference is a local
// truncation-error estimate.  It also removes the old post-step Christoffel
// evaluation (7 evaluations per accepted attempt instead of 8).
//
// Stage 1 (k[0]) only depends on x,v, which are unchanged across rejected
// attempts within a step and across the FSAL boundary between consecutive
// accepted steps (c7=1, so the last stage of an accepted step is already the
// derivative at the new x,v).  deriv_x/deriv_v/have_deriv let the caller
// carry that derivative in and out, so it is computed at most once per
// accepted step instead of once per attempt.
template <class FP>
inline void dopri54_step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de, FP* const __restrict__ deriv_x, FP* const __restrict__ deriv_v, bool& have_deriv)
{
    FP const max_step = hole.de0;
    FP const min_step = fmax(FP(1e-8), max_step * FP(1e-5));
    FP const tolerance = fmax(hole.errormax, FP(1e-8));
    FP h = fmin(max_step, fmax(min_step, de));
    FP kx[7][D], kv[7][D], trial_x[D], trial_v[D], fifth_x[D], fifth_v[D], fourth_x[D], fourth_v[D], acceleration[D];

    if (have_deriv)
    {
        for (int n = 0; n < D; ++n) { kx[0][n] = deriv_x[n]; kv[0][n] = deriv_v[n]; }
    }
    else
    {
        christoffel(hole, x, v, acceleration);
        for (int n = 0; n < D; ++n) { kx[0][n] = v[n]; kv[0][n] = acceleration[n]; }
    }

    for (int attempt = 0; attempt < 10; ++attempt)
    {
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

        // Keep the controller in direct step-size arithmetic.  NVHPC 26.3's
        // OpenACC LLVM lowering can emit an undefined temporary for a local
        // "factor" variable in this inlined routine.
        FP next_h = h * FP(5);
        if (error_norm > FP(0))
        {
            next_h = h * FP(0.9) * pow(error_norm, FP(-0.2));
        }
        next_h = fmin(h * FP(5), fmax(h * FP(0.2), next_h));
        next_h = fmin(max_step, fmax(min_step, next_h));
        if (error_norm <= FP(1) || h <= min_step || attempt == 9)
        {
            for (int n = 0; n < D; ++n) { x[n] = fifth_x[n]; v[n] = fifth_v[n]; }
            de = next_h;
            for (int n = 0; n < D; ++n) { deriv_x[n] = kx[6][n]; deriv_v[n] = kv[6][n]; }
            have_deriv = true;
            return;
        }
        h = next_h;
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

// --- BEGIN GENERATED christoffel: tools/gen_christoffel.py ---
// christoffel_general: temps=52 add=50 mul=160 div=4 trig=2 total=216
// christoffel_static: temps=10 add=10 mul=36 div=4 trig=2 total=52
template <class FP>
inline void christoffel_general(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ ch)
{
    FP a = hole.a;
    FP Q = hole.Q;
    FP rs = hole.rs;

    FP y = x[2];

    FP x0 = pown(x[1], 2);
    FP x1 = pown(a, 2);
    FP x2 = x0 + x1;
    FP x3 = pown(Q, 2) - rs*x[1];
    FP x4 = x2 + x3;
    FP x5 = sin(y);
    FP x6 = pown(x5, 2);
    FP x7 = x1*x6;
    FP x8 = x4*x7;
    FP x9 = -pown(x2, 2) + x8;
    FP x10 = 2*x[1];
    FP x11 = cos(y);
    FP x12 = x0 + x1*pown(x11, 2);
    FP x13 = pown(x12, -1);
    FP x14 = x13*x4;
    FP x15 = x13*x7;
    FP x16 = -rs - x10*x14 + x10*x15 + x10;
    FP x17 = x10*x13;
    FP x18 = rs + x17*x3;
    FP x19 = x18*x3*x7;
    FP x20 = v[0]*v[1];
    FP x21 = rs - x10;
    FP x22 = 2*x1;
    FP x23 = pown(a, 4) + x0*x22 + pown(x[1], 4);
    FP x24 = x23 - x8;
    FP x25 = 4*x1*x[1] + 4*pown(x[1], 3);
    FP x26 = -x17*x24 + x21*x7 + x25;
    FP x27 = v[1]*v[3];
    FP x28 = a*x6;
    FP x29 = x15 + 1;
    FP x30 = x29*pown(x3, 2);
    FP x31 = -x14 + x29;
    FP x32 = x22*x5;
    FP x33 = v[2]*x11;
    FP x34 = x23 - 2*x8;
    FP x35 = x15*x24 + x34;
    FP x36 = 2*v[3];
    FP x37 = a*x3;
    FP x38 = pown(x4, -1);
    FP x39 = x38/pown(x12, 2);
    FP x40 = pown(v[2], 2);
    FP x41 = pown(v[0], 2);
    FP x42 = (1.0/2.0)*x13;
    FP x43 = pown(v[1], 2);
    FP x44 = pown(v[3], 2);
    FP x45 = x1*x5;
    FP x46 = -x24;
    FP x47 = x11*x13;
    FP x48 = x47*x5;
    FP x49 = x31*x45;
    FP x50 = x4 - x7;

    // The phi equation carries an explicit 1/sin(theta).  This is a coordinate
    // artifact of Boyer-Lindquist-type charts at the polar axis (theta=0,pi), not a
    // physical curvature singularity: every numerator multiplying it is itself
    // proportional to v[3]=dphi/de, which vanishes on the axis along with
    // sin(theta), so the true ratio stays finite.  A ray integrated exactly through
    // the axis would otherwise hit a 0/0-like blow-up.  Flooring |sin(theta)|
    // regularizes the chart without perturbing any ray not already on the axis (a
    // set of measure zero).
    FP const x5_pole_safe = (fabs(x5) > FP(1e-6)) ? x5 : copysign(FP(1e-6), x5);

    FP x51 = x50/x5_pole_safe;

    ch[0] = -x39*(v[0]*x32*x33*(x30 + x31*x9) - x20*(x16*x9 + x19) + x27*x28*(-x18*x9 + x26*x3) + x33*x36*x37*x5*(x29*x9 + x35));
    ch[1] = x14*(-v[0]*v[3]*x13*x28*(rs + x17*x3) + v[1]*v[2]*x11*x32*x38 - x16*x41*x42 + x26*x42*x44*x6 - 1.0/2.0*x38*x43*(x10 + x12*x21*x38) + x40*x[1]);
    ch[2] = -x13*(-v[0]*x29*x36*x37*x48 + 2*v[1]*v[2]*x[1] + x1*x11*x38*x43*x5 - x11*x40*x45 - x41*x47*x49 - x44*x48*(-x15*x46 + x34));
    ch[3] = -x39*(2*a*v[0]*v[2]*x11*x3*(x29*x51 + x49) - a*x20*(x16*x3 + x18*x50) + 2*v[2]*v[3]*x11*(x30*x45 + x35*x51) - x27*(x19 - x50*(x17*x46 + x21*x7 + x25)));
}


template <class FP>
inline void christoffel_static(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ ch)
{
    FP Q = hole.Q;
    FP rs = hole.rs;

    FP y = x[2];

    FP x0 = pown(Q, 2) - rs*x[1] + pown(x[1], 2);
    FP x1 = pown(x0, -1);
    FP x2 = pown(x[1], -1);
    FP x3 = x0*x2;
    FP x4 = rs - 2*x[1];
    FP x5 = 2*x3 + x4;
    FP x6 = pown(v[3], 2);
    FP x7 = sin(y);
    FP x8 = cos(y);
    FP x9 = v[1]*x2;

    ch[0] = v[0]*v[1]*x1*x5;
    ch[1] = x3*((1.0/2.0)*pown(v[0], 2)*x5/pown(x[1], 3) - 1.0/2.0*pown(v[1], 2)*x1*(x1*x4*x[1] + 2) + pown(v[2], 2) + x6*pown(x7, 2));
    ch[2] = -2*v[2]*x9 + x6*x7*x8;

    // The phi equation carries an explicit 1/sin(theta).  This is a coordinate
    // artifact of Boyer-Lindquist-type charts at the polar axis (theta=0,pi), not a
    // physical curvature singularity: every numerator multiplying it is itself
    // proportional to v[3]=dphi/de, which vanishes on the axis along with
    // sin(theta), so the true ratio stays finite.  A ray integrated exactly through
    // the axis would otherwise hit a 0/0-like blow-up.  Flooring |sin(theta)|
    // regularizes the chart without perturbing any ray not already on the axis (a
    // set of measure zero).
    FP const x7_pole_safe = (fabs(x7) > FP(1e-6)) ? x7 : copysign(FP(1e-6), x7);

    ch[3] = -2*v[3]*(v[2]*x8/x7_pole_safe + x9);
}


template <class FP>
inline void christoffel(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ ch)
{
    // a=0 collapses Kerr-Newman to the spherically symmetric, non-frame-
    // dragging case (Reissner-Nordstrom, or Schwarzschild when Q is also 0):
    // no t/phi coupling and no theta-dependence in rho2, so christoffel_static
    // is ~4x fewer ops than christoffel_general (52 vs 216). A Q=0-only
    // specialization was measured and rejected - Q never couples coordinates,
    // so dropping it alone saves nothing. hole.a is the same value for every
    // thread in a render, so this branch is warp-uniform: no GPU divergence.
    if (hole.a == FP(0))
    {
        christoffel_static(hole, x, v, ch);
    }
    else
    {
        christoffel_general(hole, x, v, ch);
    }
}
// --- END GENERATED christoffel ---

// Superseded reference form (sympy simplify() instead of factor_terms -
// see tools/gen_christoffel.py for why that pipeline was dropped).
    /*FP a = hole.a;
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
*/
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

        // Sample at pixel centers (i+0.5), not pixel corners (i).  Corner
        // sampling puts one whole column/row of rays exactly on the camera's
        // meridian-symmetry axis (z=0 or y=0), where the photon's angular
        // momentum about that axis is exactly zero.  Those rays are a
        // measure-zero degenerate family that passes through the coordinate
        // pole and can hit/miss the disk differently than their immediate
        // neighbours, showing up as a visible hairline through the image
        // centre.  Centering the sample removes any ray from that axis.
        //
        // Deliberately double, even when FP is float: this is a difference of
        // two indices that are both of order SZELESregi, so its relative
        // precision is what limits how far the zoom window can be narrowed
        // before neighbouring pixels start collapsing onto the same ray.  In
        // float that ceiling is ~2^24/SZELES, i.e. a few thousand times
        // magnification; in double it is far beyond anything the geodesic
        // integrator itself can resolve.  The cost is a handful of double
        // operations per ray, against thousands of integration steps.
        double ir = double(ikezd) + ((double(i) + 0.5) / double(SZELES)) * (double(iveg) - double(ikezd));
        double jr = double(jkezd) + ((double(j) + 0.5) / double(SZELES)) * (double(iveg) - double(ikezd));//igen igy jo csak bele kell gondolni SZELES/MAGAS=(iveg-ikezd)/(jveg-jkezd)

        double y = (kemernyo_high / double(MAGASregi)) * (double(MAGASregi) / 2 - jr);
        double z = (kemernyo_high / double(MAGASregi)) * (ir - double(SZELESregi) / 2);

        double norm = sqrt(x * x + y * y + z * z);


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
    // The sin(phi) terms belong outside the (1 - cos(phi)) factor, and the
    // x[3] coefficient is a minus.  Folded in the way they were, this row was
    // only a rotation when u[0]*sin(phi) and u[2]*sin(phi) both vanished; at
    // 1.2 rad about a general axis the matrix had det = 0.444 and
    // |R R^T - I| = 0.54, so it sheared the camera basis rather than rotating
    // it, and the ray directions it produced were not unit vectors.
    //
    // That degeneracy is why this survived: web_server.py starts from a
    // 180-degree flip about the local up axis, and both of the ways out of it
    // stay in the safe set - panning composes to another rotation about up
    // (u[0] = u[2] = 0), while a tilt or roll on its own tips the axis but
    // leaves the angle at exactly pi (sin(phi) = 0).  Only a pan combined with
    // a tilt or roll leaves both, which is where the sheared frames came from.
    // See tests/test_camera_rotation.cpp.
    FP x2_tmp = (u[1] * u[0] * (1 - cos(phi)) + u[2] * sin(phi)) * x[1] + (cos(phi) + u[1] * u[1] * (1 - cos(phi))) * x[2] + (u[1] * u[2] * (1 - cos(phi)) - u[0] * sin(phi)) * x[3];
    FP x3_tmp = (u[0] * u[2] * (1 - cos(phi)) - u[1] * sin(phi)) * x[1] + (u[1] * u[2] * (1 - cos(phi)) + u[0] * sin(phi)) * x[2] + (cos(phi) + u[2] * u[2] * (1 - cos(phi))) * x[3];

    x[1]=x1_tmp;
    x[2]=x2_tmp;
    x[3]=x3_tmp;

    FP r_0 = hole.r_0;
    FP theta_0 = hole.theta_0;
    FP rs = hole.rs;
    FP a = hole.a;
    FP Q = hole.Q;

    // Standard Kerr-Newman Delta = r^2 - rs*r + a^2 + Q^2 (same convention as
    // christoffel's `Q^2 - rs*r` term and the horizon radius formula in
    // ray_step). A stray factor of 4 here made this go negative - and the
    // sqrt() below NaN - well outside the true horizon (e.g. r0=0.2 instead
    // of ~0.05 for the default rs=0.05, a=Q=0), silently turning every
    // camera ray into a NaN direction and rendering as a blank starfield.
    FP delta = r_0 * r_0 - rs * r_0 + a * a + Q * Q;
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

#if defined(USE_OPENACC)
// Omega is D-1 = 3 elements, not 4 - see its declarations in ffi_bridge.cpp and
// main.cpp, and kerr_black_hole's constructor, which reads only Omega[0..2].
#pragma acc data copyin(xd[0:4],Omega[0:3]) copyout(szin[0:SZELES*MAGAS])
{
#pragma acc parallel loop collapse(2)
#elif defined(USE_OPENMP)
#pragma omp parallel for collapse(2)
#endif
    for(uint64_t j=0; j<MAGAS; j++)
    {
        for(uint64_t i=0; i<SZELES; i++)
        {
            FP x[D] = { hole.t_0,hole.r_0,hole.theta_0,hole.phi_0 };;
            FP v[D];
            FP x_le[D]; //lemaradó hely koordináták
            FP de = de0;
            FP deriv_x[D], deriv_v[D];
            bool have_deriv = false;

            for (int k = 0; k < D; ++k)
            {
                v[k] = ijk_to_vec_zoom(i, j, k, hole, SZELESregi, MAGASregi, ikezd, jkezd, iveg);
            }




            for (int k = 0; k < D; ++k)
            {
                x_le[k] = x[k];
            }
            step(hole, x, v, de, deriv_x, deriv_v, have_deriv);

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
                step(hole, x, v, de, deriv_x, deriv_v, have_deriv);




            }



            //printf("%d\t%d\t%f\n", i, j, de);
        }
    }
#if defined(USE_OPENACC)
}
#endif
}

template <class FP>
inline void ray_step_T(FP* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const xd, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki_in, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg)//kernel
{
    kerr_black_hole<FP> hole(SZELES, MAGAS, xd, Omega, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki_in, gyuru_sugar_kicsi, gyuru_sugar_nagy);

#if defined(USE_OPENACC)
// Omega[0:3], not [0:4]. It is declared FP[D-1] = 3 elements by every caller,
// and the constructor reads only Omega[0..2], so claiming a fourth walks off
// the end of the array. That is not a harmless over-copy: the callers declare
// Omega and x as adjacent stack locals, and at FP=double the extra 8 bytes land
// exactly on top of xd, which OpenACC then rejects at runtime with
//
//   FATAL ERROR: variable in data clause is partially present on the device:
//   name=Omega[:4]
//
// because the region overlaps one already in the present table. At FP=float the
// layout happened not to collide, which is why only the f64 path ever failed -
// the f32 renders had been running against this the whole time.
#pragma acc data copyin(xd[0:4],Omega[0:3]) copyout(szin[0:RAY_CHANNELS*SZELES*MAGAS])
{
#pragma acc parallel loop collapse(2)
#elif defined(USE_OPENMP)
#pragma omp parallel for collapse(2)
#endif
    for(uint64_t j=0; j<MAGAS; j++)
    {
        for(uint64_t i=0; i<SZELES; i++)
        {

            FP x[D] = { hole.t_0,hole.r_0,hole.theta_0,hole.phi_0 };;
            FP v[D];
            FP x_le[D]; //lemaradó hely koordináták
            FP v_le[D]; // and the matching velocity, for disk_crossing
            FP de = de0;
            FP deriv_x[D], deriv_v[D];
            bool have_deriv = false;

            for (int k = 0; k < D; ++k)
            {
                v[k] = ijk_to_vec_zoom(i, j, k, hole, SZELESregi, MAGASregi, ikezd, jkezd, iveg);
            }




            for (int k = 0; k < D; ++k)
            {
                x_le[k] = x[k];
                v_le[k] = v[k];
            }
            step(hole, x, v, de, deriv_x, deriv_v, have_deriv);

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

            // All RAY_CHANNELS outputs for this pixel live contiguously here.
            FP* const px = szin + RAY_CHANNELS * (i * MAGAS + j);
            FP sky_theta, sky_phi;

            while (fut)
            {
                //0 fekete, -1 hiba kezeléses piros, egyébként meg egy FP ami reprezental egy szint
                if (gomb_be(sugar_be, x))
                {
                    px[0] = 0;//;
                    px[1] = 0;
                    px[2] = 0;
                    px[3] = x[0];
                    px[4] = 0;
                    px[5] = 0;
                    fut = false;
                }
                else if (gomb_ki(sugar_ki, x))
                {
                    // The ray reached the outer integration sphere: it
                    // escaped to the sky, rather than being captured.
                    FP x_hit[D], v_hit[D];
                    sphere_crossing(sugar_ki, x_le, v_le, x, v, x_hit, v_hit);
                    sky_direction(x_hit, v_hit, sky_theta, sky_phi);
                    px[0] = -1;
                    px[1] = 0;
                    px[2] = 0;
                    px[3] = x_hit[0];
                    px[4] = sky_theta;
                    px[5] = sky_phi;
                    fut = false;
                }
                else if (disk1(sugar_kicsi, sugar_nagy, x, x_le) || disk2(sugar_kicsi, sugar_nagy, x, x_le))
                {
                    FP x_hit[D], v_hit[D];
                    disk_crossing(x_le, v_le, x, v, x_hit, v_hit);
                    px[0] = x_hit[1];
                    px[1] = disk_redshift(hole, x_hit, v_hit);
                    px[2] = x_hit[3];
                    px[3] = x_hit[0];
                    px[4] = 0;
                    px[5] = 0;
                    fut = false;
                }
                else if (isnan(x[0]) || isnan(x[1]) || isnan(x[2]) || isnan(x[3]))
                {
                    //printf("%d\t%d\t%f\tnan\n", i, j, de);
                    // x has already gone non-finite, so the sky direction has
                    // to come from the last good state; sky_direction falls
                    // back to a fixed direction if even that is unusable.
                    sky_direction(x_le, v, sky_theta, sky_phi);
                    px[0] = -1;
                    px[1] = 0;
                    px[2] = 0;
                    px[3] = x_le[0];
                    px[4] = sky_theta;
                    px[5] = sky_phi;
                    fut = false;
                }
                else if (isinf(x[0]) || isinf(x[1]) || isinf(x[2]) || isinf(x[3]))
                {
                    //printf("%d\t%d\t%f\tinf\n", i, j, de);
                    sky_direction(x_le, v, sky_theta, sky_phi);
                    px[0] = -1;
                    px[1] = 0;
                    px[2] = 0;
                    px[3] = x_le[0];
                    px[4] = sky_theta;
                    px[5] = sky_phi;
                    fut = false;
                }
                else
                {

                }

                ++idokorlat;
                if (idokorlat >= int(1.0 / errormax))//if (idokorlat >= int(1.0 / errormax))
                {
                    //printf("%d\t%d\t%f\tmegunta\n", i, j, de);
                    sky_direction(x, v, sky_theta, sky_phi);
                    px[0] = -1;
                    px[1] = 0;
                    px[2] = 0;
                    px[3] = x[0];
                    px[4] = sky_theta;
                    px[5] = sky_phi;
                    fut = false;
                }


                for (int k = 0; k < D; ++k)
                {
                    x_le[k] = x[k];
                    v_le[k] = v[k];
                }
                step(hole, x, v, de, deriv_x, deriv_v, have_deriv);




            }



            //printf("%d\t%d\t%f\n", i, j, de);

        }
    }
#if defined(USE_OPENACC)
}
#endif
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


template <class FP>
inline void interpolate_state(FP const frac, FP const* const x_le, FP const* const v_le,
                              FP const* const x, FP const* const v, FP* const x_hit, FP* const v_hit)
{
    // The collision tests all fire on the first step that lands past the
    // boundary, so the state recorded there is up to one whole adaptive step
    // beyond where the ray actually reached it - and the step length varies
    // from ray to ray.  Reporting it unmodified makes the recorded hit jitter
    // by that much between neighbouring pixels, which is a real error in where
    // the ray landed, and which showed up downstream as moire banding across
    // the disk and a ragged background sky.  Interpolating between the two
    // bracketing states puts the reported point back on the boundary.  The
    // trajectory itself is not touched; only which point along it is reported.
    FP const clamped = fmin(FP(1), fmax(FP(0), frac));
    for (int k = 0; k < D; ++k)
    {
        x_hit[k] = x_le[k] + clamped * (x[k] - x_le[k]);
        v_hit[k] = v_le[k] + clamped * (v[k] - v_le[k]);
    }
}


template <class FP>
inline void disk_crossing(FP const* const x_le, FP const* const v_le, FP const* const x, FP const* const v,
                          FP* const x_hit, FP* const v_hit)
{
    FP const plane = asin(FP(1.0));
    FP const before = x_le[2] - plane;
    FP const span = before - (x[2] - plane);
    interpolate_state((span != FP(0)) ? before / span : FP(1), x_le, v_le, x, v, x_hit, v_hit);
    x_hit[2] = plane;
}


template <class FP>
inline void sphere_crossing(FP const sugar, FP const* const x_le, FP const* const v_le,
                            FP const* const x, FP const* const v, FP* const x_hit, FP* const v_hit)
{
    FP const span = x[1] - x_le[1];
    interpolate_state((span != FP(0)) ? (sugar - x_le[1]) / span : FP(1), x_le, v_le, x, v, x_hit, v_hit);
    x_hit[1] = sugar;
}


template <class FP>
inline void sky_direction(FP const* const x, FP const* const v, FP& sky_theta, FP& sky_phi)
{
    // Read at the outer sphere, where rs/r is small enough that the
    // Boyer-Lindquist coordinate basis is close to an orthonormal spherical
    // frame, so the outgoing direction can be assembled from the coordinate
    // velocity components directly.  Nothing about the geodesic changes here;
    // this only names the direction the ray was already travelling in.
    FP const r = x[1];
    FP const sin_theta = sin(x[2]);
    FP const cos_theta = cos(x[2]);
    FP const sin_phi = sin(x[3]);
    FP const cos_phi = cos(x[3]);

    FP const n_r = v[1];
    FP const n_theta = r * v[2];
    FP const n_phi = r * sin_theta * v[3];

    // Spherical unit vectors expanded onto a fixed Cartesian sky frame.
    FP const dx = n_r * sin_theta * cos_phi + n_theta * cos_theta * cos_phi - n_phi * sin_phi;
    FP const dy = n_r * sin_theta * sin_phi + n_theta * cos_theta * sin_phi + n_phi * cos_phi;
    FP const dz = n_r * cos_theta - n_theta * sin_theta;

    FP const norm = sqrt(dx * dx + dy * dy + dz * dz);
    if (!(norm > FP(0)))
    {
        // Only reachable for a ray that already went non-finite; pick a fixed
        // direction rather than propagating a NaN into the shading layer.
        sky_theta = FP(asin(1.0));
        sky_phi = 0;
        return;
    }
    sky_theta = acos(fmin(FP(1), fmax(FP(-1), dz / norm)));
    sky_phi = atan2(dy, dx);
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
