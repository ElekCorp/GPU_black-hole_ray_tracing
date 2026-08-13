#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <math.h>

#include "black_hole.h"
#include "pown.h"

// The geodesic equation in canonical (Hamiltonian) form, as the symplectic
// integrators in symplectic.h need it.
//
//     H(x, p) = 1/2 g^{mu nu}(x) p_mu p_nu
//
//     dx^mu / dlambda =  dH/dp_mu = g^{mu nu} p_nu
//     dp_mu / dlambda = -dH/dx^mu = -1/2 (d_mu g^{alpha beta}) p_alpha p_beta
//
// This is the same geodesic, in the same affine parameter, as the second-order
// form cuda_ray.h integrates (x'' = -Gamma^i_jk v^j v^k, with v^mu = dx^mu/
// dlambda).  The difference is only which variables carry the state, and it
// matters for two reasons:
//
//   * (x, p) is a canonical pair, so a symplectic method really is symplectic
//     on it.  On (x, v) the pulled-back symplectic form is g_{mu nu}(x) dx ^
//     dv, which is not canonical, and the same method is only symmetric.
//   * t and phi are cyclic in Kerr-Newman, so dp[0] and dp[3] are identically
//     zero.  Photon energy and axial angular momentum are conserved to the
//     last bit rather than to within the integrator's tolerance - see
//     conserved_energy / conserved_angular_momentum below.
//
// H itself is the null constraint: g^{mu nu} p_mu p_nu = 0 is what makes the
// worldline a photon, so |H| is a directly meaningful error measure that costs
// nothing to evaluate (hamiltonian_value below).
//
// The generated block is emitted by tools/gen_hamiltonian.py, which imports
// the metric from tools/gen_christoffel.py so that both generated headers
// derive from one definition of Kerr-Newman.  Regenerate with
//
//     pixi run hamiltonian --write
//
// and check it with --verify (numeric, including a cross-check against the
// christoffel generator) and --check (drift against this file).

template <class FP>
inline void hamiltonian_rhs(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ p, FP* const __restrict__ dx, FP* const __restrict__ dp);
// v^mu = g^{mu nu} p_nu.  The dx half of hamiltonian_rhs on its own, for the
// conversions at the ends of an integration.
template <class FP>
inline void raise_momentum(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ p, FP* const __restrict__ v);
// p_mu = g_{mu nu} v^nu.  The entry point from the renderer's (x, v) state.
template <class FP>
inline void lower_velocity(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ p);


// H = 1/2 g^{mu nu} p_mu p_nu = 1/2 p_mu v^mu, i.e. the null constraint: 0 for
// a photon, and (with the (+,-,-,-) signature the generators use) +1/2 m^2 for
// a massive particle whose momentum is normalised to unit mass.  Conserved
// exactly by the flow, so its numerical drift is a pure integration error.
template <class FP>
inline FP hamiltonian_value(kerr_black_hole<FP> const& hole, FP const* const x, FP const* const p)
{
    FP v[D];
    raise_momentum(hole, x, p, v);
    FP h = FP(0);
    for (int n = 0; n < D; ++n) h += p[n] * v[n];
    return FP(0.5) * h;
}

// The two Killing-vector constants.  Both are components of p, so in these
// variables they are conserved by construction rather than approximately; they
// are here to be *measured* against the (x, v) integration, which conserves
// them only to its tolerance.
template <class FP>
inline FP conserved_energy(FP const* const p) { return p[0]; }
template <class FP>
inline FP conserved_angular_momentum(FP const* const p) { return -p[3]; }

// Carter's constant, the hidden invariant that makes Kerr geodesics integrable
// and the only one of the four that no integrator gets for free:
//
//     C = p_theta^2 + cos^2(theta) [ L_z^2 / sin^2(theta) - a^2 (E^2 - m^2) ]
//
// Only squares of E and L_z appear, so it is insensitive to the overall sign
// convention of the metric.  m^2 is read off the Hamiltonian rather than
// assumed, so this is valid for the timelike test orbits too, not just photons.
template <class FP>
inline FP carter_constant(kerr_black_hole<FP> const& hole, FP const* const x, FP const* const p)
{
    FP const sin_theta = sin(x[2]);
    FP const cos_theta = cos(x[2]);
    FP const s2 = fmax(sin_theta * sin_theta, FP(1e-12));
    FP const m2 = FP(2) * hamiltonian_value(hole, x, p);
    return p[2] * p[2] + cos_theta * cos_theta
        * (p[3] * p[3] / s2 - hole.a * hole.a * (p[0] * p[0] - m2));
}


// --- BEGIN GENERATED hamiltonian: tools/gen_hamiltonian.py ---


template <class FP>
inline void hamiltonian_rhs_general(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ p, FP* const __restrict__ dx, FP* const __restrict__ dp)
{
    FP const a = hole.a;
    FP const Q = hole.Q;
    FP const rs = hole.rs;

    FP const sy = sin(x[2]);
    FP const cy = cos(x[2]);

    // g^{phi phi} and its theta-derivative carry 1/sin^2 and 1/sin^3(theta).  Like
    // the 1/sin in christoffel(), this is a coordinate artifact of the Boyer-
    // Lindquist chart at the polar axis, not a curvature singularity: every term
    // carrying it is proportional to p[3] = p_phi, which vanishes on the axis
    // along with sin(theta), so the true ratio stays finite.  Flooring |sin| in
    // the reciprocal - and only there, sin as a factor is left exact - regularizes
    // the chart without perturbing any ray not already on the axis.
    FP const sy_pole_safe = (fabs(sy) > FP(1e-6)) ? sy : copysign(FP(1e-6), sy);
    FP const isy = FP(1) / sy_pole_safe;

    FP const x0 = pown(x[1], 2);
    FP const x1 = pown(a, 2);
    FP const x2 = x0 + x1;
    FP const x3 = pown(x2, 2);
    FP const x4 = pown(Q, 2) - rs*x[1];
    FP const x5 = x2 + x4;
    FP const x6 = pown(sy, 2)*x1;
    FP const x7 = x3 - x5*x6;
    FP const x8 = a*x4;
    FP const x9 = (FP(1)/(pown(cy, 2)*x1 + x0));
    FP const x10 = (FP(1)/(x5));
    FP const x11 = x10*x9;
    FP const x12 = x5*x9;
    FP const x13 = pown(isy, 2);
    FP const x14 = x5 - x6;
    FP const x15 = x13*x14;
    FP const x16 = pown(p[1], 2);
    FP const x17 = 2*x[1];
    FP const x18 = -x5;
    FP const x19 = x17*x9;
    FP const x20 = pown(p[2], 2)*x9;
    FP const x21 = pown(p[0], 2);
    FP const x22 = rs - x17;
    FP const x23 = x10*x22;
    FP const x24 = (1.0/2.0)*x10;
    FP const x25 = pown(p[3], 2);
    FP const x26 = p[0]*p[3];
    FP const x27 = sy*x1;

    dx[0] = -x11*(-p[0]*x7 + p[3]*x8);
    dx[1] = -p[1]*x12;
    dx[2] = -p[2]*x9;
    dx[3] = -x11*(p[0]*x8 + p[3]*x15);
    dp[0] = 0;
    dp[1] = -x9*(a*x10*x26*(rs + x19*x4 - x23*x4) + x13*x24*x25*(x14*x17*x9 - x14*x23 + x22) - 1.0/2.0*x16*(-rs + x17 + x18*x19) + x20*x[1] + x21*x24*(-x19*x7 + 4*x2*x[1] + x22*x6 + x23*x7));
    dp[2] = cy*x9*(2*pown(a, 3)*sy*x11*x26*x4 + isy*x10*x25*(x1*x14*x9 - x1 - x15) + x12*x16*x27 + x20*x27 + x21*x27*(-x11*(x18*x6 + x3) + 1));
    dp[3] = 0;
}


template <class FP>
inline void hamiltonian_rhs_static(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ p, FP* const __restrict__ dx, FP* const __restrict__ dp)
{
    FP const Q = hole.Q;
    FP const rs = hole.rs;

    FP const sy = sin(x[2]);
    FP const cy = cos(x[2]);

    // g^{phi phi} and its theta-derivative carry 1/sin^2 and 1/sin^3(theta).  Like
    // the 1/sin in christoffel(), this is a coordinate artifact of the Boyer-
    // Lindquist chart at the polar axis, not a curvature singularity: every term
    // carrying it is proportional to p[3] = p_phi, which vanishes on the axis
    // along with sin(theta), so the true ratio stays finite.  Flooring |sin| in
    // the reciprocal - and only there, sin as a factor is left exact - regularizes
    // the chart without perturbing any ray not already on the axis.
    FP const sy_pole_safe = (fabs(sy) > FP(1e-6)) ? sy : copysign(FP(1e-6), sy);
    FP const isy = FP(1) / sy_pole_safe;

    FP const x0 = pown(x[1], 2);
    FP const x1 = pown(Q, 2) - rs*x[1] + x0;
    FP const x2 = (FP(1)/(x1));
    FP const x3 = (FP(1)/(x0));
    FP const x4 = pown(isy, 2);
    FP const x5 = (FP(1)/pown(x[1], 3));
    FP const x6 = pown(p[3], 2);
    FP const x7 = rs - 2*x[1];
    FP const x8 = x2*x[1];

    dx[0] = p[0]*x0*x2;
    dx[1] = -p[1]*x1*x3;
    dx[2] = -p[2]*x3;
    dx[3] = -p[3]*x3*x4;
    dp[0] = 0;
    dp[1] = -1.0/2.0*pown(p[0], 2)*x8*(x7*x8 + 2) - 1.0/2.0*pown(p[1], 2)*x3*(2*x1/x[1] + x7) - pown(p[2], 2)*x5 - x4*x5*x6;
    dp[2] = -cy*pown(isy, 3)*x3*x6;
    dp[3] = 0;
}


template <class FP>
inline void hamiltonian_rhs(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ p, FP* const __restrict__ dx, FP* const __restrict__ dp)
{
    // a=0 collapses Kerr-Newman to the spherically symmetric, non-frame-
    // dragging case, which is a much shorter block - the same split, and the
    // same warp-uniform branch, as christoffel() in cuda_ray.h.
    if (hole.a == FP(0))
    {
        hamiltonian_rhs_static(hole, x, p, dx, dp);
    }
    else
    {
        hamiltonian_rhs_general(hole, x, p, dx, dp);
    }
}


template <class FP>
inline void raise_momentum_general(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ p, FP* const __restrict__ v)
{
    FP const a = hole.a;
    FP const Q = hole.Q;
    FP const rs = hole.rs;

    FP const sy = sin(x[2]);
    FP const cy = cos(x[2]);

    // g^{phi phi} and its theta-derivative carry 1/sin^2 and 1/sin^3(theta).  Like
    // the 1/sin in christoffel(), this is a coordinate artifact of the Boyer-
    // Lindquist chart at the polar axis, not a curvature singularity: every term
    // carrying it is proportional to p[3] = p_phi, which vanishes on the axis
    // along with sin(theta), so the true ratio stays finite.  Flooring |sin| in
    // the reciprocal - and only there, sin as a factor is left exact - regularizes
    // the chart without perturbing any ray not already on the axis.
    FP const sy_pole_safe = (fabs(sy) > FP(1e-6)) ? sy : copysign(FP(1e-6), sy);
    FP const isy = FP(1) / sy_pole_safe;

    FP const x0 = pown(x[1], 2);
    FP const x1 = pown(a, 2);
    FP const x2 = x0 + x1;
    FP const x3 = pown(Q, 2) - rs*x[1];
    FP const x4 = x2 + x3;
    FP const x5 = pown(sy, 2)*x1;
    FP const x6 = a*x3;
    FP const x7 = (FP(1)/(pown(cy, 2)*x1 + x0));
    FP const x8 = x7/x4;

    v[0] = -x8*(-p[0]*(pown(x2, 2) - x4*x5) + p[3]*x6);
    v[1] = -p[1]*x4*x7;
    v[2] = -p[2]*x7;
    v[3] = -x8*(pown(isy, 2)*p[3]*(x4 - x5) + p[0]*x6);
}


template <class FP>
inline void raise_momentum_static(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ p, FP* const __restrict__ v)
{
    FP const Q = hole.Q;
    FP const rs = hole.rs;

    FP const sy = sin(x[2]);

    // g^{phi phi} and its theta-derivative carry 1/sin^2 and 1/sin^3(theta).  Like
    // the 1/sin in christoffel(), this is a coordinate artifact of the Boyer-
    // Lindquist chart at the polar axis, not a curvature singularity: every term
    // carrying it is proportional to p[3] = p_phi, which vanishes on the axis
    // along with sin(theta), so the true ratio stays finite.  Flooring |sin| in
    // the reciprocal - and only there, sin as a factor is left exact - regularizes
    // the chart without perturbing any ray not already on the axis.
    FP const sy_pole_safe = (fabs(sy) > FP(1e-6)) ? sy : copysign(FP(1e-6), sy);
    FP const isy = FP(1) / sy_pole_safe;

    FP const x0 = pown(x[1], 2);
    FP const x1 = pown(Q, 2) - rs*x[1] + x0;
    FP const x2 = (FP(1)/(x0));

    v[0] = p[0]*x0/x1;
    v[1] = -p[1]*x1*x2;
    v[2] = -p[2]*x2;
    v[3] = -pown(isy, 2)*p[3]*x2;
}


template <class FP>
inline void raise_momentum(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ p, FP* const __restrict__ v)
{
    // a=0 collapses Kerr-Newman to the spherically symmetric, non-frame-
    // dragging case, which is a much shorter block - the same split, and the
    // same warp-uniform branch, as christoffel() in cuda_ray.h.
    if (hole.a == FP(0))
    {
        raise_momentum_static(hole, x, p, v);
    }
    else
    {
        raise_momentum_general(hole, x, p, v);
    }
}


template <class FP>
inline void lower_velocity_general(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ p)
{
    FP const a = hole.a;
    FP const Q = hole.Q;
    FP const rs = hole.rs;

    FP const sy = sin(x[2]);
    FP const cy = cos(x[2]);

    FP const x0 = pown(x[1], 2);
    FP const x1 = pown(a, 2);
    FP const x2 = pown(cy, 2)*x1 + x0;
    FP const x3 = (FP(1)/(x2));
    FP const x4 = pown(sy, 2);
    FP const x5 = x1*x4;
    FP const x6 = pown(Q, 2) - rs*x[1];
    FP const x7 = x0 + x1 + x6;
    FP const x8 = a*x6;

    p[0] = -x3*(-v[0]*(-x5 + x7) + v[3]*x4*x8);
    p[1] = -v[1]*x2/x7;
    p[2] = -v[2]*x2;
    p[3] = -x3*x4*(v[0]*x8 + v[3]*(pown(a, 4) + 2*x0*x1 - x5*x7 + pown(x[1], 4)));
}


template <class FP>
inline void lower_velocity_static(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ p)
{
    FP const Q = hole.Q;
    FP const rs = hole.rs;

    FP const sy = sin(x[2]);

    FP const x0 = pown(x[1], 2);
    FP const x1 = pown(Q, 2) - rs*x[1] + x0;

    p[0] = v[0]*x1/x0;
    p[1] = -v[1]*x0/x1;
    p[2] = -v[2]*x0;
    p[3] = -pown(sy, 2)*v[3]*x0;
}


template <class FP>
inline void lower_velocity(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ p)
{
    // a=0 collapses Kerr-Newman to the spherically symmetric, non-frame-
    // dragging case, which is a much shorter block - the same split, and the
    // same warp-uniform branch, as christoffel() in cuda_ray.h.
    if (hole.a == FP(0))
    {
        lower_velocity_static(hole, x, v, p);
    }
    else
    {
        lower_velocity_general(hole, x, v, p);
    }
}


// --- END GENERATED hamiltonian ---

#endif // HAMILTONIAN_H
