#ifndef SYMPLECTIC_H
#define SYMPLECTIC_H

#include <math.h>

#include "black_hole.h"
#include "hamiltonian.h"
#include "pown.h"

// Implicit, symplectic integrators for the geodesic equation, as an
// alternative to the explicit Dormand-Prince 5(4) in cuda_ray.h.
//
// The family is Gauss-Legendre collocation: the s-stage member is implicit,
// has order 2s, and is both symplectic and symmetric (Hairer, Lubich & Wanner,
// *Geometric Numerical Integration*, II.1.3 and VI.4.1).  Symplectic is the
// property worth paying for here.  A photon is defined by the null constraint
// H = 1/2 g^{mu nu} p_mu p_nu = 0, and H is exactly the Hamiltonian; a
// symplectic method run at constant step size does not conserve H exactly, but
// it conserves the Hamiltonian of a nearby *modified* system, which pins |H|
// to a bounded O(h^{2s}) oscillation for exponentially many steps instead of
// letting it walk away linearly in the step count.  Dormand-Prince has no such
// backward-error structure and its constraint error accumulates.
//
// Three things make this worth doing in these variables specifically:
//
//   * symplecticity is a statement about a canonical pair, so the state has to
//     be (x, p), not (x, v) - see hamiltonian.h;
//   * energy and axial angular momentum are then components of the state that
//     the equations never touch, so they are conserved bit-exactly;
//   * one right-hand side evaluation is *cheaper* than one christoffel()
//     evaluation (128 vs 216 flops in the general case), because contracting
//     g^{mu nu} with p twice is less work than contracting Gamma with v twice.
//
// What it costs: the stage equations are implicit and are solved here by fixed
// point iteration, so a step is several right-hand side evaluations rather
// than one, and the step size must stay constant for the backward-error
// argument to hold.  A constant step is a poor fit for a ray that starts a
// thousand M out and grazes the photon sphere, so the Sundman/Poincare time
// transformation below is provided as the symplectic way to vary it.

// Gauss-Legendre tableaus.  Coefficients are written out to more digits than
// double carries rather than computed with sqrt(), which is not constexpr.
// Getting one wrong is not a silent error: tests/test_symplectic.cpp measures
// the observed convergence order, which drops immediately if the tableau is
// not the collocation method it claims to be.
template <class FP, int S>
struct gauss_legendre;

// s=1: the implicit midpoint rule.  Order 2.
template <class FP>
struct gauss_legendre<FP, 1>
{
    static constexpr int order = 2;
    static constexpr FP A[1][1] = {{FP(0.5)}};
    static constexpr FP b[1] = {FP(1.0)};
    static constexpr FP c[1] = {FP(0.5)};
};

// s=2: order 4.  c = 1/2 -+ sqrt(3)/6.
template <class FP>
struct gauss_legendre<FP, 2>
{
    static constexpr FP R3 = FP(0.28867513459481288225L);   // sqrt(3)/6
    static constexpr int order = 4;
    static constexpr FP A[2][2] = {{FP(0.25),      FP(0.25) - R3},
                                   {FP(0.25) + R3, FP(0.25)}};
    static constexpr FP b[2] = {FP(0.5), FP(0.5)};
    static constexpr FP c[2] = {FP(0.5) - R3, FP(0.5) + R3};
};

// s=3: order 6.  c = 1/2, 1/2 -+ sqrt(15)/10.
template <class FP>
struct gauss_legendre<FP, 3>
{
    static constexpr FP R10 = FP(0.38729833462074168852L);  // sqrt(15)/10
    static constexpr FP R15 = FP(0.25819888974716112568L);  // sqrt(15)/15
    static constexpr FP R24 = FP(0.16137430609197570355L);  // sqrt(15)/24
    static constexpr FP R30 = FP(0.12909944487358056284L);  // sqrt(15)/30
    static constexpr int order = 6;
    static constexpr FP A[3][3] = {
        {FP(5.0/36.0),       FP(2.0/9.0) - R15, FP(5.0/36.0) - R30},
        {FP(5.0/36.0) + R24, FP(2.0/9.0),       FP(5.0/36.0) - R24},
        {FP(5.0/36.0) + R30, FP(2.0/9.0) + R15, FP(5.0/36.0)},
    };
    static constexpr FP b[3] = {FP(5.0/18.0), FP(4.0/9.0), FP(5.0/18.0)};
    static constexpr FP c[3] = {FP(0.5) - R10, FP(0.5), FP(0.5) + R10};
};


// The Sundman (Poincare) time transformation, i.e. the symplectic way to take
// short steps near the hole and long ones far from it.
//
// Rescaling the step directly would destroy the backward-error argument that
// is the entire reason to be here: the modified Hamiltonian a symplectic
// method conserves depends on h, so changing h changes which quantity is being
// conserved, and the constraint drifts again.  The fix is to keep the step
// constant and change the independent variable instead.  Integrating
//
//     K(x, p) = s(x) [ H(x, p) - H0 ],       s(x) = (r / r_ref)^kappa
//
// in a fictitious time tau, at fixed step in tau, is still a canonical
// Hamiltonian system - so the method is still symplectic - and on the level
// set K = 0 (which is where the ray starts, H0 being its initial H) its orbits
// are the orbits of H with dlambda = s(x) dtau.
//
//     dx/dtau =  s dH/dp
//     dp/dtau = -s dH/dx - (H - H0) grad s
//
// The second term vanishes on the exact solution but is what makes the system
// Hamiltonian, so it is kept.  grad s has only a radial component, ds/dr =
// kappa s / r.
//
// kappa = 0 recovers the plain, untransformed system with s = 1.  kappa = 2 is
// the natural choice for light: the deflection accumulated per unit affine
// parameter goes as M/r^2, so a step in tau covers roughly constant deflection.
// `factor` returns s(x) at the evaluation point, which the caller needs to
// accumulate the affine parameter.
template <class FP>
inline void sundman_rhs(kerr_black_hole<FP> const& hole, int const kappa, FP const r_ref, FP const H0,
                        FP const* const __restrict__ x, FP const* const __restrict__ p,
                        FP* const __restrict__ dx, FP* const __restrict__ dp, FP& factor)
{
    hamiltonian_rhs(hole, x, p, dx, dp);
    if (kappa == 0)
    {
        factor = FP(1);
        return;
    }
    FP hamiltonian = FP(0);
    for (int n = 0; n < D; ++n) hamiltonian += p[n] * dx[n];
    hamiltonian *= FP(0.5);

    FP const s = pown(x[1] / r_ref, kappa);
    for (int n = 0; n < D; ++n) { dx[n] *= s; dp[n] *= s; }
    dp[1] -= (hamiltonian - H0) * FP(kappa) * s / x[1];
    factor = s;
}


// One Gauss-Legendre step of size h, in place.
//
// The stage system is solved by fixed point iteration.  That is the right
// choice for this problem and not a shortcut: geodesics are non-stiff, so the
// iteration contracts at a rate ~ h * |df/dz| and needs no Jacobian, no linear
// solve, and no storage beyond the stages - which also keeps it usable inside
// a per-pixel GPU kernel, where a Newton solve on an 8s-dimensional system
// would not be.
//
// The iteration is run to `tol` on the stage increment rather than stopped at
// the truncation error, because an inexact stage solve perturbs the map by the
// residual and reintroduces exactly the secular drift the method exists to
// avoid (Hairer, Lubich & Wanner, VIII.6.2).  The residual is measured before
// the stages are updated, so the derivatives used for the final update are the
// ones evaluated at the converged stage values.
//
// Returns true if the iteration converged.  `evals` accumulates right-hand
// side evaluations, so callers can compare cost against the explicit methods
// on an equal footing, and `dlambda` reports the affine parameter actually
// advanced (h itself unless a time transformation is in use).
template <class FP, int S>
inline bool gauss_legendre_step(kerr_black_hole<FP> const& hole,
                                FP* const __restrict__ x, FP* const __restrict__ p,
                                FP const h, FP const tol, int const max_iter,
                                int const kappa, FP const r_ref, FP const H0,
                                FP& dlambda, int& evals)
{
    using GL = gauss_legendre<FP, S>;
    FP stage_x[S][D], stage_p[S][D], fx[S][D], fp[S][D], factor[S];

    // Constant predictor.  A stage predictor extrapolated from the previous
    // step would start closer, but it has to be carried across steps and the
    // measured iteration counts (see docs/integrator-comparison.md) did not
    // justify the extra state.
    for (int i = 0; i < S; ++i)
    {
        for (int n = 0; n < D; ++n) { stage_x[i][n] = x[n]; stage_p[i][n] = p[n]; }
    }

    bool converged = false;
    for (int iter = 0; iter < max_iter && !converged; ++iter)
    {
        for (int i = 0; i < S; ++i)
        {
            sundman_rhs(hole, kappa, r_ref, H0, stage_x[i], stage_p[i], fx[i], fp[i], factor[i]);
        }
        evals += S;

        // Residual of the stage equations at the point the derivatives were
        // just evaluated, scaled the same way the Dormand-Prince controller
        // scales its error estimate.
        FP residual = FP(0);
        for (int i = 0; i < S; ++i)
        {
            for (int n = 0; n < D; ++n)
            {
                FP sum_x = FP(0), sum_p = FP(0);
                for (int j = 0; j < S; ++j)
                {
                    sum_x += GL::A[i][j] * fx[j][n];
                    sum_p += GL::A[i][j] * fp[j][n];
                }
                FP const next_x = x[n] + h * sum_x;
                FP const next_p = p[n] + h * sum_p;
                residual = fmax(residual, fabs(next_x - stage_x[i][n]) / (FP(1) + fabs(next_x)));
                residual = fmax(residual, fabs(next_p - stage_p[i][n]) / (FP(1) + fabs(next_p)));
                stage_x[i][n] = next_x;
                stage_p[i][n] = next_p;
            }
        }
        converged = (residual <= tol);
    }

    for (int n = 0; n < D; ++n)
    {
        FP sum_x = FP(0), sum_p = FP(0);
        for (int i = 0; i < S; ++i)
        {
            sum_x += GL::b[i] * fx[i][n];
            sum_p += GL::b[i] * fp[i][n];
        }
        x[n] += h * sum_x;
        p[n] += h * sum_p;
    }

    // Affine parameter advanced, by the same quadrature as the state: exact
    // for the untransformed system, and the consistent estimate under a time
    // transformation.
    dlambda = h;
    if (kappa != 0)
    {
        FP sum = FP(0);
        for (int i = 0; i < S; ++i) sum += GL::b[i] * factor[i];
        dlambda = h * sum;
    }
    return converged;
}


// The same step with the plain (untransformed) Hamiltonian, for callers that
// do not want the time transformation.
template <class FP, int S>
inline bool gauss_legendre_step(kerr_black_hole<FP> const& hole,
                                FP* const __restrict__ x, FP* const __restrict__ p,
                                FP const h, FP const tol, int const max_iter, int& evals)
{
    FP dlambda;
    return gauss_legendre_step<FP, S>(hole, x, p, h, tol, max_iter, 0, FP(1), FP(0), dlambda, evals);
}


// Step the renderer's (x, v) state instead of (x, p).
//
// The conversion at each end is exact algebra (v^mu = g^{mu nu} p_nu and its
// inverse), so nothing about the integration changes; the integrator still
// runs on the canonical pair.  This exists so the symplectic methods can be
// dropped into code written against the second-order form - the collision
// tests, the disk crossing and the redshift in cuda_ray.h all want v.
//
// Converting per step costs one extra raise and one extra lower.  Code that
// takes many steps in a row should keep p and convert once at each end.
template <class FP, int S>
inline bool gauss_legendre_step_v(kerr_black_hole<FP> const& hole,
                                  FP* const __restrict__ x, FP* const __restrict__ v,
                                  FP const h, FP const tol, int const max_iter, int& evals)
{
    FP p[D];
    lower_velocity(hole, x, v, p);
    bool const ok = gauss_legendre_step<FP, S>(hole, x, p, h, tol, max_iter, evals);
    raise_momentum(hole, x, p, v);
    return ok;
}


// Fixed point iteration is a contraction only for h small enough against the
// local curvature scale, and a ray that is about to be captured can cross that
// threshold.  This retries a failed step as 2, 4, ... substeps.
//
// Halving the step breaks the constant-step condition the backward-error
// argument needs, so this is a safety net, not a step size controller: on any
// ray where it fires on more than a few steps, the step size is simply wrong.
template <class FP, int S>
inline bool gauss_legendre_step_safe(kerr_black_hole<FP> const& hole,
                                     FP* const __restrict__ x, FP* const __restrict__ p,
                                     FP const h, FP const tol, int const max_iter,
                                     int const kappa, FP const r_ref, FP const H0,
                                     FP& dlambda, int& evals)
{
    FP const x_in[D] = {x[0], x[1], x[2], x[3]};
    FP const p_in[D] = {p[0], p[1], p[2], p[3]};

    for (int split = 1; split <= 8; split *= 2)
    {
        FP advanced = FP(0);
        bool ok = true;
        for (int k = 0; k < split && ok; ++k)
        {
            FP piece;
            ok = gauss_legendre_step<FP, S>(hole, x, p, h / FP(split), tol, max_iter,
                                            kappa, r_ref, H0, piece, evals);
            advanced += piece;
        }
        if (ok)
        {
            dlambda = advanced;
            return true;
        }
        for (int n = 0; n < D; ++n) { x[n] = x_in[n]; p[n] = p_in[n]; }
    }
    dlambda = FP(0);
    return false;
}

#endif // SYMPLECTIC_H
