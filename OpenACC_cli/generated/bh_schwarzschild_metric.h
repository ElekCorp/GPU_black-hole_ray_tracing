// ===========================================================================
//  GENERATED FILE - DO NOT EDIT.
//
//  Metric  : schwarzschild
//  Chart   : Schwarzschild (t, r, theta, phi); signature (+,-,-,-).
//
//  Regenerate with:  make -C codegen
//  Source of truth:  codegen/metrics/schwarzschild.py
//
//  Verification at generation time:
//    [PROVED                ] geodesic RHS: Christoffel formula == Euler-Lagrange
//    [PROVED                ] camera tetrad: e^T g e == eta
//    [PROVED                ] CSE (8 temporaries) preserves the expression
//    [PROVED                ] nonzero: expr[0]: 1/(x1)
//    [PROVED                ] nonzero: expr[0]: 1/(rs - x1)
//    [PROVED                ] nonzero: expr[3]: 1/(sin(x2))
//
//  Code size:
//   metric            2 CSE temporaries
//   geodesic_rhs      8 CSE temporaries
//   observer          0 CSE temporaries
//   event_fields      2 CSE temporaries
//   event_guards      3 CSE temporaries
// ===========================================================================
#pragma once

#include "bh_metric_api.h"

namespace bh_schwarzschild
{

constexpr int N_PARAM = 5;
constexpr int N_EVENT = 4;
constexpr int N_GUARD = 4;

inline constexpr char const* PARAM_NAME[N_PARAM] = { "rs", "eps_h", "R_sky", "r_in", "r_out" };
inline constexpr double PARAM_DEFAULT[N_PARAM] = { 0.05, 0.001, 1.01, 0.1, 0.5 };
inline constexpr char const* PARAM_DOC[N_PARAM] = { "Schwarzschild radius, rs = 2M", "stop this far outside the horizon", "sky sphere radius", "disk inner radius", "disk outer radius" };

inline constexpr char const* EVENT_NAME[N_EVENT] = { "HORIZON", "SKY", "DISK_TOP", "DISK_BOTTOM" };
inline constexpr bh_event_kind EVENT_KIND[N_EVENT] = { EVENT_REGION_NEG, EVENT_REGION_NEG, EVENT_CROSS_NEG_TO_POS, EVENT_CROSS_POS_TO_NEG };
inline constexpr int EVENT_N_GUARD[N_EVENT] = { 0, 0, 2, 2 };
inline constexpr int EVENT_GUARD_OFFSET[N_EVENT] = { 0, 0, 0, 2 };

//! Signature convention: +1 means (+,-,-,-), -1 means (-,+,+,+).
constexpr int SIGNATURE = 1;
constexpr int EV_HORIZON = 0;
constexpr int EV_SKY = 1;
constexpr int EV_DISK_TOP = 2;
constexpr int EV_DISK_BOTTOM = 3;

/**
 * The metric g_{mu nu} at x.  Symmetric; only the upper triangle is
 * computed and the rest is mirrored.
 */
template <class FP>
inline void metric(FP const* p, FP const* x, FP g[4][4])
{
    FP const t0 = p[0]/x[1] - 1;
    FP const t1 = pown(x[1], 2);

    g[0][0] = -t0;
    g[0][1] = 0;
    g[0][2] = 0;
    g[0][3] = 0;
    g[1][1] = (FP(1) / t0);
    g[1][2] = 0;
    g[1][3] = 0;
    g[2][2] = -t1;
    g[2][3] = 0;
    g[3][3] = -t1*pown(sin(x[2]), 2);

    g[1][0] = g[0][1];
    g[2][0] = g[0][2];
    g[3][0] = g[0][3];
    g[2][1] = g[1][2];
    g[3][1] = g[1][3];
    g[3][2] = g[2][3];
}

/**
 * Geodesic acceleration acc^mu = -Gamma^mu_{alpha beta} v^alpha v^beta,
 * i.e. the right-hand side of dv^mu / de = acc^mu.

 * The Christoffel symbols are never formed: the full double contraction is
 * simplified and common-subexpression eliminated as one block.
 */
template <class FP>
inline void geodesic_rhs(FP const* p, FP const* x, FP const* v, FP* acc)
{
    FP const t0 = (FP(1) / x[1]);
    FP const t1 = p[0] - x[1];
    FP const t2 = (FP(1) / t1);
    FP const t3 = pown(x[1], 3);
    FP const t4 = pown(t1, 2);
    FP const t5 = pown(v[3], 2);
    FP const t6 = sin(x[2]);
    FP const t7 = 2*t0;

    acc[0] = p[0]*t0*t2*v[0]*v[1];
    acc[1] = -FP(1.0 / 2.0)*t2*(-p[0]*t4*pown(v[0], 2) + p[0]*pown(v[1], 2)*pown(x[1], 2) + 2*t3*t4*(t5*pown(t6, 2) + pown(v[2], 2)))/t3;
    acc[2] = (FP(1.0 / 2.0))*t5*sin(2*x[2]) - t7*v[1]*v[2];
    acc[3] = -t7*v[3]*(v[1] + v[2]*x[1]*cos(x[2])/t6);
}

/**
 * The camera's 4-velocity at x, up to normalisation - the time leg the
 * frame builder is seeded with.  Only its direction matters.
 */
template <class FP>
inline void observer(FP const* p, FP const* x, FP* u)
{
    u[0] = 1;
    u[1] = 0;
    u[2] = 0;
    u[3] = 0;
}
/**
 * The scalar field of every event, evaluated at x.  An event fires on the
 * sign or sign change of its scalar; see EVENT_KIND.
 */
template <class FP>
inline void event_fields(FP const* p, FP const* x, FP* f)
{
    FP const t0 = -x[1];
    FP const t1 = x[2] - FP(1.0 / 2.0)*FP(3.14159265358979323846);

    f[0] = -p[0] - p[1] - t0;
    f[1] = p[2] + t0;
    f[2] = t1;
    f[3] = t1;
}

/**
 * Guard scalars, concatenated over all events in order.  Event i owns the
 * slice [EVENT_GUARD_OFFSET[i], EVENT_GUARD_OFFSET[i] + EVENT_N_GUARD[i]);
 * it can only fire while every one of them is positive.
 */
template <class FP>
inline void event_guards(FP const* p, FP const* x, FP* gd)
{
    FP const t0 = -x[1];
    FP const t1 = -p[3] - t0;
    FP const t2 = p[4] + t0;

    gd[0] = t1;
    gd[1] = t2;
    gd[2] = t1;
    gd[3] = t2;
}

} // namespace bh_schwarzschild
