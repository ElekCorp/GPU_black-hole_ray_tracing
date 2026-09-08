// ===========================================================================
//  GENERATED FILE - DO NOT EDIT.
//
//  Metric  : kerr_newman
//  Chart   : Boyer-Lindquist (t, r, theta, phi); signature (+,-,-,-).
//
//  Regenerate with:  make -C codegen
//  Source of truth:  codegen/metrics/kerr_newman.py
//
//  Verification at generation time:
//    [PROVED                ] geodesic RHS: Christoffel formula == Euler-Lagrange
//    [PROVED                ] camera tetrad: e^T g e == eta
//    [PROVED                ] CSE (64 temporaries) preserves the expression
//    [PROVED                ] nonzero: expr[0]: 1/(a**2*cos(x2)**2 + x1**2)
//    [PROVED                ] nonzero: expr[0]: 1/(Q**2 + a**2 - rs*x1 + x1**2)
//    [PROVED                ] nonzero: expr[3]: 1/(sin(x2))
//    [NO_BOUND              ] round-off double acc[0]
//    [NO_BOUND              ] round-off double acc[1]
//    [NO_BOUND              ] round-off double acc[2]
//    [NO_BOUND              ] round-off double acc[3]
//    [NO_BOUND              ] round-off float acc[0]
//    [NO_BOUND              ] round-off float acc[1]
//    [NO_BOUND              ] round-off float acc[2]
//    [NO_BOUND              ] round-off float acc[3]
//
//  Code size:
//   metric            9 CSE temporaries
//   geodesic_rhs     64 CSE temporaries
//   observer          0 CSE temporaries
//   event_fields      1 CSE temporaries
//   event_guards      3 CSE temporaries
// ===========================================================================
#pragma once

#include "bh_metric_api.h"

namespace bh_kerr_newman
{

constexpr int N_PARAM = 7;
constexpr int N_EVENT = 4;
constexpr int N_GUARD = 4;

inline constexpr char const* PARAM_NAME[N_PARAM] = { "a", "Q", "rs", "eps_h", "R_sky", "r_in", "r_out" };
inline constexpr double PARAM_DEFAULT[N_PARAM] = { 0.0, 0.0, 0.05, 0.001, 1.01, 0.1, 0.5 };
inline constexpr char const* PARAM_DOC[N_PARAM] = { "spin, a = J/M", "electric charge", "Schwarzschild radius, rs = 2M", "stop this far outside the horizon, in units of Delta", "sky sphere radius; beyond it the ray has escaped", "accretion disk inner radius", "accretion disk outer radius" };

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
    FP const t0 = pown(x[1], 2);
    FP const t1 = pown(p[0], 2);
    FP const t2 = t0 + t1*pown(cos(x[2]), 2);
    FP const t3 = (FP(1) / t2);
    FP const t4 = pown(sin(x[2]), 2);
    FP const t5 = t1*t4;
    FP const t6 = pown(p[1], 2) - p[2]*x[1];
    FP const t7 = t0 + t1 + t6;
    FP const t8 = t3*t4;

    g[0][0] = t3*(-t5 + t7);
    g[0][1] = 0;
    g[0][2] = 0;
    g[0][3] = -p[0]*t6*t8;
    g[1][1] = -t2/t7;
    g[1][2] = 0;
    g[1][3] = 0;
    g[2][2] = -t2;
    g[2][3] = 0;
    g[3][3] = -t8*(pown(p[0], 4) + 2*t0*t1 - t5*t7 + pown(x[1], 4));

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
    FP const t0 = pown(x[1], 4);
    FP const t1 = p[2]*t0;
    FP const t2 = pown(x[1], 3);
    FP const t3 = pown(p[1], 2);
    FP const t4 = 2*t3;
    FP const t5 = t2*t4;
    FP const t6 = pown(p[0], 2);
    FP const t7 = cos(x[2]);
    FP const t8 = pown(t7, 2);
    FP const t9 = t6*t8;
    FP const t10 = pown(x[1], 2);
    FP const t11 = p[2]*t10;
    FP const t12 = t11*t9;
    FP const t13 = pown(p[0], 4);
    FP const t14 = p[2]*t13;
    FP const t15 = t14*t8;
    FP const t16 = -t11*t6;
    FP const t17 = t3*t6;
    FP const t18 = 2*x[1];
    FP const t19 = t17*t18;
    FP const t20 = t15 + t16 + t19;
    FP const t21 = t10 + t9;
    FP const t22 = pown(t21, 2);
    FP const t23 = t22*v[1];
    FP const t24 = t23*v[0];
    FP const t25 = 2*x[2];
    FP const t26 = sin(t25);
    FP const t27 = 2*t10;
    FP const t28 = pown(t27 + t6*cos(t25) + t6, 2);
    FP const t29 = p[2]*x[1];
    FP const t30 = -p[2]*t2 + t10*t3;
    FP const t31 = pown(p[1], 4) + pown(p[2], 2)*t10 + t17 - t29*t4 - t29*t6 + t30;
    FP const t32 = sin(x[2]);
    FP const t33 = pown(t32, 2);
    FP const t34 = t33*t6;
    FP const t35 = t10 + t6;
    FP const t36 = 2*v[3];
    FP const t37 = t7*v[2];
    FP const t38 = t36*t37;
    FP const t39 = 4*t2;
    FP const t40 = -t12;
    FP const t41 = t18*t3;
    FP const t42 = -t29 + t3;
    FP const t43 = t35 + t42;
    FP const t44 = (FP(1) / t43);
    FP const t45 = t44/pown(t21, 4);
    FP const t46 = pown(v[1], 2);
    FP const t47 = p[2] - 2*x[1];
    FP const t48 = t26*t6;
    FP const t49 = t43*v[2];
    FP const t50 = pown(v[0], 2);
    FP const t51 = p[2]*t9 - t11 + t41;
    FP const t52 = t36*v[0];
    FP const t53 = pown(v[3], 2);
    FP const t54 = t34*t43;
    FP const t55 = t0 + t13 + t27*t6;
    FP const t56 = -t54 + t55;
    FP const t57 = 4*x[1];
    FP const t58 = t44/pown(t21, 3);
    FP const t59 = t35*t42;
    FP const t60 = t21*(-2*t54 + t55) + t34*t56;
    FP const t61 = t13*t18;
    FP const t62 = pown(t8 - 1, 2);
    FP const t63 = -t29*t9 + t3*t9 + t30;

    acc[0] = -t45*(-pown(p[0], 3)*t31*pown(t32, 3)*t38*pown(-t34 + t35, 2) + p[0]*t22*t33*v[1]*v[3]*(-3*t1 + t20 + t3*t39 + t40 + t41*t9) - t24*(-t1 + t12 + t20 + t5) + (FP(1.0 / 4.0))*t26*t28*t31*t6*v[0]*v[2]);
    acc[1] = (FP(1.0 / 2.0))*t58*(2*t22*t43*v[2]*(t48*v[1] + t49*x[1]) - t22*t46*(t18*t43 + t21*t47) - pown(t43, 2)*(p[0]*t33*t51*t52 + t33*t53*(t18*t56 - t21*(t34*t47 + t39 + t57*t6)) - t50*t51));
    acc[2] = t58*((FP(1.0 / 2.0))*t22*t49*(t48*v[2] - t57*v[1]) - FP(1.0 / 8.0)*t28*t46*t48 + t32*t43*t7*(p[0]*t52*t59 - t42*t50*t6 + t53*t60));
    acc[3] = -t45*(2*p[0]*t22*t31*t37*v[0] - p[0]*t24*t32*t51 + t23*t32*v[3]*(-2*t1 + t13*t57*t8 - t14*t62 + t14 - t15 + t16 + t19 + t39*t9 + t40 + t5 + t61*t62 - t61 + 2*pown(x[1], 5)) + t38*(t34*t59*t63 + t60*(t0 + t13*pown(t7, 4) + t27*t9 + t63)))/t32;
}

/**
 * The camera's 4-velocity at x, up to normalisation - the time leg the
 * frame builder is seeded with.  Only its direction matters.
 */
template <class FP>
inline void observer(FP const* p, FP const* x, FP* u)
{
    u[0] = pown(p[0], 2) + pown(x[1], 2);
    u[1] = 0;
    u[2] = 0;
    u[3] = p[0];
}
/**
 * The scalar field of every event, evaluated at x.  An event fires on the
 * sign or sign change of its scalar; see EVENT_KIND.
 */
template <class FP>
inline void event_fields(FP const* p, FP const* x, FP* f)
{
    FP const t0 = x[2] - FP(1.0 / 2.0)*FP(3.14159265358979323846);

    f[0] = pown(p[0], 2) + pown(p[1], 2) - p[2]*x[1] - p[3] + pown(x[1], 2);
    f[1] = p[4] - x[1];
    f[2] = t0;
    f[3] = t0;
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
    FP const t1 = -p[5] - t0;
    FP const t2 = p[6] + t0;

    gd[0] = t1;
    gd[1] = t2;
    gd[2] = t1;
    gd[3] = t2;
}

} // namespace bh_kerr_newman
