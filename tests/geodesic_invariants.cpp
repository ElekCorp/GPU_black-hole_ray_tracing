// Regression test for the ray tracer's mathematics.
//
// Two things must hold for any correct null-geodesic integrator in Kerr-Newman:
//
//   1. The ray leaves the camera NULL:  g_{mu nu} v^mu v^nu == 0.
//      This is a pure statement about the camera tetrad in ijk_to_vec_zoom().
//
//   2. Along the ray, these stay constant to within integration error:
//        - the norm      g_{mu nu} v^mu v^nu            (still 0)
//        - the energy    E = -g_{t mu} v^mu             (d/dt is Killing)
//        - the ang. mom. L =  g_{phi mu} v^mu           (d/dphi is Killing)
//
// Run it against the OpenACC header, which is plain host C++.

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <algorithm>

#include "black_hole.h"
#include "cuda_ray.h"

typedef double FP;

static void metric(FP a, FP Q, FP rs, FP r, FP th, FP g[4][4])
{
    FP const S = r * r + a * a * cos(th) * cos(th);
    FP const D_ = r * r - rs * r + a * a + Q * Q;
    FP const A = rs * r - Q * Q;
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) g[i][j] = 0.0;
    g[0][0] = -(1.0 - A / S);
    g[0][3] = g[3][0] = -A * a * sin(th) * sin(th) / S;
    g[1][1] = S / D_;
    g[2][2] = S;
    g[3][3] = (r * r + a * a + A * a * a * sin(th) * sin(th) / S) * sin(th) * sin(th);
}

static FP dot(FP const g[4][4], FP const* u, FP const* w)
{
    FP s = 0.0;
    for (int i = 0; i < 4; ++i) for (int j = 0; j < 4; ++j) s += g[i][j] * u[i] * w[j];
    return s;
}

int main()
{
    FP const a = 0.30, Q = 0.10, rs = 1.00;
    FP const x0[4] = { 0.0, 12.0, 1.30, 0.0 };
    FP const Om[3] = { 0.37, 2.10, -0.55 };          // a general, non-degenerate orientation

    kerr_black_hole<FP> hole(/*SZELES*/ 64, /*MAGAS*/ 32, x0, Om,
                             a, Q, rs, /*errormax*/ 1e-10, /*de0*/ 0.5, /*max_steps*/ 200000,
                             /*kepernyo_high*/ 0.5, /*kepernyo_tav*/ 0.75,
                             /*sugar_ki*/ 1e3, /*kicsi*/ 3.0, /*nagy*/ 10.0);

    FP worst_null0 = 0.0, worst_null = 0.0, worst_E = 0.0, worst_L = 0.0;
    int rays = 0;

    for (uint64_t j = 4; j < 32; j += 7)
    {
        for (uint64_t i = 4; i < 64; i += 11)
        {
            FP x[4] = { hole.t_0, hole.r_0, hole.theta_0, hole.phi_0 };
            FP v[4];
            ijk_to_vec_zoom(i, j, hole, 64, 32, 0, 0, 64, v);

            FP g[4][4];
            metric(a, Q, rs, x[1], x[2], g);
            FP const n0 = dot(g, v, v);
            FP const E0 = -(g[0][0] * v[0] + g[0][3] * v[3]);
            FP const L0 = g[3][0] * v[0] + g[3][3] * v[3];
            worst_null0 = std::max(worst_null0, fabs(n0));

            FP de = hole.de0;
            for (int n = 0; n < 4000; ++n)
            {
                step(hole, x, v, de);
                if (!std::isfinite(x[1]) || x[1] < 1.2 || x[1] > 400.0) break;
            }
            if (!std::isfinite(x[1])) continue;

            metric(a, Q, rs, x[1], x[2], g);
            FP const n1 = dot(g, v, v);
            FP const E1 = -(g[0][0] * v[0] + g[0][3] * v[3]);
            FP const L1 = g[3][0] * v[0] + g[3][3] * v[3];

            worst_null = std::max(worst_null, fabs(n1));
            worst_E = std::max(worst_E, fabs(E1 - E0) / fabs(E0));
            worst_L = std::max(worst_L, fabs(L1 - L0) / std::max(fabs(L0), 1e-3));
            ++rays;
        }
    }

    printf("rays traced                             : %d\n", rays);
    printf("max |g(v,v)| at the camera              : %.3e\n", worst_null0);
    printf("max |g(v,v)| after integration          : %.3e\n", worst_null);
    printf("max relative drift in energy E          : %.3e\n", worst_E);
    printf("max relative drift in ang. momentum L   : %.3e\n", worst_L);

    bool const ok = worst_null0 < 1e-12 && worst_null < 1e-7 && worst_E < 1e-7 && worst_L < 1e-7;
    printf("\n%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
