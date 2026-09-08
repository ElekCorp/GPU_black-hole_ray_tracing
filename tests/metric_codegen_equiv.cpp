// Bridge test: the GENERATED Kerr-Newman geodesic right-hand side must agree
// with the hand-generated christoffel() block the tracer currently ships.
//
// This is the load-bearing check on the whole codegen pipeline.  The shipped
// block was independently verified against symbolically derived Kerr-Newman
// Christoffel symbols, so agreement here transfers that confidence to the
// generated code - and any disagreement localises immediately to the generator.
//
// It also checks the generated metric() against the closed-form Boyer-Lindquist
// metric, and that the geodesic RHS is consistent with it via g(v, dv/de) = 0
// along the flow (the norm of a geodesic's tangent is conserved).

#include <cstdio>
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <random>

#include "black_hole.h"
#include "cuda_ray.h"
#include "generated/bh_kerr_newman_metric.h"

typedef double FP;
namespace KN = bh_kerr_newman;

int main()
{
    std::mt19937_64 rng(20260907);
    auto U = [&](double lo, double hi) {
        return std::uniform_real_distribution<double>(lo, hi)(rng);
    };

    double worst_rhs = 0.0, worst_metric = 0.0, worst_norm = 0.0;
    int samples = 0;

    for (int trial = 0; trial < 2000; ++trial)
    {
        FP const a = U(-0.6, 0.6), Q = U(0.0, 0.4), rs = U(0.01, 1.0);
        FP const r = U(2.0, 30.0), th = U(0.15, 3.0);

        FP const Delta = r * r - rs * r + a * a + Q * Q;
        if (Delta < 1e-2) continue;                 // outside the horizon only
        ++samples;

        FP x[4] = { U(-5, 5), r, th, U(-3, 3) };
        FP v[4] = { U(-2, 2), U(-2, 2), U(-2, 2), U(-2, 2) };
        FP p[KN::N_PARAM] = { a, Q, rs, 1e-3, 1e3, 3.0, 10.0 };

        // ---- generated RHS vs the shipped hand-generated block --------------
        FP gen[4];
        KN::geodesic_rhs<FP>(p, x, v, gen);

        FP const Omega[3] = { 0, 0, 0 };
        kerr_black_hole<FP> hole(64, 32, x, Omega, a, Q, rs,
                                 1e-10, 0.01, 1000, 0.5, 0.75, 1e3, 3.0, 10.0);
        FP ref[4];
        christoffel(hole, x, v, ref);

        for (int i = 0; i < 4; ++i)
        {
            FP const scale = std::max({std::fabs(gen[i]), std::fabs(ref[i]), FP(1e-8)});
            worst_rhs = std::max<double>(worst_rhs, std::fabs(gen[i] - ref[i]) / scale);
        }

        // ---- generated metric vs the closed form ---------------------------
        FP g[4][4];
        KN::metric<FP>(p, x, g);

        FP const S = r * r + a * a * std::cos(th) * std::cos(th);
        FP const s2 = std::sin(th) * std::sin(th);
        FP cf[4][4] = {};
        cf[0][0] = -a * a * s2 / S + Delta / S;
        cf[1][1] = -S / Delta;
        cf[2][2] = -S;
        cf[3][3] = (Delta * a * a * s2 * s2 - a * a * a * a * s2
                    - 2 * a * a * r * r * s2 - r * r * r * r * s2) / S;
        cf[0][3] = cf[3][0] = (a * a * a * s2 + a * r * r * s2 - Delta * a * s2) / S;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
            {
                FP const scale = std::max<FP>(std::fabs(cf[i][j]), 1e-8);
                worst_metric = std::max<double>(worst_metric, std::fabs(g[i][j] - cf[i][j]) / scale);
            }

        // ---- d/de of g(v,v) must vanish: 2 g(v, acc) + (d_c g) v v v^c = 0 ---
        // A cheap independent consistency check between metric() and
        // geodesic_rhs(): differentiate the norm numerically along the flow.
        FP dot_acc = 0.0;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                dot_acc += 2.0 * g[i][j] * v[i] * gen[j];

        FP dg_term = 0.0;
        FP const h = 1e-6;
        for (int c = 0; c < 4; ++c)
        {
            FP xp[4], xm[4], gp[4][4], gm[4][4];
            for (int k = 0; k < 4; ++k) { xp[k] = x[k]; xm[k] = x[k]; }
            xp[c] += h; xm[c] -= h;
            KN::metric<FP>(p, xp, gp);
            KN::metric<FP>(p, xm, gm);
            for (int i = 0; i < 4; ++i)
                for (int j = 0; j < 4; ++j)
                    dg_term += (gp[i][j] - gm[i][j]) / (2 * h) * v[i] * v[j] * v[c];
        }
        FP const mag = std::fabs(dot_acc) + std::fabs(dg_term) + 1e-8;
        worst_norm = std::max<double>(worst_norm, std::fabs(dot_acc + dg_term) / mag);
    }

    printf("samples                                          : %d\n", samples);
    printf("max rel. diff, generated vs shipped geodesic RHS : %.3e\n", worst_rhs);
    printf("max rel. diff, generated vs closed-form metric   : %.3e\n", worst_metric);
    printf("max rel. residual of d/de g(v,v) = 0             : %.3e\n", worst_norm);

    bool ok = worst_rhs < 1e-10 && worst_metric < 1e-12 && worst_norm < 1e-5;
    printf("\n%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
