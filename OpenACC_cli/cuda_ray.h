#ifndef CUDA_RAY_H
#define CUDA_RAY_H

#include <iostream>
#include <math.h>

#include "black_hole.h"

//#include "debugmalloc.h"


template <class FP>
inline void step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de);
template <class FP>
inline void step_size(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP& de);

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
    //RK38(hole, x, v, de);//RK38 vagy 6
    RK6(hole, x, v, de);
    step_size(hole, x, v, de);
}

template <class FP>
inline void step_size(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)
{
    FP ch[D];
    FP de0 = hole.de0;

    christoffel(hole, x, v, ch);

    FP err = hole.errormax;
    FP sum = 0.0;
    for (size_t i = 0; i < D; ++i)
    {
        sum += fabs(ch[i]);
    }

    de = sqrt(err / sum);

    if (de > de0)
    {
        de = de0;
    }
    else if (isnan(de))
    {
        de = de0 / 10;
    }
    else
    {

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
    
    FP y = x[2];
    
    FP x0 = pown(Q, 2) - rs*x[1];
    FP x1 = -x0;
    FP x2 = pown(a, 2);
    FP x3 = sin(y);
    FP x4 = pown(x3, 2);
    FP x5 = x2*x4;
    FP x6 = pown(x[1], 2);
    FP x7 = x0 + x2 + x6;
    FP x8 = x2*x4 - x7;
    FP x9 = x5*x7;
    FP x10 = 2*x2;
    FP x11 = x2*x2 + pown(x[1], 4) + x10*x6;
    FP x12 = x11 - x9;
    FP x13 = pown(-pown(x1, 2)*x5 + x12*x8, -1);
    FP x14 = cos(y);
    FP x15 = pown(x14, 2)*x2 + x6;
    FP x16 = pown(x15, -1);
    FP x17 = -rs + 2*x[1]*x1*x16;
    FP x18 = v[3]*x4;
    FP x19 = a*x17*x18;
    FP x20 = 2*x[1];
    FP x21 = x16*x5;
    FP x22 = -rs - 2*x[1]*x16*x7 + x20*x21 + x20;
    FP x23 = -x22;
    FP x24 = x21 + 1;
    FP x25 = x1*x24;
    FP x26 = x16*x7;
    FP x27 = x24 - x26;
    FP x28 = a*v[0];
    FP x29 = -v[3]*x25 + x27*x28;
    FP x30 = v[2]*x14;
    FP x31 = 2*a*x3*x30;
    FP x32 = rs - x20;
    FP x33 = 4*pown(x[1], 3) + 4*x[1]*x2;
    FP x34 = v[1]*x3;
    FP x35 = x11 - 2*x9;
    FP x36 = 2*v[2]*x14*(-v[3]*(x12*x21 + x35) + x25*x28) - x34*(v[3]*(-x12*x16*x20 + x32*x5 + x33) + x17*x28);
    FP x37 = a*x1;
    FP x38 = pown(v[2], 2);
    FP x39 = x18*x28;
    FP x40 = pown(x7, -1);
    FP x41 = pown(x15, -2);
    FP x42 = x1*x41;
    FP x43 = (1.0/2.0)*x16;
    FP x44 = pown(v[3], 2);
    FP x45 = x4*x44;
    FP x46 = pown(v[1], 2)*x40;
    FP x47 = x12*x41;
    FP x48 = pown(v[0], 2);
    FP x49 = x14*x3;
    FP x50 = x2*x49;
    FP x51 = x16*x49;
    FP x52 = 2*v[3];
    FP x53 = x14*pown(x3, 3);

    ch[0] = x13*(-x12*(v[1]*(v[0]*x23 + x19) + x29*x31) + x3*x36*x37);
    ch[1] = x26*(-rs*x16*x39 + x[1]*x38 - x[1]*x45*x47 + x10*x30*x34*x40 + x20*x39*x42 + x23*x43*x48 + x43*x45*(x32*x5 + x33) - 1.0/2.0*x46*(x15*x32*x40 + x20));
    ch[2] = -x16*(pown(a, 3)*v[0]*x42*x52*x53 + v[1]*v[2]*x20 + x1*x28*x51*x52 - x16*x27*x48*x50 - x2*x44*x47*x53 - x35*x44*x51 - x38*x50 + x46*x50);
    ch[3] = x13*(x37*(v[1]*(v[0]*x22 - x19) - x29*x31) + x36*x8/x3);
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
                //0 fekete, -1 hiba kezeléses piros, egyébként meg egy FP ami reprezental egy szint
                if (gomb_be(sugar_be, x))
                {
                    szin[i * MAGAS + j] = 0;//;
                    fut = false;
                }
                else if (gomb_ki(sugar_ki, x))
                {
                    szin[i * MAGAS + j] = 0;
                    fut = false;
                }
                else if (disk1(sugar_kicsi, sugar_nagy, x, x_le))
                {
                    szin[i * MAGAS + j] = x[1];
                    fut = false;
                }
                else if (disk2(sugar_kicsi, sugar_nagy, x, x_le))
                {
                    szin[i * MAGAS + j] = x[1];
                    fut = false;
                }
                else if (isnan(x[0]) || isnan(x[1]) || isnan(x[2]) || isnan(x[3]))
                {
                    //printf("%d\t%d\t%f\tnan\n", i, j, de);
                    szin[i * MAGAS + j] = -1;
                    fut = false;
                }
                else if (isinf(x[0]) || isinf(x[1]) || isinf(x[2]) || isinf(x[3]))
                {
                    //printf("%d\t%d\t%f\tinf\n", i, j, de);
                    szin[i * MAGAS + j] = -1;
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
