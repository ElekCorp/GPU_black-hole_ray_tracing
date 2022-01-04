#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "black_hole.cuh"

#include <iostream>

template <class FP>
__device__ void step(kerr_black_hole<FP>& hole, FP* x, FP* v, FP& de);
template <class FP>
__device__ void step_size(kerr_black_hole<FP>& hole, FP* x, FP* v, FP& de);

template <class FP>
__device__ FP ijk_to_vec_mink_zoom(int i, int j, int k, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg, kerr_black_hole<FP>& hole);
template <class FP>
__device__ FP ijk_to_vec_zoom(int i, int j, int k, kerr_black_hole<FP>& hole, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg);



template <class FP>
__device__ void RK38(kerr_black_hole<FP>& hole, FP* x, FP* v, FP de);//4-ed foku legpontosabb
template <class FP>
__device__ void RK6(kerr_black_hole<FP>& hole, FP* x, FP* v, FP de);

template <class FP>
__device__ void christoffel(kerr_black_hole<FP>& hole, FP* x, FP* v, FP* ch);

template <class FP>
__global__ void ray_step(int8_t* szin, int SZELES, int MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg);
template <class FP>
__global__ void ray_step_T(FP* szin, int SZELES, int MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg);


template <class FP>
__device__ bool gomb_be(FP sugar, FP* x);
template <class FP>
__device__ bool gomb_ki(FP sugar, FP* x);
template <class FP>
__device__ bool disk(FP sugar_kicsi, FP sugar_nagy, FP* x1, FP* x2);

template <class FP>
__device__ bool disk1(FP sugar_kicsi, FP sugar_nagy, FP* x1, FP* x2);
template <class FP>
__device__ bool disk2(FP sugar_kicsi, FP sugar_nagy, FP* x1, FP* x2);

template <class FP>
__host__ __device__ int ijk_to_n(int i, int j, int k, kerr_black_hole<FP>& hole);

template <class FP>
__device__ FP pown(FP x, int n);





template <class FP>
inline __device__ void step(kerr_black_hole<FP>& hole, FP* x, FP* v, FP& de)//adaptiv step size
{
    //RK38(hole, x, v, de);//RK38 vagy 6
    RK6(hole, x, v, de);
    step_size(hole, x, v, de);
}

template <class FP>
inline __device__ void step_size(kerr_black_hole<FP>& hole, FP* x, FP* v, FP& de)
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
inline __device__ void RK38(kerr_black_hole<FP>& hole, FP* x, FP* v, FP de)
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
inline __device__ void RK6(kerr_black_hole<FP>& hole, FP* x, FP* v, FP de)
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
inline __device__ void christoffel(kerr_black_hole<FP>& hole, FP* x, FP* v, FP* ch)
{
    FP a = hole.a;
    FP Q = hole.Q;
    FP rs = hole.rs;

    //original
    /*
    ch[0] = 1.0 / (((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (a * v[1] * (sin(x[2])) * (sin(x[2])) * (a * rs * (2 * a * v[3] - v[0]) * (x[1] * x[1]) - rs * v[0] * a * a * a - rs * v[3] * a * a * a * a + 3 * rs * v[3] * (x[1] * x[1] * x[1] * x[1]) - 4 * v[3] * x[1] * Q * Q * a * a - 4 * v[3] * Q * Q * x[1] * x[1] * x[1]) + rs * v[0] * v[1] * (a * a * a * a) - rs * v[0] * v[1] * x[1] * x[1] * x[1] * x[1] + 2 * v[0] * v[1] * x[1] * (Q * Q) * (a * a) + 2 * v[0] * v[1] * (Q * Q) * (x[1] * x[1] * x[1]) + v[0] * v[2] * sin(2 * x[2]) * (a * a) * (rs * x[1] * (2 * (Q * Q) + a * a) + rs * (x[1] * x[1] * x[1]) - Q * Q * (Q * Q + a * a) - x[1] * x[1] * (Q * Q + rs * rs)) + v[1] * v[3] * (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * (a * a * a) * (rs * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - 2 * v[2] * v[3] * (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * cos(x[2]) * a * a * a * (rs * x[1] * (2 * (Q * Q) + a * a) + rs * (x[1] * x[1] * x[1]) - Q * Q * (Q * Q + a * a) - x[1] * x[1] * (Q * Q + rs * rs)));

    ch[1] = (1.0 / 2.0) * 1.0 / (((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-rs * v[0] * v[0] * x[1] * x[1] * (6 * (Q * Q) * (a * a) + 5 * (Q * Q * Q * Q) + a * a * a * a) - 4 * rs * v[2] * v[2] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] - rs * x[1] * x[1] * x[1] * x[1] * ((a * a) * (v[1] * v[1]) + (v[0] * v[0]) * (6 * (Q * Q) + 2 * (a * a) + rs * rs)) + rs * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) * (-v[0] * v[0] + v[1] * v[1] - 4 * v[2] * v[2] * (Q * Q + a * a)) - 1.0 / 2.0 * v[1] * v[2] * sin(2 * x[2]) * a * a * (cos(2 * x[2]) * (a * a) + a * a + 2 * (x[1] * x[1])) * (cos(2 * x[2]) * (a * a) + a * a + 2 * (x[1] * x[1])) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * x[1] * (Q * Q) * (v[0] * v[0]) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) + (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * (a * a * a * a) * (v[3] * v[3]) * (-rs * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) - rs * x[1] * x[1] * (6 * (Q * Q) + 6 * (a * a) + rs * rs) - 5 * rs * x[1] * x[1] * x[1] * x[1] + 2 * x[1] * (2 * (Q * Q) * (a * a) + (rs * rs) * (Q * Q + a * a) + Q * Q * Q * Q + a * a * a * a) + 4 * (x[1] * x[1] * x[1]) * (Q * Q + a * a + rs * rs) + 2 * (x[1] * x[1] * x[1] * x[1] * x[1])) + (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * (a * a) * (v[3] * v[3]) * (rs * (a * a) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) + rs * (x[1] * x[1]) * (4 * (Q * Q) * (a * a) + (a * a) * (rs * rs) - 5 * Q * Q * Q * Q + 9 * (a * a * a * a)) + rs * (x[1] * x[1] * x[1] * x[1]) * (2 * (Q * Q) + 15 * (a * a) - rs * rs) + 7 * rs * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) - 2 * x[1] * (3 * (Q * Q) * (a * a * a * a) + (a * a) * (rs * rs) * (Q * Q + a * a) - Q * Q * Q * Q * Q * Q + 2 * (a * a * a * a * a * a)) + 4 * (x[1] * x[1] * x[1]) * (-3 * Q * Q * a * a + (rs * rs) * (Q * Q - a * a) - 3 * a * a * a * a) - 2 * x[1] * x[1] * x[1] * x[1] * x[1] * (3 * (Q * Q) + 6 * (a * a) + rs * rs) - 4 * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) + (sin(x[2])) * (sin(x[2])) * (2 * a * rs * v[3] * (x[1] * x[1]) * (v[0] * (6 * (Q * Q) * (a * a) + 5 * (Q * Q * Q * Q) + a * a * a * a) - 2 * v[3] * a * a * a * (Q * Q + a * a)) + a * rs * (x[1] * x[1] * x[1] * x[1]) * (a * (v[1] * v[1]) - 4 * a * v[3] * v[3] * (2 * (Q * Q) + 3 * (a * a)) + 2 * v[0] * v[3] * (6 * (Q * Q) + 2 * (a * a) + rs * rs)) - 2 * a * v[3] * x[1] * (2 * v[0] * (Q * Q) - v[3] * a * a * a) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) + 2 * a * v[3] * (x[1] * x[1] * x[1]) * (a * v[3] * (6 * (Q * Q) * (a * a) + (a * a) * (rs * rs) + 2 * (Q * Q * Q * Q) + 4 * (a * a * a * a)) - 2 * v[0] * (2 * (Q * Q) * (a * a) + (rs * rs) * (2 * (Q * Q) + a * a) + 2 * (Q * Q * Q * Q))) + 2 * rs * v[3] * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) * (a * v[0] - 2 * v[3] * (Q * Q + 3 * (a * a))) - 4 * rs * v[3] * v[3] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] + (rs - 2 * x[1]) * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (a * a * a * a * a * a) * (v[1] * v[1]) + 2 * (cos(x[2])) * (cos(x[2])) * (a * a * a) * (-rs * v[0] * v[3] * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) - rs * v[0] * v[3] * x[1] * x[1] * x[1] * x[1] + rs * (x[1] * x[1]) * (a * (v[1] * v[1]) - v[0] * v[3] * (2 * (Q * Q) + 2 * (a * a) + rs * rs)) + 2 * v[0] * v[3] * x[1] * (rs * rs) * (Q * Q + a * a) + 2 * (x[1] * x[1] * x[1]) * (-a * v[1] * v[1] + v[0] * v[3] * (rs * rs))) + 2 * (v[3] * v[3]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) * (2 * (Q * Q) + 4 * (a * a) + rs * rs) + 2 * (v[3] * v[3]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) + 2 * (x[1] * x[1] * x[1] * x[1] * x[1]) * (-2 * a * v[0] * v[3] * (Q * Q + rs * rs) - a * a * v[1] * v[1] + (v[3] * v[3]) * (6 * (Q * Q) * (a * a) + 2 * (a * a) * (rs * rs) + Q * Q * Q * Q + 6 * (a * a * a * a)))) + (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (a * a * a * a) * (-rs * a * a * v[1] * v[1] - 4 * rs * v[2] * v[2] * x[1] * x[1] * x[1] * x[1] + rs * (x[1] * x[1]) * (v[1] * v[1] - 4 * v[2] * v[2] * (Q * Q + a * a)) - 2 * x[1] * ((Q * Q) * (v[1] * v[1]) - v[2] * v[2] * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a)) + 2 * (v[2] * v[2]) * (x[1] * x[1] * x[1]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs) + 2 * (v[2] * v[2]) * (x[1] * x[1] * x[1] * x[1] * x[1])) + (cos(x[2])) * (cos(x[2])) * (a * a) * (rs * (v[0] * v[0]) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) - 8 * rs * v[2] * v[2] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] + rs * (x[1] * x[1]) * (-2 * a * a * v[1] * v[1] + (v[0] * v[0]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs)) + rs * (x[1] * x[1] * x[1] * x[1]) * (v[0] * v[0] + 2 * (v[1] * v[1]) - 8 * v[2] * v[2] * (Q * Q + a * a)) - 2 * x[1] * rs * rs * v[0] * v[0] * (Q * Q + a * a) + 4 * (v[2] * v[2]) * (x[1] * x[1] * x[1] * x[1] * x[1]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs) + 4 * (v[2] * v[2]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) + 2 * (x[1] * x[1] * x[1]) * (-2 * Q * Q * v[1] * v[1] - rs * rs * v[0] * v[0] + 2 * (v[2] * v[2]) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a))) + 2 * (v[0] * v[0]) * (x[1] * x[1] * x[1]) * (2 * (Q * Q) * (a * a) + (rs * rs) * (2 * (Q * Q) + a * a) + 2 * (Q * Q * Q * Q)) + 2 * (v[2] * v[2]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs) + 2 * (v[2] * v[2]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) + 2 * (x[1] * x[1] * x[1] * x[1] * x[1]) * (-Q * Q * v[1] * v[1] + (v[0] * v[0]) * (Q * Q + rs * rs) + (v[2] * v[2]) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a)));

    ch[2] = 1.0 / (((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1])) / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (2 * rs * v[1] * v[2] * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) + 2 * v[1] * v[2] * x[1] * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (a * a * a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 4 * v[1] * v[2] * (cos(x[2])) * (cos(x[2])) * (a * a) * (x[1] * x[1] * x[1]) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) - 2 * v[1] * v[2] * x[1] * x[1] * x[1] * x[1] * x[1] * (Q * Q + a * a) - 2 * v[1] * v[2] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] - sin(x[2]) * cos(x[2]) * (a * rs * (x[1] * x[1] * x[1]) * (-a * v[0] * v[0] - a * v[3] * v[3] * (4 * (Q * Q) + 3 * (a * a)) + 4 * v[0] * v[3] * (Q * Q + a * a)) + a * rs * (x[1] * x[1] * x[1] * x[1] * x[1]) * (a * (v[2] * v[2]) - a * v[3] * v[3] + 2 * v[0] * v[3]) + a * (x[1] * x[1]) * (a * (v[0] * v[0]) * (Q * Q + rs * rs) + a * (v[3] * v[3]) * (3 * (Q * Q) * (a * a) + (a * a) * (rs * rs) + 2 * (Q * Q * Q * Q)) - 2 * v[0] * v[3] * (2 * (Q * Q) * (a * a) + (a * a) * (rs * rs) + Q * Q * Q * Q)) + a * (x[1] * x[1] * x[1] * x[1]) * (a * (v[1] * v[1]) - a * v[2] * v[2] * (Q * Q + a * a) + a * (v[3] * v[3]) * (Q * Q - a * a + 2 * (rs * rs)) - 2 * v[0] * v[3] * (Q * Q + rs * rs)) - rs * x[1] * a * a * (2 * (Q * Q) + a * a) * (-2 * a * v[0] * v[3] + (a * a) * (v[3] * v[3]) + v[0] * v[0]) + rs * (v[3] * v[3]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) - 2 * v[0] * v[3] * Q * Q * a * a * a * (Q * Q + a * a) + (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (a * a * a * a) * (rs * x[1] * ((a * a) * (v[2] * v[2]) + 2 * (v[3] * v[3]) * (Q * Q + a * a)) + 2 * rs * (v[3] * v[3]) * (x[1] * x[1] * x[1]) + (a * a) * (v[1] * v[1]) - a * a * v[2] * v[2] * (Q * Q + a * a) - v[3] * v[3] * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) - v[3] * v[3] * x[1] * x[1] * x[1] * x[1] - x[1] * x[1] * ((a * a) * (v[2] * v[2]) + (v[3] * v[3]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs))) + 2 * (cos(x[2])) * (cos(x[2])) * (a * a) * (x[1] * x[1]) * (rs * x[1] * ((a * a) * (v[2] * v[2]) + 2 * (v[3] * v[3]) * (Q * Q + a * a)) + 2 * rs * (v[3] * v[3]) * (x[1] * x[1] * x[1]) + (a * a) * (v[1] * v[1]) - a * a * v[2] * v[2] * (Q * Q + a * a) - v[3] * v[3] * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) - v[3] * v[3] * x[1] * x[1] * x[1] * x[1] - x[1] * x[1] * ((a * a) * (v[2] * v[2]) + (v[3] * v[3]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs))) + (Q * Q) * (a * a) * (v[0] * v[0]) * (Q * Q + a * a) + (Q * Q) * (a * a * a * a) * (v[3] * v[3]) * (Q * Q + a * a) - v[3] * v[3] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] - x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * ((a * a) * (v[2] * v[2]) + (v[3] * v[3]) * (Q * Q + 2 * (a * a)))));

    ch[3] = 1.0 / (((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1])) / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-v[1] * sin(x[2]) * (a * rs * (-a * v[3] + v[0]) * (x[1] * x[1]) + 2 * a * x[1] * (a * v[3] - v[0]) * (Q * Q) - 2 * rs * v[3] * x[1] * x[1] * x[1] * x[1] + v[3] * (-rs + 2 * x[1]) * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (a * a * a * a) + 2 * v[3] * (Q * Q) * (x[1] * x[1] * x[1]) + 2 * v[3] * (x[1] * x[1] * x[1] * x[1] * x[1]) + (cos(x[2])) * (cos(x[2])) * (a * a) * (-a * rs * v[0] + rs * v[3] * (a * a) - rs * v[3] * x[1] * x[1] + 4 * v[3] * (x[1] * x[1] * x[1]))) + 2 * v[2] * v[3] * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (a * a * a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * v[2] * v[3] * (cos(x[2])) * (cos(x[2])) * (cos(x[2])) * (a * a) * (rs * x[1] * (2 * (Q * Q) + a * a) + 3 * rs * (x[1] * x[1] * x[1]) - Q * Q * (Q * Q + a * a) - x[1] * x[1] * (3 * (Q * Q) + 2 * (a * a) + rs * rs) - 2 * x[1] * x[1] * x[1] * x[1]) - 2 * v[2] * cos(x[2]) * (-a * rs * x[1] * (-a * v[3] + v[0]) * (2 * (Q * Q) + a * a) + a * rs * (a * v[3] - v[0]) * (x[1] * x[1] * x[1]) + a * v[0] * (Q * Q) * (Q * Q + a * a) + a * (-a * v[3] + v[0]) * (x[1] * x[1]) * (Q * Q + rs * rs) - rs * v[3] * x[1] * x[1] * x[1] * x[1] * x[1] - v[3] * Q * Q * a * a * (Q * Q + a * a) + v[3] * (x[1] * x[1] * x[1] * x[1]) * (Q * Q + a * a) + v[3] * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1]))) / sin(x[2]);
    */

    //ez a gyorsabb
    //-use_fast_math flag esetén nan-ba vagy inf-be megy néhány pont
/*
    ch[0] = 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0f / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (a * v[1] * v[3] * (sinf(x[2])) * (sinf(x[2])) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) * (x[1] * x[1]) - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * a * a + rs * (a * a) * (x[1] * x[1]) + 3 * rs * (x[1] * x[1] * x[1] * x[1]) - 2 * x[1] * (cosf(x[2])) * (cosf(x[2])) * Q * Q * a * a - 2 * x[1] * Q * Q * a * a - 4 * Q * Q * x[1] * x[1] * x[1]) + a * v[2] * v[3] * (-rs * x[1] + Q * Q) * (-2 * sinf(x[2]) * cosf(x[2]) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) + sinf(2 * x[2]) * (a * a + x[1] * x[1]) * (rs * x[1] * (sinf(x[2])) * (sinf(x[2])) * (a * a) - (sinf(x[2])) * (sinf(x[2])) * Q * Q * a * a - (sinf(x[2])) * (sinf(x[2])) * a * a * x[1] * x[1] - (sinf(x[2])) * (sinf(x[2])) * a * a * a * a + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) + v[0] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * (sinf(x[2])) * (sinf(x[2])) * a * a * x[1] * x[1] - rs * (sinf(x[2])) * (sinf(x[2])) * a * a * a * a + rs * (a * a * a * a) - rs * x[1] * x[1] * x[1] * x[1] + 2 * x[1] * (Q * Q) * (a * a) + 2 * (Q * Q) * (x[1] * x[1] * x[1])) - v[0] * v[2] * sinf(2 * x[2]) * a * a * (-rs * x[1] + Q * Q) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]));

    ch[1] = (1.0f / 2.0f) * 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (2 * v[1] * v[2] * sinf(2 * x[2]) * (a * a) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + 2 * x[1] * (v[2] * v[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-2 * a * v[0] * v[3] * (sinf(x[2])) * (sinf(x[2])) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - (sinf(x[2])) * (sinf(x[2])) * v[3] * v[3] * (2 * x[1] * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) - ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (4 * x[1] * (a * a) + (rs - 2 * x[1]) * (sinf(x[2])) * (sinf(x[2])) * (a * a) + 4 * (x[1] * x[1] * x[1]))) + (v[0] * v[0]) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q))) - v[1] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * x[1] * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (rs - 2 * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])));

    ch[2] = (1.0f / 2.0f) * 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-4 * v[1] * v[2] * x[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) - sinf(2 * x[2]) * a * a * v[1] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) + sinf(2 * x[2]) * (a * a) * (v[2] * v[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (2 * a * v[0] * v[3] * sinf(2 * x[2]) * (-rs * x[1] + Q * Q) * (a * a + x[1] * x[1]) + 2 * sinf(x[2]) * cosf(x[2]) * (v[3] * v[3]) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) - sinf(2 * x[2]) * a * a * v[0] * v[0] * (-rs * x[1] + Q * Q)) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]));

    ch[3] = 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-2 * a * v[0] * v[2] * (-rs * x[1] + Q * Q) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + v[1] * tan(x[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (a * v[0] * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - v[3] * (-rs * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * a * a * a * a - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * x[1] * x[1] - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * a * a - rs * a * a * x[1] * x[1] + rs * (a * a * a * a) - 2 * rs * x[1] * x[1] * x[1] * x[1] + 2 * x[1] * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (a * a * a * a) + 4 * x[1] * (cosf(x[2])) * (cosf(x[2])) * (a * a * a * a) + 2 * x[1] * (Q * Q) * (a * a) - 2 * x[1] * a * a * a * a + 4 * (cosf(x[2])) * (cosf(x[2])) * (a * a) * (x[1] * x[1] * x[1]) + 2 * (Q * Q) * (x[1] * x[1] * x[1]) + 2 * (x[1] * x[1] * x[1] * x[1] * x[1]))) - v[2] * v[3] * (sinf(2 * x[2]) * tan(x[2]) * (a * a) * (-rs * x[1] + Q * Q) * (-rs * x[1] + Q * Q) * (a * a + x[1] * x[1]) + 2 * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) * (-rs * x[1] + (cosf(x[2])) * (cosf(x[2])) * (a * a) + Q * Q + x[1] * x[1]))) / tan(x[2]);
    */

    /*ch[1] = (1.0 / 2.0) * 1.0 / (((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (2 * v[1] * v[2] * sin(2 * x[2]) * (a * a) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + 2 * x[1] * (v[2] * v[2]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-2 * a * v[0] * v[3] * (sin(x[2])) * (sin(x[2])) * (rs * (cos(x[2])) * (cos(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - (sin(x[2])) * (sin(x[2])) * v[3] * v[3] * (2 * x[1] * ((sin(x[2])) * (sin(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) - ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * (4 * x[1] * (a * a) + (rs - 2 * x[1]) * (sin(x[2])) * (sin(x[2])) * (a * a) + 4 * (x[1] * x[1] * x[1]))) + (v[0] * v[0]) * (rs * (cos(x[2])) * (cos(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q))) - v[1] * v[1] * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * (2 * x[1] * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (rs - 2 * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1])));

    ch[2] = (1.0 / 2.0) * 1.0 / (((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-4 * v[1] * v[2] * x[1] * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) - sin(2 * x[2]) * a * a * v[1] * v[1] * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) + sin(2 * x[2]) * (a * a) * (v[2] * v[2]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (2 * a * v[0] * v[3] * sin(2 * x[2]) * (-rs * x[1] + Q * Q) * (a * a + x[1] * x[1]) + 2 * sin(x[2]) * cos(x[2]) * (v[3] * v[3]) * ((sin(x[2])) * (sin(x[2])) * (a * a) * ((sin(x[2])) * (sin(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sin(x[2])) * (sin(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) - sin(2 * x[2]) * a * a * v[0] * v[0] * (-rs * x[1] + Q * Q)) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]));

    ch[3] = 1.0 / (((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-2 * a * v[0] * v[2] * (-rs * x[1] + Q * Q) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + v[1] * tan(x[2]) * ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * (a * v[0] * (rs * (cos(x[2])) * (cos(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - v[3] * (-rs * (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * a * a * a * a - rs * (cos(x[2])) * (cos(x[2])) * a * a * x[1] * x[1] - rs * (cos(x[2])) * (cos(x[2])) * a * a * a * a - rs * a * a * x[1] * x[1] + rs * (a * a * a * a) - 2 * rs * x[1] * x[1] * x[1] * x[1] + 2 * x[1] * (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * (sin(x[2])) * (a * a * a * a) + 4 * x[1] * (cos(x[2])) * (cos(x[2])) * (a * a * a * a) + 2 * x[1] * (Q * Q) * (a * a) - 2 * x[1] * a * a * a * a + 4 * (cos(x[2])) * (cos(x[2])) * (a * a) * (x[1] * x[1] * x[1]) + 2 * (Q * Q) * (x[1] * x[1] * x[1]) + 2 * (x[1] * x[1] * x[1] * x[1] * x[1]))) - v[2] * v[3] * (sin(2 * x[2]) * tan(x[2]) * (a * a) * (-rs * x[1] + Q * Q) * (-rs * x[1] + Q * Q) * (a * a + x[1] * x[1]) + 2 * ((sin(x[2])) * (sin(x[2])) * (a * a) * ((sin(x[2])) * (sin(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cos(x[2])) * (cos(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sin(x[2])) * (sin(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) * (-rs * x[1] + (cos(x[2])) * (cos(x[2])) * (a * a) + Q * Q + x[1] * x[1]))) / tan(x[2]);
    */



    //double esetén ez a leg gyorsabb
    
    ch[0] = pown(pown(cos(x[2]), 2) * (a * a) + x[1] * x[1], -2) * pown(-rs * x[1] + Q * Q + a * a + x[1] * x[1], -1) * (a * v[1] * pown(sin(x[2]), 2) * (a * rs * (2 * a * v[3] - v[0]) * (x[1] * x[1]) - rs * v[0] * a * a * a - rs * v[3] * a * a * a * a + 3 * rs * v[3] * (x[1] * x[1] * x[1] * x[1]) - 4 * v[3] * x[1] * Q * Q * a * a - 4 * v[3] * Q * Q * x[1] * x[1] * x[1]) + rs * v[0] * v[1] * (a * a * a * a) - rs * v[0] * v[1] * x[1] * x[1] * x[1] * x[1] + 2 * v[0] * v[1] * x[1] * (Q * Q) * (a * a) + 2 * v[0] * v[1] * (Q * Q) * (x[1] * x[1] * x[1]) + v[0] * v[2] * sin(2 * x[2]) * (a * a) * (rs * x[1] * (2 * (Q * Q) + a * a) + rs * (x[1] * x[1] * x[1]) - Q * Q * (Q * Q + a * a) - x[1] * x[1] * (Q * Q + rs * rs)) + v[1] * v[3] * pown(sin(x[2]), 4) * (a * a * a) * (rs * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - 2 * v[2] * v[3] * pown(sin(x[2]), 3) * cos(x[2]) * a * a * a * (rs * x[1] * (2 * (Q * Q) + a * a) + rs * (x[1] * x[1] * x[1]) - Q * Q * (Q * Q + a * a) - x[1] * x[1] * (Q * Q + rs * rs)));

    ch[1] = (1.0 / 2.0) * pown(pown(cos(x[2]), 2) * (a * a) + x[1] * x[1], -3) * pown(-rs * x[1] + Q * Q + a * a + x[1] * x[1], -1) * (-rs * v[0] * v[0] * x[1] * x[1] * (6 * (Q * Q) * (a * a) + 5 * (Q * Q * Q * Q) + a * a * a * a) - 4 * rs * v[2] * v[2] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] - rs * x[1] * x[1] * x[1] * x[1] * ((a * a) * (v[1] * v[1]) + (v[0] * v[0]) * (6 * (Q * Q) + 2 * (a * a) + rs * rs)) + rs * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) * (-v[0] * v[0] + v[1] * v[1] - 4 * v[2] * v[2] * (Q * Q + a * a)) - 1.0 / 2.0 * v[1] * v[2] * sin(2 * x[2]) * a * a * pown(cos(2 * x[2]) * (a * a) + a * a + 2 * (x[1] * x[1]), 2) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * x[1] * (Q * Q) * (v[0] * v[0]) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) + pown(sin(x[2]), 6) * (a * a * a * a) * (v[3] * v[3]) * (-rs * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) - rs * x[1] * x[1] * (6 * (Q * Q) + 6 * (a * a) + rs * rs) - 5 * rs * x[1] * x[1] * x[1] * x[1] + 2 * x[1] * (2 * (Q * Q) * (a * a) + (rs * rs) * (Q * Q + a * a) + Q * Q * Q * Q + a * a * a * a) + 4 * (x[1] * x[1] * x[1]) * (Q * Q + a * a + rs * rs) + 2 * (x[1] * x[1] * x[1] * x[1] * x[1])) + pown(sin(x[2]), 4) * (a * a) * (v[3] * v[3]) * (rs * (a * a) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) + rs * (x[1] * x[1]) * (4 * (Q * Q) * (a * a) + (a * a) * (rs * rs) - 5 * Q * Q * Q * Q + 9 * (a * a * a * a)) + rs * (x[1] * x[1] * x[1] * x[1]) * (2 * (Q * Q) + 15 * (a * a) - rs * rs) + 7 * rs * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) - 2 * x[1] * (3 * (Q * Q) * (a * a * a * a) + (a * a) * (rs * rs) * (Q * Q + a * a) - Q * Q * Q * Q * Q * Q + 2 * (a * a * a * a * a * a)) + 4 * (x[1] * x[1] * x[1]) * (-3 * Q * Q * a * a + (rs * rs) * (Q * Q - a * a) - 3 * a * a * a * a) - 2 * x[1] * x[1] * x[1] * x[1] * x[1] * (3 * (Q * Q) + 6 * (a * a) + rs * rs) - 4 * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) + pown(sin(x[2]), 2) * (2 * a * rs * v[3] * (x[1] * x[1]) * (v[0] * (6 * (Q * Q) * (a * a) + 5 * (Q * Q * Q * Q) + a * a * a * a) - 2 * v[3] * a * a * a * (Q * Q + a * a)) + a * rs * (x[1] * x[1] * x[1] * x[1]) * (a * (v[1] * v[1]) - 4 * a * v[3] * v[3] * (2 * (Q * Q) + 3 * (a * a)) + 2 * v[0] * v[3] * (6 * (Q * Q) + 2 * (a * a) + rs * rs)) - 2 * a * v[3] * x[1] * (2 * v[0] * (Q * Q) - v[3] * a * a * a) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) + 2 * a * v[3] * (x[1] * x[1] * x[1]) * (a * v[3] * (6 * (Q * Q) * (a * a) + (a * a) * (rs * rs) + 2 * (Q * Q * Q * Q) + 4 * (a * a * a * a)) - 2 * v[0] * (2 * (Q * Q) * (a * a) + (rs * rs) * (2 * (Q * Q) + a * a) + 2 * (Q * Q * Q * Q))) + 2 * rs * v[3] * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) * (a * v[0] - 2 * v[3] * (Q * Q + 3 * (a * a))) - 4 * rs * v[3] * v[3] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] + (rs - 2 * x[1]) * pown(cos(x[2]), 4) * (a * a * a * a * a * a) * (v[1] * v[1]) + 2 * pown(cos(x[2]), 2) * (a * a * a) * (-rs * v[0] * v[3] * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) - rs * v[0] * v[3] * x[1] * x[1] * x[1] * x[1] + rs * (x[1] * x[1]) * (a * (v[1] * v[1]) - v[0] * v[3] * (2 * (Q * Q) + 2 * (a * a) + rs * rs)) + 2 * v[0] * v[3] * x[1] * (rs * rs) * (Q * Q + a * a) + 2 * (x[1] * x[1] * x[1]) * (-a * v[1] * v[1] + v[0] * v[3] * (rs * rs))) + 2 * (v[3] * v[3]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) * (2 * (Q * Q) + 4 * (a * a) + rs * rs) + 2 * (v[3] * v[3]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) + 2 * (x[1] * x[1] * x[1] * x[1] * x[1]) * (-2 * a * v[0] * v[3] * (Q * Q + rs * rs) - a * a * v[1] * v[1] + (v[3] * v[3]) * (6 * (Q * Q) * (a * a) + 2 * (a * a) * (rs * rs) + Q * Q * Q * Q + 6 * (a * a * a * a)))) + pown(cos(x[2]), 4) * (a * a * a * a) * (-rs * a * a * v[1] * v[1] - 4 * rs * v[2] * v[2] * x[1] * x[1] * x[1] * x[1] + rs * (x[1] * x[1]) * (v[1] * v[1] - 4 * v[2] * v[2] * (Q * Q + a * a)) - 2 * x[1] * ((Q * Q) * (v[1] * v[1]) - v[2] * v[2] * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a)) + 2 * (v[2] * v[2]) * (x[1] * x[1] * x[1]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs) + 2 * (v[2] * v[2]) * (x[1] * x[1] * x[1] * x[1] * x[1])) + pown(cos(x[2]), 2) * (a * a) * (rs * (v[0] * v[0]) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) - 8 * rs * v[2] * v[2] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] + rs * (x[1] * x[1]) * (-2 * a * a * v[1] * v[1] + (v[0] * v[0]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs)) + rs * (x[1] * x[1] * x[1] * x[1]) * (v[0] * v[0] + 2 * (v[1] * v[1]) - 8 * v[2] * v[2] * (Q * Q + a * a)) - 2 * x[1] * rs * rs * v[0] * v[0] * (Q * Q + a * a) + 4 * (v[2] * v[2]) * (x[1] * x[1] * x[1] * x[1] * x[1]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs) + 4 * (v[2] * v[2]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) + 2 * (x[1] * x[1] * x[1]) * (-2 * Q * Q * v[1] * v[1] - rs * rs * v[0] * v[0] + 2 * (v[2] * v[2]) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a))) + 2 * (v[0] * v[0]) * (x[1] * x[1] * x[1]) * (2 * (Q * Q) * (a * a) + (rs * rs) * (2 * (Q * Q) + a * a) + 2 * (Q * Q * Q * Q)) + 2 * (v[2] * v[2]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs) + 2 * (v[2] * v[2]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) + 2 * (x[1] * x[1] * x[1] * x[1] * x[1]) * (-Q * Q * v[1] * v[1] + (v[0] * v[0]) * (Q * Q + rs * rs) + (v[2] * v[2]) * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a)));

    ch[2] = pown(pown(cos(x[2]), 2) * (a * a) + x[1] * x[1], -3) * pown(-rs * x[1] + Q * Q + a * a + x[1] * x[1], -1) * (2 * rs * v[1] * v[2] * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) + 2 * v[1] * v[2] * x[1] * pown(cos(x[2]), 4) * (a * a * a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 4 * v[1] * v[2] * pown(cos(x[2]), 2) * (a * a) * (x[1] * x[1] * x[1]) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) - 2 * v[1] * v[2] * x[1] * x[1] * x[1] * x[1] * x[1] * (Q * Q + a * a) - 2 * v[1] * v[2] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] - sin(x[2]) * cos(x[2]) * (a * rs * (x[1] * x[1] * x[1]) * (-a * v[0] * v[0] - a * v[3] * v[3] * (4 * (Q * Q) + 3 * (a * a)) + 4 * v[0] * v[3] * (Q * Q + a * a)) + a * rs * (x[1] * x[1] * x[1] * x[1] * x[1]) * (a * (v[2] * v[2]) - a * v[3] * v[3] + 2 * v[0] * v[3]) + a * (x[1] * x[1]) * (a * (v[0] * v[0]) * (Q * Q + rs * rs) + a * (v[3] * v[3]) * (3 * (Q * Q) * (a * a) + (a * a) * (rs * rs) + 2 * (Q * Q * Q * Q)) - 2 * v[0] * v[3] * (2 * (Q * Q) * (a * a) + (a * a) * (rs * rs) + Q * Q * Q * Q)) + a * (x[1] * x[1] * x[1] * x[1]) * (a * (v[1] * v[1]) - a * v[2] * v[2] * (Q * Q + a * a) + a * (v[3] * v[3]) * (Q * Q - a * a + 2 * (rs * rs)) - 2 * v[0] * v[3] * (Q * Q + rs * rs)) - rs * x[1] * a * a * (2 * (Q * Q) + a * a) * (-2 * a * v[0] * v[3] + (a * a) * (v[3] * v[3]) + v[0] * v[0]) + rs * (v[3] * v[3]) * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1]) - 2 * v[0] * v[3] * Q * Q * a * a * a * (Q * Q + a * a) + pown(cos(x[2]), 4) * (a * a * a * a) * (rs * x[1] * ((a * a) * (v[2] * v[2]) + 2 * (v[3] * v[3]) * (Q * Q + a * a)) + 2 * rs * (v[3] * v[3]) * (x[1] * x[1] * x[1]) + (a * a) * (v[1] * v[1]) - a * a * v[2] * v[2] * (Q * Q + a * a) - v[3] * v[3] * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) - v[3] * v[3] * x[1] * x[1] * x[1] * x[1] - x[1] * x[1] * ((a * a) * (v[2] * v[2]) + (v[3] * v[3]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs))) + 2 * pown(cos(x[2]), 2) * (a * a) * (x[1] * x[1]) * (rs * x[1] * ((a * a) * (v[2] * v[2]) + 2 * (v[3] * v[3]) * (Q * Q + a * a)) + 2 * rs * (v[3] * v[3]) * (x[1] * x[1] * x[1]) + (a * a) * (v[1] * v[1]) - a * a * v[2] * v[2] * (Q * Q + a * a) - v[3] * v[3] * (2 * (Q * Q) * (a * a) + Q * Q * Q * Q + a * a * a * a) - v[3] * v[3] * x[1] * x[1] * x[1] * x[1] - x[1] * x[1] * ((a * a) * (v[2] * v[2]) + (v[3] * v[3]) * (2 * (Q * Q) + 2 * (a * a) + rs * rs))) + (Q * Q) * (a * a) * (v[0] * v[0]) * (Q * Q + a * a) + (Q * Q) * (a * a * a * a) * (v[3] * v[3]) * (Q * Q + a * a) - v[3] * v[3] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * x[1] - x[1] * x[1] * x[1] * x[1] * x[1] * x[1] * ((a * a) * (v[2] * v[2]) + (v[3] * v[3]) * (Q * Q + 2 * (a * a)))));

    ch[3] = pown(pown(cos(x[2]), 2) * (a * a) + x[1] * x[1], -2) * pown(-rs * x[1] + Q * Q + a * a + x[1] * x[1], -1) * (-v[1] * sin(x[2]) * (a * rs * (-a * v[3] + v[0]) * (x[1] * x[1]) + 2 * a * x[1] * (a * v[3] - v[0]) * (Q * Q) - 2 * rs * v[3] * x[1] * x[1] * x[1] * x[1] + v[3] * (-rs + 2 * x[1]) * pown(cos(x[2]), 4) * (a * a * a * a) + 2 * v[3] * (Q * Q) * (x[1] * x[1] * x[1]) + 2 * v[3] * (x[1] * x[1] * x[1] * x[1] * x[1]) + pown(cos(x[2]), 2) * (a * a) * (-a * rs * v[0] + rs * v[3] * (a * a) - rs * v[3] * x[1] * x[1] + 4 * v[3] * (x[1] * x[1] * x[1]))) - 2 * v[2] * v[3] * pown(sin(x[2]), 4) * cos(x[2]) * a * a * a * a * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + 2 * v[2] * v[3] * pown(sin(x[2]), 2) * cos(x[2]) * (a * a) * (-rs * x[1] * (2 * (Q * Q) + 3 * (a * a)) - 3 * rs * x[1] * x[1] * x[1] + 3 * (Q * Q) * (a * a) + (x[1] * x[1]) * (3 * (Q * Q) + 4 * (a * a) + rs * rs) + Q * Q * Q * Q + 2 * (a * a * a * a) + 2 * (x[1] * x[1] * x[1] * x[1])) - 2 * v[2] * cos(x[2]) * (-a * rs * x[1] * (v[0] * (2 * (Q * Q) + a * a) + v[3] * (a * a * a)) - a * rs * (2 * a * v[3] + v[0]) * x[1] * x[1] * x[1] + a * v[0] * (Q * Q) * (Q * Q + a * a) + a * (x[1] * x[1]) * (a * v[3] * (2 * (Q * Q) + 3 * (a * a)) + v[0] * (Q * Q + rs * rs)) - rs * v[3] * x[1] * x[1] * x[1] * x[1] * x[1] + v[3] * (a * a * a * a) * (Q * Q + a * a) + v[3] * (x[1] * x[1] * x[1] * x[1]) * (Q * Q + 3 * (a * a)) + v[3] * (x[1] * x[1] * x[1] * x[1] * x[1] * x[1]))) / sin(x[2]);
    





    /*
    if (isnan(ch[0]))
    {
        printf("ch[0]=nan\n");
        //ch[0] = 0;
    }
    if (isnan(ch[1]))
    {
        printf("ch[1]=nan\n");
        //ch[1] = 0;
    }
    if (isnan(ch[2]))
    {
        printf("ch[2]=nan\n");
        //ch[2] = 0;
    }
    if (isnan(ch[3]))
    {
        printf("ch[3]=nan\n");
        //ch[3] = 0;
    }*/
}

template <>
inline __device__ void christoffel<float>(kerr_black_hole<float>& hole, float* x, float* v, float* ch)
{
    float a = hole.a;
    float Q = hole.Q;
    float rs = hole.rs;

    ch[0] = 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0f / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (a * v[1] * v[3] * (sinf(x[2])) * (sinf(x[2])) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) * (x[1] * x[1]) - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * a * a + rs * (a * a) * (x[1] * x[1]) + 3 * rs * (x[1] * x[1] * x[1] * x[1]) - 2 * x[1] * (cosf(x[2])) * (cosf(x[2])) * Q * Q * a * a - 2 * x[1] * Q * Q * a * a - 4 * Q * Q * x[1] * x[1] * x[1]) + a * v[2] * v[3] * (-rs * x[1] + Q * Q) * (-2 * sinf(x[2]) * cosf(x[2]) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) + sinf(2 * x[2]) * (a * a + x[1] * x[1]) * (rs * x[1] * (sinf(x[2])) * (sinf(x[2])) * (a * a) - (sinf(x[2])) * (sinf(x[2])) * Q * Q * a * a - (sinf(x[2])) * (sinf(x[2])) * a * a * x[1] * x[1] - (sinf(x[2])) * (sinf(x[2])) * a * a * a * a + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) + v[0] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * (sinf(x[2])) * (sinf(x[2])) * a * a * x[1] * x[1] - rs * (sinf(x[2])) * (sinf(x[2])) * a * a * a * a + rs * (a * a * a * a) - rs * x[1] * x[1] * x[1] * x[1] + 2 * x[1] * (Q * Q) * (a * a) + 2 * (Q * Q) * (x[1] * x[1] * x[1])) - v[0] * v[2] * sinf(2 * x[2]) * a * a * (-rs * x[1] + Q * Q) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]));

    ch[1] = (1.0f / 2.0f) * 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (2 * v[1] * v[2] * sinf(2 * x[2]) * (a * a) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + 2 * x[1] * (v[2] * v[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-2 * a * v[0] * v[3] * (sinf(x[2])) * (sinf(x[2])) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - (sinf(x[2])) * (sinf(x[2])) * v[3] * v[3] * (2 * x[1] * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) - ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (4 * x[1] * (a * a) + (rs - 2 * x[1]) * (sinf(x[2])) * (sinf(x[2])) * (a * a) + 4 * (x[1] * x[1] * x[1]))) + (v[0] * v[0]) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q))) - v[1] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * x[1] * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (rs - 2 * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])));

    ch[2] = (1.0f / 2.0f) * 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-4 * v[1] * v[2] * x[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) - sinf(2 * x[2]) * a * a * v[1] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) + sinf(2 * x[2]) * (a * a) * (v[2] * v[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (2 * a * v[0] * v[3] * sinf(2 * x[2]) * (-rs * x[1] + Q * Q) * (a * a + x[1] * x[1]) + 2 * sinf(x[2]) * cosf(x[2]) * (v[3] * v[3]) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) - sinf(2 * x[2]) * a * a * v[0] * v[0] * (-rs * x[1] + Q * Q)) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]));

    ch[3] = 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-2 * a * v[0] * v[2] * (-rs * x[1] + Q * Q) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + v[1] * tan(x[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (a * v[0] * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - v[3] * (-rs * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * a * a * a * a - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * x[1] * x[1] - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * a * a - rs * a * a * x[1] * x[1] + rs * (a * a * a * a) - 2 * rs * x[1] * x[1] * x[1] * x[1] + 2 * x[1] * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (a * a * a * a) + 4 * x[1] * (cosf(x[2])) * (cosf(x[2])) * (a * a * a * a) + 2 * x[1] * (Q * Q) * (a * a) - 2 * x[1] * a * a * a * a + 4 * (cosf(x[2])) * (cosf(x[2])) * (a * a) * (x[1] * x[1] * x[1]) + 2 * (Q * Q) * (x[1] * x[1] * x[1]) + 2 * (x[1] * x[1] * x[1] * x[1] * x[1]))) - v[2] * v[3] * (sinf(2 * x[2]) * tan(x[2]) * (a * a) * (-rs * x[1] + Q * Q) * (-rs * x[1] + Q * Q) * (a * a + x[1] * x[1]) + 2 * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) * (-rs * x[1] + (cosf(x[2])) * (cosf(x[2])) * (a * a) + Q * Q + x[1] * x[1]))) / tan(x[2]);

}

template <class FP>
inline __device__ FP ijk_to_vec_mink_zoom(int i, int j, int k, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg, kerr_black_hole<FP>& hole)//futtassuk le zoom nélkül és válasszunk ki egy pontot ez lesz az új kép bal felsõ sarka ezek az ikezd,jkezd// a jobb alsó sarok pedig a iveg,jveg de jveg a SZELES MAGAS arányából következik
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
inline __device__ FP ijk_to_vec_zoom(int i, int j, int k, kerr_black_hole<FP>& hole, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg)
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



    x[1] = (cos(phi) + u[0] * u[0] * (1 - cos(phi))) * x[1] + (u[0] * u[1] * (1 - cos(phi)) - u[2] * sin(phi)) * x[2] + (u[0] * u[2] * (1 - cos(phi)) + u[1] * sin(phi)) * x[3];
    x[2] = (u[0] * u[1] * (1 - cos(phi) + u[2] * sin(phi))) * x[1] + (cos(phi) + u[1] * u[1] * (1 - cos(phi))) * x[2] + (u[1] * u[2] * (1 - cos(phi) + u[0] * sin(phi))) * x[3];
    x[3] = (u[0] * u[2] * (1 - cos(phi)) - u[1] * sin(phi)) * x[1] + (u[1] * u[2] * (1 - cos(phi)) + u[0] * sin(phi)) * x[2] + (cos(phi) + u[2] * u[2] * (1 - cos(phi))) * x[3];


    FP r_0 = hole.r_0;
    FP theta_0 = hole.theta_0;
    FP rs = hole.rs;
    FP a = hole.a;
    FP Q = hole.Q;

    FP delta = r_0 * r_0 - 4 * rs * r_0 + a * a + Q * Q;
    FP rho = sqrt(r_0 * r_0 + a * a * cos(theta_0) * cos(theta_0));

    x[0] = x[0] * (a * a + r_0 * r_0) * rho / ((a * a * cos(theta_0) * cos(theta_0) + r_0 * r_0) * sqrt(delta)) + x[3] * a * rho / (sqrt(delta) * (a * a * cos(theta_0) * cos(theta_0) + r_0 * r_0));
    x[1] = sqrt(delta) / rho * x[1];
    x[2] = x[2] / rho;
    x[3] = x[0] * a * rho * sin(theta_0) / (a * a * cos(theta_0) * cos(theta_0) + r_0 * r_0) + x[3] * rho / (sin(theta_0) * (r_0 * r_0 + a * a * cos(theta_0) * cos(theta_0)));



    return x[k];
}

template <class FP>
inline __global__ void ray_step(int8_t* szin, int SZELES, int MAGAS, FP* xd, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg)//kernel
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < SZELES && j < MAGAS)
    {

        kerr_black_hole<FP> hole(SZELES, MAGAS, xd, Omega, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy);
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

template <class FP>
inline __global__ void ray_step_T(FP* szin, int SZELES, int MAGAS, FP* xd, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg)//kernel
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < SZELES && j < MAGAS)
    {

        kerr_black_hole<FP> hole(SZELES, MAGAS, xd, Omega, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy);
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
        {//0 fekete, -1 hiba kezeléses piros, egyébként meg egy FP ami reprezental egy szint
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

//ha a gömbön belül van akkor igaz
template <class FP>
inline __device__ bool gomb_be(FP sugar, FP* x)
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
inline __device__ bool gomb_ki(FP sugar, FP* x)
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
inline __device__ bool disk(FP sugar_kicsi, FP sugar_nagy, FP* x1, FP* x2)
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
inline __device__ bool disk1(FP sugar_kicsi, FP sugar_nagy, FP* x1, FP* x2)
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
inline __device__ bool disk2(FP sugar_kicsi, FP sugar_nagy, FP* x1, FP* x2)
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
inline __host__ __device__ int ijk_to_n(int i, int j, int k, kerr_black_hole<FP>& hole)
{
    return i * hole.MAGAS * D + j * D + k;
}

template <class FP>
inline __device__ FP pown(FP x, int n)
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
