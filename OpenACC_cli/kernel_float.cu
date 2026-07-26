
typedef float FP;
typedef unsigned long long uint64_t;
typedef signed char int8_t;
#ifndef BLACK_HOLE_H
#define BLACK_HOLE_H

//#include <math.h>


//#include "debugmalloc.h"

#define D 4


class kerr_black_hole
{
public:

	uint64_t SZELES;// 512//1024//1808//1024//8192
	uint64_t MAGAS;//512//608//512//8192 esetén még mûködik


	FP errormax;
	FP de0;

	FP rs;//rs=2*m


	FP a;
	FP Q;



	FP t_0;
	FP r_0;
	FP theta_0;//ne legyen nulla
	FP phi_0;

	FP kepernyo_high;
	FP kepernyo_tav;//0.4;//0.75
	FP arany;

	FP sugar_ki;

	FP sugar_kicsi;
	FP sugar_nagy;

	FP Omega_1;
	FP Omega_2;
	FP Omega_3;






	

	__device__ kerr_black_hole(uint64_t const SZELES, uint64_t const MAGAS, FP const* const x, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy)
	{
		this->SZELES = SZELES;
		this->MAGAS = MAGAS;//512//608//512//8192 esetén még mûködik


		this->errormax = errormax;
		this->de0 = de0;

		this->rs = rs;//rs=2*m


		this->a = a;
		this->Q = Q;



		this->t_0 = x[0];
		this->r_0 = x[1];
		this->theta_0 = x[2];//ne legyen nulla
		this->phi_0 = x[3];

		this->kepernyo_high = kepernyo_high;
		this->kepernyo_tav = kepernyo_tav;//0.4;//0.75
		this->arany = kepernyo_high / MAGAS;

		this->sugar_ki = sugar_ki;

		this->sugar_kicsi = gyuru_sugar_kicsi;
		this->sugar_nagy = gyuru_sugar_nagy;

		this->Omega_1 = Omega[0];
		this->Omega_2 = Omega[1];
		this->Omega_3 = Omega[2];

	}

	/*__host__ __device__ __device__ kerr_black_hole(const kerr_black_hole& x)
	{

	}*/


	/*~__device__ kerr_black_hole()
	{

	}*/

};


#endif // BLACK_HOLE_H
#ifndef CUDA_RAY_H
#define CUDA_RAY_H



//#include "debugmalloc.h"



__device__ inline void step(kerr_black_hole const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de);
// step_size removed


__device__ inline FP ijk_to_vec_mink_zoom(uint64_t const i, uint64_t const j, uint64_t const k, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg, kerr_black_hole const& hole);

__device__ inline FP ijk_to_vec_zoom(uint64_t const i, uint64_t const j, uint64_t const k, kerr_black_hole const& hole, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg);




__device__ inline void RK38(kerr_black_hole const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de);//4-ed foku legpontosabb

__device__ inline FP RKDP45_core(kerr_black_hole const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ x_new, FP* const __restrict__ v_new, FP const de);


__device__ inline void christoffel(kerr_black_hole const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ ch);


__device__ inline void ray_step(int8_t* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const x, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg);

__device__ inline void ray_step_T(FP* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const x, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg);





__device__ inline bool gomb_be(FP const sugar, FP const* const x);

__device__ inline bool gomb_ki(FP const sugar, FP const* const x);

__device__ inline bool disk(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2);


__device__ inline bool disk1(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2);

__device__ inline bool disk2(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2);


__device__ inline int ijk_to_n(uint64_t const i, uint64_t const j, uint64_t const k, kerr_black_hole const& hole);


__device__ inline FP pown(FP const x, int const n);

__device__ inline FP pown_gen(FP const x, int const n);

__device__ inline FP pown_rec(FP const x, int const n);



__device__ inline void step(kerr_black_hole const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)
{
    bool accepted = false;
    FP x_new[D];
    FP v_new[D];
    FP de_new = de;
    int max_retries = 10;
    int retries = 0;

    while (!accepted && retries < max_retries)
    {
        FP err = RKDP45_core(hole, x, v, x_new, v_new, de_new);
        
        if (err < 1e-20) err = 1e-20;
        
        FP factor = 0.9 * pow(hole.errormax / err, FP(0.2));
        
        if (factor > 5.0) factor = 5.0;
        if (factor < 0.1) factor = 0.1;

        if (err <= hole.errormax) {
            accepted = true;
        }

        FP de_next = de_new * factor;
        if (de_next > hole.de0) de_next = hole.de0;

        if (accepted) {
            for (int i = 0; i < D; ++i) {
                x[i] = x_new[i];
                v[i] = v_new[i];
            }
            de = de_next; 
        } else {
            de_new = de_next; 
            retries++;
        }
    }
    
    if (!accepted) {
        for (int i = 0; i < D; ++i) {
            x[i] = x_new[i];
            v[i] = v_new[i];
        }
        de = de_new;
    }
}




__device__ inline void RK38(kerr_black_hole const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)
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


__device__ inline FP RKDP45_core(kerr_black_hole const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ x_new, FP* const __restrict__ v_new, FP const de)
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



    FP err = 0.0;
    for (int i = 0; i < D; ++i)
    {
        x_new[i] = x[i] + (kx1[i] * (35.0 / 384.0) + kx3[i] * (500.0 / 1113.0) + kx4[i] * (125.0 / 192.0) - kx5[i] * (2187.0 / 6784.0) + kx6[i] * (11.0 / 84.0)) * de;
        v_new[i] = v[i] + (kv1[i] * (35.0 / 384.0) + kv3[i] * (500.0 / 1113.0) + kv4[i] * (125.0 / 192.0) - kv5[i] * (2187.0 / 6784.0) + kv6[i] * (11.0 / 84.0)) * de;
        
        FP x_4th = x[i] + (kx1[i] * (5179.0 / 57600.0) + kx3[i] * (7571.0 / 16695.0) + kx4[i] * (393.0 / 640.0) - kx5[i] * (92097.0 / 339200.0) + kx6[i] * (187.0 / 2100.0) + kx7[i] * (1.0 / 40.0)) * de;
        FP v_4th = v[i] + (kv1[i] * (5179.0 / 57600.0) + kv3[i] * (7571.0 / 16695.0) + kv4[i] * (393.0 / 640.0) - kv5[i] * (92097.0 / 339200.0) + kv6[i] * (187.0 / 2100.0) + kv7[i] * (1.0 / 40.0)) * de;

        FP diff_x = x_new[i] - x_4th;
        FP diff_v = v_new[i] - v_4th;
        err += diff_x * diff_x + diff_v * diff_v;
    }
    return sqrt(err);
}


__device__ inline void christoffel(kerr_black_hole const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ ch)
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
__device__ inline void christoffel<float>(kerr_black_hole<float>& hole, float* x, float* v, float* ch)
{
    float a = hole.a;
    float Q = hole.Q;
    float rs = hole.rs;

    ch[0] = 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0f / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (a * v[1] * v[3] * (sinf(x[2])) * (sinf(x[2])) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) * (x[1] * x[1]) - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * a * a + rs * (a * a) * (x[1] * x[1]) + 3 * rs * (x[1] * x[1] * x[1] * x[1]) - 2 * x[1] * (cosf(x[2])) * (cosf(x[2])) * Q * Q * a * a - 2 * x[1] * Q * Q * a * a - 4 * Q * Q * x[1] * x[1] * x[1]) + a * v[2] * v[3] * (-rs * x[1] + Q * Q) * (-2 * sinf(x[2]) * cosf(x[2]) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) + sinf(2 * x[2]) * (a * a + x[1] * x[1]) * (rs * x[1] * (sinf(x[2])) * (sinf(x[2])) * (a * a) - (sinf(x[2])) * (sinf(x[2])) * Q * Q * a * a - (sinf(x[2])) * (sinf(x[2])) * a * a * x[1] * x[1] - (sinf(x[2])) * (sinf(x[2])) * a * a * a * a + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) + v[0] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * (sinf(x[2])) * (sinf(x[2])) * a * a * x[1] * x[1] - rs * (sinf(x[2])) * (sinf(x[2])) * a * a * a * a + rs * (a * a * a * a) - rs * x[1] * x[1] * x[1] * x[1] + 2 * x[1] * (Q * Q) * (a * a) + 2 * (Q * Q) * (x[1] * x[1] * x[1])) - v[0] * v[2] * sinf(2 * x[2]) * a * a * (-rs * x[1] + Q * Q) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]));

    ch[1] = (1.0f / 2.0f) * 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (2 * v[1] * v[2] * sinf(2 * x[2]) * (a * a) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + 2 * x[1] * (v[2] * v[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-2 * a * v[0] * v[3] * (sinf(x[2])) * (sinf(x[2])) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - (sinf(x[2])) * (sinf(x[2])) * v[3] * v[3] * (2 * x[1] * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) - ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (4 * x[1] * (a * a) + (rs - 2 * x[1]) * (sinf(x[2])) * (sinf(x[2])) * (a * a) + 4 * (x[1] * x[1] * x[1]))) + (v[0] * v[0]) * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q))) - v[1] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * x[1] * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (rs - 2 * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])));

    ch[2] = (1.0f / 2.0f) * 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-4 * v[1] * v[2] * x[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) - sinf(2 * x[2]) * a * a * v[1] * v[1] * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) + sinf(2 * x[2]) * (a * a) * (v[2] * v[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + (2 * a * v[0] * v[3] * sinf(2 * x[2]) * (-rs * x[1] + Q * Q) * (a * a + x[1] * x[1]) + 2 * sinf(x[2]) * cosf(x[2]) * (v[3] * v[3]) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) - sinf(2 * x[2]) * a * a * v[0] * v[0] * (-rs * x[1] + Q * Q)) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]));

    ch[3] = 1.0f / (((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1])) * 1.0 / (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) * (-2 * a * v[0] * v[2] * (-rs * x[1] + Q * Q) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (-rs * x[1] + Q * Q + a * a + x[1] * x[1]) + v[1] * tan(x[2]) * ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (a * v[0] * (rs * (cosf(x[2])) * (cosf(x[2])) * (a * a) - rs * x[1] * x[1] + 2 * x[1] * (Q * Q)) - v[3] * (-rs * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * a * a * a * a - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * x[1] * x[1] - rs * (cosf(x[2])) * (cosf(x[2])) * a * a * a * a - rs * a * a * x[1] * x[1] + rs * (a * a * a * a) - 2 * rs * x[1] * x[1] * x[1] * x[1] + 2 * x[1] * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (sinf(x[2])) * (a * a * a * a) + 4 * x[1] * (cosf(x[2])) * (cosf(x[2])) * (a * a * a * a) + 2 * x[1] * (Q * Q) * (a * a) - 2 * x[1] * a * a * a * a + 4 * (cosf(x[2])) * (cosf(x[2])) * (a * a) * (x[1] * x[1] * x[1]) + 2 * (Q * Q) * (x[1] * x[1] * x[1]) + 2 * (x[1] * x[1] * x[1] * x[1] * x[1]))) - v[2] * v[3] * (sinf(2 * x[2]) * tan(x[2]) * (a * a) * (-rs * x[1] + Q * Q) * (-rs * x[1] + Q * Q) * (a * a + x[1] * x[1]) + 2 * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * ((sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1]) + ((cosf(x[2])) * (cosf(x[2])) * (a * a) + x[1] * x[1]) * (2 * (sinf(x[2])) * (sinf(x[2])) * (a * a) * (rs * x[1] - Q * Q - a * a - x[1] * x[1]) + 2 * (a * a) * (x[1] * x[1]) + a * a * a * a + x[1] * x[1] * x[1] * x[1])) * (-rs * x[1] + (cosf(x[2])) * (cosf(x[2])) * (a * a) + Q * Q + x[1] * x[1]))) / tan(x[2]);

}*/


__device__ inline FP ijk_to_vec_mink_zoom(uint64_t i, uint64_t j, uint64_t k, uint64_t SZELESregi, uint64_t MAGASregi, uint64_t ikezd, uint64_t jkezd, uint64_t iveg, kerr_black_hole& hole)//futtassuk le zoom nélkül és válasszunk ki egy pontot ez lesz az új kép bal felsõ sarka ezek az ikezd,jkezd// a jobb alsó sarok pedig a iveg,jveg de jveg a SZELES MAGAS arányából következik
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


__device__ inline FP ijk_to_vec_zoom(uint64_t i, uint64_t j, uint64_t k, kerr_black_hole& hole, uint64_t SZELESregi, uint64_t MAGASregi, uint64_t ikezd, uint64_t jkezd, uint64_t iveg)
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


//ha a gömbön belül van akkor igaz

__device__ inline bool gomb_be(FP const sugar, FP const* const x)
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

__device__ inline bool gomb_ki(FP const sugar, FP const* const x)
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

__device__ inline bool disk(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2)
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


__device__ inline bool disk1(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2)
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


__device__ inline bool disk2(FP const sugar_kicsi, FP const sugar_nagy, FP const* const __restrict__ x1, FP const* const __restrict__ x2)
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


__device__ inline int ijk_to_n(uint64_t const i, uint64_t const j, uint64_t const k, kerr_black_hole const& hole)
{
    return i * hole.MAGAS * D + j * D + k;
}


__device__ inline FP pown(FP const x, int const n)
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

__device__ inline FP pown_gen(FP const x, int const n) 
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


__device__ inline FP powni_rec(FP const x, int const n)
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

extern "C" __global__ void ray_step_kernel(int8_t* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const xd, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki_in, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg) {
    uint64_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= SZELES || j >= MAGAS) return;

    kerr_black_hole hole(SZELES, MAGAS, xd, Omega, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki_in, gyuru_sugar_kicsi, gyuru_sugar_nagy);

    FP x[D] = { hole.t_0, hole.r_0, hole.theta_0, hole.phi_0 };
    FP v[D];
    FP x_le[D];
    FP de = de0;

    for (int k = 0; k < D; ++k) {
        v[k] = ijk_to_vec_zoom(i, j, k, hole, SZELESregi, MAGASregi, ikezd, jkezd, iveg);
    }

    for (int k = 0; k < D; ++k) {
        x_le[k] = x[k];
    }
    step(hole, x, v, de);

    bool fut = true;
    FP sugar_be = 0;

    if ((rs * rs - 4 * (a * a + Q * Q)) > 0.0) {
        sugar_be = (rs + sqrt(rs * rs - 4 * (a * a + Q * Q))) / 2 + 0.001;
    }

    FP sugar_ki = hole.sugar_ki;
    FP sugar_kicsi = hole.sugar_kicsi;
    FP sugar_nagy = hole.sugar_nagy;
    int idokorlat = 0;

    while (fut) {
        if (gomb_be(sugar_be, x)) {
            szin[i * MAGAS + j] = 1;
            fut = false;
        } else if (gomb_ki(sugar_ki, x)) {
            szin[i * MAGAS + j] = 1;
            fut = false;
        } else if (disk1(sugar_kicsi, sugar_nagy, x, x_le)) {
            szin[i * MAGAS + j] = 0;
            fut = false;
        } else if (disk2(sugar_kicsi, sugar_nagy, x, x_le)) {
            szin[i * MAGAS + j] = 3;
            fut = false;
        } else if (isnan(x[0]) || isnan(x[1]) || isnan(x[2]) || isnan(x[3])) {
            szin[i * MAGAS + j] = 2;
            fut = false;
        } else if (isinf(x[0]) || isinf(x[1]) || isinf(x[2]) || isinf(x[3])) {
            szin[i * MAGAS + j] = 2;
            fut = false;
        }

        ++idokorlat;
        if (idokorlat >= int(20.0 / errormax)) {
            szin[i * MAGAS + j] = -1;
            fut = false;
        }

        if (!fut) break;

        for (int k = 0; k < D; ++k) {
            x_le[k] = x[k];
        }
        step(hole, x, v, de);
    }
}

extern "C" __global__ void ray_step_kernel_T(FP* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const xd, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki_in, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg) {
    uint64_t i = blockIdx.x * blockDim.x + threadIdx.x;
    uint64_t j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= SZELES || j >= MAGAS) return;

    kerr_black_hole hole(SZELES, MAGAS, xd, Omega, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki_in, gyuru_sugar_kicsi, gyuru_sugar_nagy);

    FP x[D] = { hole.t_0, hole.r_0, hole.theta_0, hole.phi_0 };
    FP v[D];
    FP x_le[D];
    FP de = de0;

    for (int k = 0; k < D; ++k) {
        v[k] = ijk_to_vec_zoom(i, j, k, hole, SZELESregi, MAGASregi, ikezd, jkezd, iveg);
    }

    for (int k = 0; k < D; ++k) {
        x_le[k] = x[k];
    }
    step(hole, x, v, de);

    bool fut = true;
    FP sugar_be = 0;

    if ((rs * rs - 4 * (a * a + Q * Q)) > 0.0) {
        sugar_be = (rs + sqrt(rs * rs - 4 * (a * a + Q * Q))) / 2 + 0.001;
    }

    FP sugar_ki = hole.sugar_ki;
    FP sugar_kicsi = hole.sugar_kicsi;
    FP sugar_nagy = hole.sugar_nagy;
    int idokorlat = 0;

    while (fut) {
        if (gomb_be(sugar_be, x)) {
            szin[i * MAGAS + j] = 0;
            fut = false;
        } else if (gomb_ki(sugar_ki, x)) {
            szin[i * MAGAS + j] = 0;
            fut = false;
        } else if (disk1(sugar_kicsi, sugar_nagy, x, x_le)) {
            szin[i * MAGAS + j] = x[1];
            fut = false;
        } else if (disk2(sugar_kicsi, sugar_nagy, x, x_le)) {
            szin[i * MAGAS + j] = x[1];
            fut = false;
        } else if (isnan(x[0]) || isnan(x[1]) || isnan(x[2]) || isnan(x[3])) {
            szin[i * MAGAS + j] = -1;
            fut = false;
        } else if (isinf(x[0]) || isinf(x[1]) || isinf(x[2]) || isinf(x[3])) {
            szin[i * MAGAS + j] = -1;
            fut = false;
        }

        ++idokorlat;
        if (idokorlat >= int(20.0 / errormax)) {
            szin[i * MAGAS + j] = -1;
            fut = false;
        }

        if (!fut) break;

        for (int k = 0; k < D; ++k) {
            x_le[k] = x[k];
        }
        step(hole, x, v, de);
    }
}
