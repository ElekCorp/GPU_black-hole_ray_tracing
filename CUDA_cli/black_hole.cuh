#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <cstdio>
//#include <cmath>

#include <iostream>

// Spacetime dimension; x and v are always D-vectors.
#define D 4

/**
 * Every parameter one ray needs, bundled so it can be built once per thread.
 *
 * Geometry (geometric units, G = c = 1):
 *   rs                  Schwarzschild radius, rs = 2M
 *   a                   spin parameter, a = J/M
 *   Q                   electric charge
 *   Sigma = r^2 + a^2 cos^2(theta),  Delta = r^2 - rs*r + a^2 + Q^2
 *
 * Camera:
 *   t_0,r_0,theta_0,phi_0   Boyer-Lindquist position.  theta_0 must not be 0
 *                           or pi (the metric is singular on the axis) and r_0
 *                           must be outside the outer horizon, i.e. Delta > 0 -
 *                           the tetrad in ijk_to_vec_zoom() takes sqrt(Delta).
 *   Omega_1..3          orientation as a rotation vector; direction is the
 *                       axis, length is the angle in radians
 *   kepernyo_tav        distance from pinhole to image plane
 *   kepernyo_high       total height of the image plane
 *
 * Integration:
 *   de0                 largest affine step allowed
 *   errormax            per-step local error tolerance for the DOPRI5 pair
 *   max_steps           step budget per ray; a ray that exceeds it is marked
 *                       as a failure.  This used to be derived from errormax
 *                       as int(1/errormax), which both conflated two unrelated
 *                       knobs (tightening the tolerance shrank the budget) and
 *                       overflowed int for errormax below about 1e-9.
 *
 * Scene:
 *   sugar_ki            sky sphere; a ray beyond it has escaped
 *   sugar_kicsi/nagy    inner/outer radius of the equatorial accretion disk
 */
template <class FP>
class kerr_black_hole
{
public:

	uint64_t SZELES;// 512//1024//1808//1024//8192
	uint64_t MAGAS;//512//608//512//8192 esetén még mûködik


	FP errormax;
	FP de0;
	uint64_t max_steps;

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




   

	__device__ void iro(void)
	{
		printf("%llu,%llu,%f", (unsigned long long)SZELES, (unsigned long long)MAGAS, double(de0));
	}
	__device__ kerr_black_hole(uint64_t SZELES, uint64_t MAGAS, FP const* x, FP const* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, uint64_t max_steps, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy)
	{
		this->SZELES = SZELES;
		this->MAGAS = MAGAS;//512//608//512//8192 esetén még mûködik


		this->errormax = errormax;
		this->de0 = de0;
		this->max_steps = max_steps;

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

	/*__host__ __device__ kerr_black_hole(const kerr_black_hole& x)
	{
		
	}*/


	__device__ ~kerr_black_hole()
	{
		
	}

};

