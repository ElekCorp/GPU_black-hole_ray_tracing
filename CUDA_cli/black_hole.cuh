#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <cstdio>
//#include <cmath>

#include <iostream>

#define D 4

template <class FP>
class kerr_black_hole
{
public:

	int SZELES;// 512//1024//1808//1024//8192
	int MAGAS;//512//608//512//8192 esetén még mûködik


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




   

	__device__ void iro(void)
	{
		printf("%d,%d,%f", SZELES, MAGAS, de0);
	}
	__device__ kerr_black_hole(int SZELES, int MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy)
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

	/*__host__ __device__ kerr_black_hole(const kerr_black_hole& x)
	{
		
	}*/


	__device__ ~kerr_black_hole()
	{
		
	}

};

