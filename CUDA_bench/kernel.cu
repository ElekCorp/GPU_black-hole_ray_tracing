#include <iostream>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//#include <fstream>

#include <string>

#include <chrono>

#include <math.h>

//#include <SDL2/SDL.h>
//#include <SDL2/SDL_ttf.h>

#include "black_hole.cuh"
#include "cuda_ray.cuh"

#include "szinsaver.h"

//#include "debugmalloc.h"

template <class FP>
int8_t* makeframe(int SZELES, int MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg);

template <class FP>
FP* makeframe_T(int SZELES, int MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg);

void device_info(void);
//void SDL_szinezo(SDL_Renderer* ren, int8_t* SZIN, int i, int j, int ki, int kj, int MAGAS);
//void fps_kiiro(int fps);

int n_oszto(int SZELES, int MAGAS, int kepernyoSZELES, int kepernyoMAGAS, int n);


int main(int argc, char* argv[])
{
    device_info();

    //auto start = std::chrono::high_resolution_clock::now();

    int kepernyoSZELES = 1024;
    int kepernyoMAGAS = 512;

    int SZELES = 128;
    int MAGAS = 64;

    if (kepernyoSZELES * MAGAS != kepernyoMAGAS * SZELES)
    {
        std::cout << "renderelendo kep es a kepernyo aranya nem azonos\n";
        return 1;
    }

    float errormax = 0.001f;
    float de0 = 0.01f;

    float rs = 0.05f;//2*rs=m

    float delta_a = 0.0001f;

    float a = 0.0f;//rs / 2 -delta_a;//0.0;
    float Q = 0.0f;



    float t_0 = 0.0f;
    float r_0 = 1.0f;
    float theta_0 = 1.57f + 0.06f;//ne legyen nulla
    float phi_0 = 0.0f;

    float kepernyo_high = 0.5f;
    float kepernyo_tav = 0.75f;//0.4;//0.75

    float sugar_ki = 1.01f;

    float gyuru_sugar_kicsi = 0.1f;
    float gyuru_sugar_nagy = 0.5f;

    float x[D] = { t_0,r_0,theta_0,phi_0 };
    float pi_cucc = float(asin(1) * 2);
    float Omega[D - 1] = { 0,pi_cucc,0 };




    int SZELESregi = SZELES;
    int MAGASregi = MAGAS;

    int ikezd = 0;
    int jkezd = 0;
    int iveg = SZELES;


    double* SZIN = NULL;


	double errormax_d = double(errormax);
    double de0_d = double(de0);

    double rs_d = double(rs);//2*rs=m

    double delta_a_d = double(delta_a);

    double a_d = double(a);//rs / 2 -delta_a;//0.0;
    double Q_d = double(Q);

    double kepernyo_high_d = double(kepernyo_high);
    double kepernyo_tav_d = double(kepernyo_tav);//0.4;//0.75

    double sugar_ki_d = double(sugar_ki);

    double gyuru_sugar_kicsi_d = double(gyuru_sugar_kicsi);
    double gyuru_sugar_nagy_d = double(gyuru_sugar_nagy);

    double x_d[D] = { double(x[0]),double(x[1]),double(x[2]),double(x[3]) };
    double Omega_d[D - 1] = { double(Omega[0]),double(Omega[1]),double(Omega[2]) };


    

    for(int i=0;i<100;++i)
    {
        a_d+=0.01*rs_d/2;

        SZIN = makeframe_T<double>(SZELES, MAGAS, x_d, Omega_d, a_d, Q_d, rs_d, errormax_d, de0_d, kepernyo_high_d, kepernyo_tav_d, sugar_ki_d, gyuru_sugar_kicsi_d, gyuru_sugar_nagy_d, SZELESregi, MAGASregi, ikezd, jkezd, iveg);

        std::string kep_double_string="./gif_material/kep_double"+std::to_string(i)+".dat";
        datasaver_T<double>(SZIN, SZELES, MAGAS, kep_double_string);
    
        free(SZIN);
    }

    return 0;
}

template <class FP>
  FP* makeframe_T(int SZELES, int MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, int SZELESregi, int MAGASregi, int ikezd, int jkezd, int iveg)//ekkor a SZIN egy FP* es a homersekletet reprezentalja
  {
      std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

      FP* x_d = NULL;
      FP* Omega_d = NULL;
      cudaError_t err = cudaSuccess;

      //auto start = std::chrono::high_resolution_clock::now();

      err = cudaMalloc(&x_d, D * sizeof(FP));
      if (err != cudaSuccess)
      {
          fprintf(stderr, "\n cudaMalloc hiba x_d foglalasakor %s\n", cudaGetErrorString(err));
          exit(EXIT_FAILURE);
      }
      err = cudaMalloc(&Omega_d, (D - 1) * sizeof(FP));
      if (err != cudaSuccess)
      {
          fprintf(stderr, "\n cudaMalloc hiba Omega_d foglalasakor %s\n", cudaGetErrorString(err));
          exit(EXIT_FAILURE);
      }

      err = cudaMemcpy(x_d, x, D * sizeof(FP), cudaMemcpyHostToDevice);
      err = cudaMemcpy(Omega_d, Omega, (D - 1) * sizeof(FP), cudaMemcpyHostToDevice);

     FP* SZIN_d = NULL;

      err = cudaMalloc(&SZIN_d, size_t(SZELES) * size_t(MAGAS) * sizeof(FP));
      if (err != cudaSuccess)
      {
          fprintf(stderr, "\n cudaMalloc hiba SZIN_d foglalasakor %s\n", cudaGetErrorString(err));
          exit(EXIT_FAILURE);
      }

      dim3 threadsPerBlock(16, 16);
      int xdim = SZELES / threadsPerBlock.x, ydim = MAGAS / threadsPerBlock.y;


      if (SZELES % threadsPerBlock.x != 0)
      {
          xdim = (SZELES + threadsPerBlock.x - 1) / threadsPerBlock.x;
      }
      if (MAGAS % threadsPerBlock.y != 0)
      {
          ydim = (MAGAS + threadsPerBlock.y - 1) / threadsPerBlock.y;
      }

      dim3 numBlocks(xdim, ydim);

      ray_step_T << <numBlocks, threadsPerBlock >> > (SZIN_d, SZELES, MAGAS, x_d, Omega_d, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy, SZELESregi, MGASregi, ikezd, jkezd, iveg);

      cudaDeviceSynchronize();

      err = cudaPeekAtLastError();
      if (err != cudaSuccess)
      {
          fprintf(stderr, "\nkernel:%s\n", cudaGetErrorString(err));
          exit(EXIT_FAILURE);
      }

      FP* SZIN = NULL;
      SZIN = (FP*)malloc(SZELES * MAGAS * sizeof(FP));

      err = cudaMemcpy(SZIN, SZIN_d, SZELES * MAGAS * sizeof(FP), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess)
      {
          fprintf(stderr, "\ncudaMemcpy SZIN,SZIN_d, DeviceToHost nal:%s\n", cudaGetErrorString(err));
          exit(EXIT_FAILURE);
      }

      cudaFree(SZIN_d);
      cudaFree(x_d);
      cudaFree(Omega_d);

      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

      std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 1000000.0 << "[s]" << std::endl;

      return SZIN;
}

void device_info(void)
{
    int num_gpus = 0;
    cudaGetDeviceCount(&num_gpus);

    if (num_gpus == 0)
    {
        std::cout << "no capable GPU\n";
        exit(EXIT_FAILURE);
    }

    //printf("number of host CPUs:\t%d\n", omp_get_num_procs());
    printf("number of CUDA devices:\t%d\n", num_gpus);

    for (int i = 0; i < num_gpus; i++)
    {
        cudaDeviceProp dprop;
        cudaGetDeviceProperties(&dprop, i);
        printf("%d: %s\n", i, dprop.name);
        printf("%d: L2 Cache Size:%d bytes\n", i, dprop.l2CacheSize);


        int driverVersion = 0, runtimeVersion = 0;

        cudaDriverGetVersion(&driverVersion);
        cudaRuntimeGetVersion(&runtimeVersion);
        printf("  CUDA Driver Version / Runtime Version          %d.%d / %d.%d\n",
               driverVersion / 1000, (driverVersion % 100) / 10,
               runtimeVersion / 1000, (runtimeVersion % 100) / 10);
        printf("  CUDA Capability Major/Minor version number:    %d.%d\n",
               dprop.major, dprop.minor);

    }
}

int n_oszto(int SZELES, int MAGAS, int kepernyoSZELES, int kepernyoMAGAS, int n)
{
    int oszto = 1;
    int kis_kep=kepernyoMAGAS;
    int num = 1;
    if (kepernyoMAGAS > kepernyoSZELES)
    {
        kis_kep = kepernyoSZELES;
    }
    while (num < n)
    {
        ++oszto;
        if (kis_kep % oszto == 0)
        {
            ++num;
        }

    }

    if ((kepernyoMAGAS % oszto != 0) || (kepernyoSZELES % oszto != 0))
    {
        std::cout << oszto << "A kepernyoMAGAS vagy kepernyoSZELES nem oszthato a viszatert oszto-val\n";
    }
    std::cout << oszto << "\n";


    return oszto;
}
