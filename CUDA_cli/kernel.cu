#include <iostream>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <string>

#include <chrono>

#include <math.h>
#include <vector_types.h>

#include "black_hole.cuh"
#include "cuda_ray.cuh"

#include "szinsaver.h"

#include "cli_parser.h"

//#include "debugmalloc.h"

template <class FP>
int8_t* makeframe(uint64_t SZELES, uint64_t MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, uint64_t SZELESregi, uint64_t MAGASregi, uint64_t ikezd, uint64_t jkezd, uint64_t iveg);

template <class FP>
FP* makeframe_T(uint64_t SZELES, uint64_t MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, uint64_t SZELESregi, uint64_t MAGASregi, uint64_t ikezd, uint64_t jkezd, uint64_t iveg);

void device_info(void);
//void SDL_szinezo(SDL_Renderer* ren, int8_t* SZIN, int i, int j, int ki, int kj, int MAGAS);
//void fps_kiiro(int fps);

uint64_t n_oszto(uint64_t SZELES, uint64_t MAGAS, uint64_t kepernyoSZELES, uint64_t kepernyoMAGAS, uint64_t n);


int main(int argc, char* argv[])
{
    device_info();
    
    Params p;
    parse_args(argc,argv,p);
    //auto start = std::chrono::high_resolution_clock::now();

    if (p.kepernyoSZELES * p.MAGAS != p.kepernyoMAGAS * p.SZELES)
    {
        std::cout << "renderelendo kep es a kepernyo aranya nem azonos\n";
        return 1;
    }

    float x[D] = { float(p.t_0),float(p.r_0),float(p.theta_0),float(p.phi_0) };


    double pi_cucc = (asin(1) * 2);
    float Omega[D - 1] = { 0,float(pi_cucc),0 };




    uint64_t SZELESregi = p.kepernyoSZELES;
    uint64_t MAGASregi = p.kepernyoMAGAS;

    uint64_t ikezd = p.ikezd;
    uint64_t jkezd = p.jkezd;
    uint64_t iveg = p.iveg;

    double* SZIN = NULL;
    float* SZIN_f=NULL;


    double x_d[D] = { p.t_0, p.r_0, p.theta_0, p.phi_0 };
    double Omega_d[D - 1] = { double(Omega[0]),double(Omega[1]),double(Omega[2]) };
    
    std::cout<<p.iveg<<"/n"<<p.ikezd<<"/n"<<p.jkezd<<std::endl;

if(p.prec==Precession::Double)
{
    SZIN = makeframe_T<double>(p.SZELES, p.MAGAS, x_d, Omega_d, p.a, p.Q, p.rs, p.errormax, p.de0, p.kepernyo_high, p.kepernyo_tav, p.sugar_ki, p.gyuru_sugar_kicsi, p.gyuru_sugar_nagy, SZELESregi, MAGASregi, ikezd, jkezd, iveg);
    std::string kep_double_string="./web_images/kep_cli.dat";
    datasaver_T<double>(SZIN, p.SZELES, p.MAGAS, kep_double_string);

    free(SZIN);
}
else
{
    SZIN_f = makeframe_T<float>(p.SZELES, p.MAGAS, x, Omega,float(p.a),float(p.Q), float(p.rs), float(p.errormax),float(p.de0), float(p.kepernyo_high), float(p.kepernyo_tav), float(p.sugar_ki), float(p.gyuru_sugar_kicsi), float(p.gyuru_sugar_nagy), SZELESregi, MAGASregi, ikezd, jkezd, iveg);
    std::string kep_string="./web_images/kep_cli.dat";
    datasaver_T<float>(SZIN_f, p.SZELES, p.MAGAS, kep_string);

    free(SZIN_f);
}


    return 0;
}

template <class FP>
FP* makeframe_T(uint64_t SZELES, uint64_t MAGAS, FP* x, FP* Omega, FP a, FP Q, FP rs, FP errormax, FP de0, FP kepernyo_high, FP kepernyo_tav, FP sugar_ki, FP gyuru_sugar_kicsi, FP gyuru_sugar_nagy, uint64_t SZELESregi, uint64_t MAGASregi, uint64_t ikezd, uint64_t jkezd, uint64_t iveg)//ekkor a SZIN egy FP* es a homersekletet reprezentalja
{

    FP* x_d = NULL;
    FP* Omega_d = NULL;
    cudaError_t err = cudaSuccess;

    auto start = std::chrono::high_resolution_clock::now();

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

    int minGridSize = 0, blockSize = 0;
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, ray_step<FP>, 0, 0);

    // Convert optimal 1D block size to 2D
    int blockDimX = sqrt(blockSize);
    int blockDimY = blockSize / blockDimX;

    dim3 threadsPerBlock(blockDimX, blockDimY);

    // Compute the required number of blocks
    int xdim = (SZELES + blockDimX - 1) / blockDimX;
    int ydim = (MAGAS + blockDimY - 1) / blockDimY;
    dim3 numBlocks(xdim, ydim);

    ray_step_T <<<numBlocks, threadsPerBlock >>> (SZIN_d, SZELES, MAGAS, x_d, Omega_d, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy, SZELESregi, MAGASregi, ikezd, jkezd, iveg);

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

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    //std::cout << "Teljes lefutasi ido\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";

    std::cout << double(duration.count()) / 1000000 << "sec\n" << (1 / (double(duration.count()) / 1000000)) << "fps\n";

    //console_kep(MAGAS, SZELES, SZIN);

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

uint64_t n_oszto(uint64_t SZELES, uint64_t MAGAS, uint64_t kepernyoSZELES, uint64_t kepernyoMAGAS, uint64_t n)
{
    uint64_t oszto = 1;
    uint64_t kis_kep=kepernyoMAGAS;
    uint64_t num = 1;
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
