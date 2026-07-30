#include <iostream>
//#include <fstream>

#include <chrono>
#include <filesystem>
#include <system_error>

//#include <vector>

#include <math.h>
//#include <vector_types.h>

#include "black_hole.h"
#include "cuda_ray.h"

#include "szinsaver.h"

#include "cli_parser.h"

//#include "debugmalloc.h"

template <class FP>
int8_t* makeframe(uint64_t const SZELES, uint64_t const MAGAS, FP const* const x, FP const* const Omega, FP a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg);

template <class FP>
FP* makeframe_T(uint64_t const SZELES, uint64_t const MAGAS, FP const* const x, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg);

void device_info(void);

uint64_t n_oszto(uint64_t const SZELES, uint64_t const MAGAS, uint64_t const kepernyoSZELES, uint64_t const kepernyoMAGAS, uint64_t const n);


int main(int argc, char* argv[])
{
    device_info();

    // Batch jobs commonly start in a fresh working directory.  The renderer
    // writes its hit buffer here before the Python compositor creates the PNG.
    std::error_code output_directory_error;
    std::filesystem::create_directories("web_images", output_directory_error);
    if (output_directory_error)
    {
        std::cerr << "Cannot create output directory ./web_images: "
                  << output_directory_error.message() << std::endl;
        return 1;
    }

    Params p;
    parse_args(argc,argv,p);
    //auto start = std::chrono::high_resolution_clock::now();

    if (p.kepernyoSZELES * p.MAGAS != p.kepernyoMAGAS * p.SZELES)
    {
        std::cout << "renderelendo kep es a kepernyo aranya nem azonos\n";
        return 1;
    }

    float const x[D] = { float(p.t_0),float(p.r_0),float(p.theta_0),float(p.phi_0) };


    double const pi_cucc = (asin(1) * 2);
    float const Omega[D - 1] = { 0,float(pi_cucc),0 };




    uint64_t SZELESregi = p.kepernyoSZELES;
    uint64_t MAGASregi = p.kepernyoMAGAS;

    uint64_t ikezd = p.ikezd;
    uint64_t jkezd = p.jkezd;
    uint64_t iveg = p.iveg;


    double const x_d[D] = { p.t_0, p.r_0, p.theta_0, p.phi_0 };
    double const Omega_d[D - 1] = { double(Omega[0]),double(Omega[1]),double(Omega[2]) };

    std::cout<<p.iveg<<"/n"<<p.ikezd<<"/n"<<p.jkezd<<"/n"<<SZELESregi<<std::endl;

if(p.prec==Precession::Double)
{
    double* SZIN = NULL;
    SZIN = makeframe_T<double>(p.SZELES, p.MAGAS, x_d, Omega_d, p.a, p.Q, p.rs, p.errormax, p.de0, p.kepernyo_high, p.kepernyo_tav, p.sugar_ki, p.gyuru_sugar_kicsi, p.gyuru_sugar_nagy, SZELESregi, MAGASregi, ikezd, jkezd, iveg);
    std::string kep_double_string="./web_images/kep_cli.dat";
    datasaver_T<double>(SZIN, p.SZELES, p.MAGAS, kep_double_string, 3);

    free(SZIN);
}
else
{
    float* SZIN_f=NULL;
    SZIN_f = makeframe_T<float>(p.SZELES, p.MAGAS, x, Omega,float(p.a),float(p.Q), float(p.rs), float(p.errormax),float(p.de0), float(p.kepernyo_high), float(p.kepernyo_tav), float(p.sugar_ki), float(p.gyuru_sugar_kicsi), float(p.gyuru_sugar_nagy), SZELESregi, MAGASregi, ikezd, jkezd, iveg);
    std::string kep_string="./web_images/kep_cli.dat";
    datasaver_T<float>(SZIN_f, p.SZELES, p.MAGAS, kep_string, 3);

    free(SZIN_f);
}


    return 0;
}

template <class FP>
int8_t* makeframe(uint64_t const SZELES, uint64_t const MAGAS, FP const* const x, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg)
{

    FP const* x_d = x;
    FP const* Omega_d = Omega;

    //auto start = std::chrono::high_resolution_clock::now();


    int8_t* SZIN_d = NULL;
    SZIN_d = (int8_t*)malloc(SZELES * MAGAS * sizeof(int8_t));


    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    ray_step(SZIN_d, SZELES, MAGAS, x_d, Omega_d, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy, SZELESregi, MAGASregi, ikezd, jkezd, iveg);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000.0 << "[s]" << std::endl;

    int8_t* SZIN = SZIN_d;



    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    //std::cout << "Teljes lefutasi ido\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";

    //std::cout << double(duration.count()) / 1000000 << "sec\n"<<(1/(double(duration.count()) / 1000000))<<"fps\n";

    //console_kep(MAGAS, SZELES, SZIN);

    return SZIN;

}


template <class FP>
FP* makeframe_T(uint64_t const SZELES, uint64_t const MAGAS, FP const* const x, FP const* const Omega, FP const a, FP const Q, FP const rs, FP const errormax, FP const de0, FP const kepernyo_high, FP const kepernyo_tav, FP const sugar_ki, FP const gyuru_sugar_kicsi, FP const gyuru_sugar_nagy, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg)//ekkor a SZIN egy FP* es a homersekletet reprezentalja
{

    FP const* x_d = x;
    FP const* Omega_d = Omega;

    //auto start = std::chrono::high_resolution_clock::now();


    FP* SZIN_d = NULL;
    // Each pixel carries disk radius, ray-traced frequency shift, and the
    // local disk azimuth used to map emissivity structure through lensing.
    SZIN_d = (FP*)malloc(3 * SZELES * MAGAS * sizeof(FP));

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    ray_step_T(SZIN_d, SZELES, MAGAS, x_d, Omega_d, a, Q, rs, errormax, de0, kepernyo_high, kepernyo_tav, sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy, SZELESregi, MAGASregi, ikezd, jkezd, iveg);

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count()/1000000.0 << "[s]" << std::endl;


    FP* SZIN = SZIN_d;

    //auto end = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    //std::cout << "Teljes lefutasi ido\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";
    //std::cout << "\n\n\n\n\n\n\n";

    //std::cout << double(duration.count()) / 1000000 << "sec\n"<<(1/(double(duration.count()) / 1000000))<<"fps\n";

    //console_kep(MAGAS, SZELES, SZIN);

    return SZIN;

}


void device_info(void)
{
    /*int num_gpus = 0;
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

    }*/
}


int n_oszto(uint64_t const SZELES, uint64_t const MAGAS, uint64_t const kepernyoSZELES, uint64_t const kepernyoMAGAS, int const n)
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
