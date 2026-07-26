import re

def create_cu_code(precision='float'):
    with open("black_hole.h", "r") as f:
        bh_code = f.read()
    with open("cuda_ray.h", "r") as f:
        cr_code = f.read()

    # Remove includes
    bh_code = '\n'.join([line for line in bh_code.split('\n') if not line.startswith('#include')])
    cr_code = '\n'.join([line for line in cr_code.split('\n') if not line.startswith('#include')])

    # Typedefs and Macros
    preamble = f"""
typedef {precision} FP;
typedef unsigned long long uint64_t;
typedef signed char int8_t;
"""

    # Process black_hole.h
    bh_code = bh_code.replace('template <class FP>', '')
    bh_code = bh_code.replace('kerr_black_hole<FP>', 'kerr_black_hole')
    # Remove iro completely to avoid printf and PRIu64
    import re
    bh_code = re.sub(r'void iro\(void\)[\s\S]*?}', '', bh_code)
    bh_code = bh_code.replace('kerr_black_hole(', '__device__ kerr_black_hole(')

    # Process cuda_ray.h - extract everything except ray_step and ray_step_T
    cr_code = cr_code.replace('template <class FP>', '')
    cr_code = cr_code.replace('inline ', '__device__ inline ')
    cr_code = cr_code.replace('kerr_black_hole<FP>', 'kerr_black_hole')
    cr_code = cr_code.replace('pow(hole.errormax / err, 0.2)', 'pow(hole.errormax / err, FP(0.2))')
    
    # Slice out ray_step and ray_step_T
    # We must match the definition, not the forward declaration. The definition has "//kernel" at the end.
    idx_start = cr_code.find('__device__ inline void ray_step(int8_t* const szin, uint64_t const SZELES, uint64_t const MAGAS, FP const* const xd')
    idx_end = cr_code.find('//ha a gömbön')
    
    if idx_start != -1 and idx_end != -1:
        cr_code = cr_code[:idx_start] + cr_code[idx_end:]

    # Now add our custom kernels
    kernels = """
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
"""

    return preamble + bh_code + cr_code + kernels

if __name__ == "__main__":
    import sys
    prec = sys.argv[1] if len(sys.argv) > 1 else 'float'
    with open(f"kernel_{prec}.cu", "w") as f:
        f.write(create_cu_code(prec))
