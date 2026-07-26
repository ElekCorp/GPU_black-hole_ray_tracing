#!/usr/bin/env python
import sys
import argparse
import numpy as np
import cupy as cp
import os
import math
from pathlib import Path

def generate_and_load_kernels():
    # Only generate if the .cu files don't exist
    precisions = ['float', 'double']
    for prec in precisions:
        cu_file = f"kernel_{prec}.cu"
        if not os.path.exists(cu_file):
            import create_cu_code
            with open(cu_file, "w") as f:
                f.write(create_cu_code.create_cu_code(prec))
                
generate_and_load_kernels()

# Load modules
with open("kernel_float.cu", "r") as f:
    mod_float = cp.RawModule(code=f.read())
ker_float = mod_float.get_function("ray_step_kernel_T")

with open("kernel_double.cu", "r") as f:
    mod_double = cp.RawModule(code=f.read())
ker_double = mod_double.get_function("ray_step_kernel_T")

def datasaver_T(szin, szeles, magas, filename):
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    with open(filename, 'wb') as f:
        f.write(np.int32(magas).tobytes())
        f.write(np.int32(szeles).tobytes())
        # The C++ code writes it as a flattened array of float32s
        f.write(szin.astype(np.float32).tobytes())

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--errormax", type=float, default=0.0001)
    parser.add_argument("--de0", type=float, default=0.01)
    parser.add_argument("--rs", type=float, default=0.05)
    parser.add_argument("--delta_a", type=float, default=0.0001)
    parser.add_argument("--a", type=float, default=0.0)
    parser.add_argument("--Q", type=float, default=0.0)
    parser.add_argument("--t0", type=float, default=0.0)
    parser.add_argument("--r0", type=float, default=1.0)
    parser.add_argument("--theta0", type=float, default=1.57 + 0.06)
    parser.add_argument("--phi0", type=float, default=0.0)
    parser.add_argument("--kepernyo_high", type=float, default=0.5)
    parser.add_argument("--kepernyo_tav", type=float, default=0.75)
    parser.add_argument("--sugar_ki", type=float, default=1.01)
    parser.add_argument("--gyuru_kicsi", type=float, default=0.1)
    parser.add_argument("--gyuru_nagy", type=float, default=0.5)
    
    parser.add_argument("--kepernyoSZELES", type=int, default=10240)
    parser.add_argument("--kepernyoMAGAS", type=int, default=5120)
    parser.add_argument("--SZELES", type=int, default=640)
    parser.add_argument("--MAGAS", type=int, default=320)
    parser.add_argument("--ikezd", type=int, default=0)
    parser.add_argument("--jkezd", type=int, default=0)
    parser.add_argument("--iveg", type=int, default=0)
    
    parser.add_argument("--double", action="store_true")
    parser.add_argument("--float", action="store_true")
    
    # Parse arguments
    # Note: cli_parser.cpp handles unknown args by skipping them. 
    # argparse might complain, so we use parse_known_args
    args, _ = parser.parse_known_args()
    
    if args.iveg == 0:
        args.iveg = args.kepernyoSZELES
        
    use_double = args.double

    if use_double:
        np_dtype = np.float64
        cp_dtype = cp.float64
        ker = ker_double
    else:
        np_dtype = np.float32
        cp_dtype = cp.float32
        ker = ker_float
        
    x = np.array([args.t0, args.r0, args.theta0, args.phi0], dtype=np_dtype)
    Omega = np.array([0.0, math.asin(1.0) * 2, 0.0], dtype=np_dtype)
    
    x_gpu = cp.asarray(x)
    Omega_gpu = cp.asarray(Omega)
    szin_gpu = cp.zeros(args.SZELES * args.MAGAS, dtype=cp_dtype)
    
    block_x, block_y = 16, 16
    grid_x = (args.SZELES + block_x - 1) // block_x
    grid_y = (args.MAGAS + block_y - 1) // block_y
    
    print("Launching GPU Kernel...")
    
    ker((grid_x, grid_y, 1), (block_x, block_y, 1), (
        szin_gpu,
        np.uint64(args.SZELES),
        np.uint64(args.MAGAS),
        x_gpu,
        Omega_gpu,
        np_dtype(args.a),
        np_dtype(args.Q),
        np_dtype(args.rs),
        np_dtype(args.errormax),
        np_dtype(args.de0),
        np_dtype(args.kepernyo_high),
        np_dtype(args.kepernyo_tav),
        np_dtype(args.sugar_ki),
        np_dtype(args.gyuru_kicsi),
        np_dtype(args.gyuru_nagy),
        np.uint64(args.kepernyoSZELES),
        np.uint64(args.kepernyoMAGAS),
        np.uint64(args.ikezd),
        np.uint64(args.jkezd),
        np.uint64(args.iveg)
    ))
    
    cp.cuda.Stream.null.synchronize()
    
    szin_cpu = szin_gpu.get()
    
    out_file = "./web_images/kep_cli.dat"
    datasaver_T(szin_cpu, args.SZELES, args.MAGAS, out_file)
    print(f"Saved to {out_file}")

if __name__ == "__main__":
    main()
