# christoffel / integrator / precision benchmark

Compares four `cuda_ray.h` build variants (current `christoffel_general` vs.
the pre-CSE 84-temp `christoffel()` from commit `9bc60f5`, and `dopri54_step`
vs. the older fixed-step `RK6`), each in fp32/fp64, in a wide (full-frame)
render and a zoomed (1e6x, centred on the shadow's photon-ring edge) render.

The two BENCH_* macros live in `cuda_ray.h` itself (see
`christoffel_old87()`, and the `#if defined(BENCH_OLD_CHRISTOFFEL)` /
`#if defined(BENCH_RK6)` switches in `christoffel()` and `step()`) so the
harness can select a variant with a `-D` flag instead of patching source.
Neither macro is defined by a normal `make` — plain builds are unaffected.

## Run it

Needs `nvc++` and a GPU (this project builds via `distrobox enter nvcc --`
on the machine this was developed on):

```
distrobox enter nvcc -- bash benchmarks/run_bench.sh
python3 benchmarks/summarize.py
```

`run_bench.sh` rebuilds `./main` four times (`make clean && make
COMMON_FLAGS="-std=c++17 -Wall <defs>"` per variant), runs it
`{float,double} x {wide,zoomed} x 7` times per build, and writes
`benchmarks/results/results.csv` (git-ignored — it's regenerated per run).
It finishes by rebuilding the plain, unbenchmarked binary so `./main` and
`libblackhole.so` are left in their normal state.

## Changing the scene

`EDGE_CY` in `run_bench.sh` is a bisected constant, not a formula — it was
found by `find_shadow_edge.py` for the specific scene run_bench.sh uses
(a=0.02, Q=0, rs=0.05, r0=1.0, theta0=1.63). If you change any of those,
rebuild the plain binary and rerun `find_shadow_edge.py` to get a new
`EDGE_CY` before rerunning the benchmark, or "zoomed" will be zooming into
whatever used to be the edge.
