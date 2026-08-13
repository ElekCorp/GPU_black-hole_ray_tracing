# GPU_black-hole_ray_tracing
Ray tracing in general space time. Currently implemented with the Kerr metric.
![blackholle07_07_2021_13_01_37](https://user-images.githubusercontent.com/93953948/206489449-9d9e8d47-fcda-4f07-a7c3-b4baf9f5e663.gif)


# Installation OpenACC version
Use the AppImage on Linux.
Also use the AppImage on Windows trough WSL2.

# Layout

| Path | Contents |
| --- | --- |
| `OpenACC_cli/` | The renderer: C++/OpenACC core, CLI, web server and movie scripts. |
| `tools/` | Code generators for the renderer. |
| `notebooks/` | Symbolic derivations and numerical validation — see `notebooks/README.md`. |
| `docs/` | Design notes and investigations. |

## Regenerating the geodesic integrator

The `christoffel()` body in `OpenACC_cli/cuda_ray.h` is generated symbolically
rather than written by hand:

```bash
pixi run christoffel
```

`tools/gen_christoffel.py` also takes `--verify` (checks both the generated form
and the block currently in `cuda_ray.h` against an untouched symbolic
reference), `--bench` (compares simplification pipelines) and `--check`.

The same geodesic in canonical (Hamiltonian) form lives in
`OpenACC_cli/hamiltonian.h`, generated the same way:

```bash
pixi run hamiltonian --write
```

It is what the implicit symplectic integrators in `OpenACC_cli/symplectic.h`
run on. `--verify` cross-checks it against the christoffel generator, so the
two derivations are held to each other rather than each to itself.

## Integrators

| | |
| --- | --- |
| `make test` | Correctness of both integrator families, including observed convergence orders. |
| `make bench` | Work-precision comparison of all of them. |

[`docs/integrator-comparison.md`](docs/integrator-comparison.md) has the
measured results and the recommendation: Dormand-Prince stays the render
default, the symplectic methods are for photon-ring work and for validation,
and the shipped controller's `1e-8` tolerance floor is where its accuracy
actually stops.




