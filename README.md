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




