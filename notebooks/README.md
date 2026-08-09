# Notebooks

Exploratory and derivation work. Nothing here is on the render path — the
shipping code is `OpenACC_cli/`. Run these from the root pixi environment
(`pixi run jupyter lab`).

## `derivation/`

Symbolic derivations of the Kerr–Newman geodesic equation, i.e. where the
`christoffel()` body in `OpenACC_cli/cuda_ray.h` came from.

| Notebook | What it does |
| --- | --- |
| `kerr_christoffel.ipynb` | sympy `diffgeom` → Christoffel symbols → contract with `v` → CSE → C. The original generator for the block in `cuda_ray.h`. |
| `christoffel_jax_verif_sympy.ipynb` | Near-duplicate of the above, kept because it carries the C-codegen cells (`pown` user-function mapping). |
| `kerr_metric_code_generator_sage_math.ipynb` | Earlier SageMath attempt at the same derivation, using an orthonormal tetrad. Superseded. |
| `jax_christoffel.ipynb` | Cross-check: the sympy-derived Christoffel symbols against JAX autodiff of the metric, over 100k random states. |

**These are superseded by [`tools/gen_christoffel.py`](../tools/gen_christoffel.py)**,
which reproduces the same result as a committed, verifiable script rather than
notebook output pasted into a header. Prefer it for regenerating the C.

One trap is preserved in `kerr_christoffel.ipynb` and worth knowing about: the
cell that calls `.simplify()` on the contracted RHS takes ~10 minutes and makes
the output *worse*, introducing `cos(4θ)` / `cos(6θ)` multiple-angle terms. That
is the origin of the large commented-out block still sitting in `cuda_ray.h`
below the live one. Use `factor_terms` instead.

## `validation/`

Numerical cross-checks of the integrator against independent ODE solvers.

| Notebook | What it does |
| --- | --- |
| `jax_blackhole.ipynb` | JAX reimplementation of the ray tracer, used as a reference for the C++ integrator. |
| `jax_blackhole-diffrax_sanity_check.ipynb` | Same geodesics through `diffrax`'s adaptive solvers, to check the hand-rolled Dormand–Prince 5(4) step. |
| `diffrax_blackhole-Copy1.ipynb` | Working copy of the above with additional step-size experiments. |

## `jax_ein.py`

Standalone helper: metric / Christoffel / geodesic RHS in JAX via `einsum`,
imported by the validation notebooks' earlier revisions. Kept as reference for
the autodiff formulation.
