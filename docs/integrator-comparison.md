# Implicit symplectic geodesic integrators, and how they compare

The renderer integrates the geodesic equation with an explicit, adaptive
Dormand–Prince 5(4) on the second-order form `x'' = -Γ(x)(v,v)`
(`OpenACC_cli/cuda_ray.h`). This note covers a second family added alongside
it — implicit, symplectic Gauss–Legendre collocation on the canonical form —
and what measuring the two against each other actually shows.

Short version: **for a single rendered ray, Dormand–Prince is still the
cheapest way to a given accuracy, by roughly 3×.** It has one hard limit that
the new methods do not — its tolerance is floored at `1e-8` inside
`dopri54_step`, which caps a 4-winding photon-ring ray at ~1e-4 of position
error no matter what is asked of it — and it drifts on long integrations, where
the symplectic methods do not drift at all.

## What was added

| File | Contents |
| --- | --- |
| `OpenACC_cli/hamiltonian.h` | The geodesic equation in canonical form, `H = ½ g^{μν}(x) p_μ p_ν`. Generated. |
| `tools/gen_hamiltonian.py` | The generator for it (`pixi run hamiltonian`), sharing its metric with `gen_christoffel.py`. |
| `OpenACC_cli/symplectic.h` | Gauss–Legendre orders 2/4/6, the Sundman time transformation, and an `(x, v)` wrapper. |
| `OpenACC_cli/tests/test_symplectic.cpp` | Correctness: `make test`. |
| `OpenACC_cli/tests/bench_integrators.cpp` | The measurements below: `make bench`. |

### Why the change of variables

Symplecticity is a property of a map on a *canonical* pair. In the renderer's
`(x, v)` state the symplectic form pulled back to the tangent bundle is
`g_{μν}(x) dx^μ ∧ dv^ν`, which is not canonical, so applying a "symplectic"
method there yields a method that is merely symmetric. Moving to `(x, p)` with
`p_μ = g_{μν} v^ν` fixes that, and pays for itself three more times over:

* **`H` is the null constraint.** `g^{μν} p_μ p_ν = 0` is what makes a
  worldline a photon, and it is exactly the Hamiltonian, so `|H|` is a free,
  physically meaningful error measure.
* **Energy and axial angular momentum become exact.** `t` and `φ` are cyclic,
  so `dp₀ = dp₃ = 0` identically in the generated code. Both are conserved to
  the last bit — by *any* integrator run in these variables, not just the
  symplectic ones. The test asserts bit-exact, not "within tolerance".
* **The right-hand side is cheaper.** Contracting `g^{μν}` with `p` twice is
  less work than contracting `Γ^i_{jk}` with `v` twice: 128 flops against 216
  in the general case, and a measured 28 ns per evaluation — including the
  implicit method's own iteration overhead — against 37 ns.

`gen_hamiltonian.py --verify` checks the emitted code against the raw symbolic
form, and independently against `gen_christoffel.py` — that
`d/dλ(g_{μν}v^ν)` computed from `x'' = -Γ(v,v)` equals `-½ ∂_μ g^{αβ} p_α p_β`.
Both agree to 5e-15 over 3000 random states per specialization.

### The methods

`gauss_legendre_step<FP, S>` is the s-stage collocation method: implicit, order
2s, symplectic and symmetric. `S = 1` is the implicit midpoint rule. The stage
equations are solved by fixed-point iteration — geodesics are non-stiff, so it
contracts at a rate ~`h·|∂f/∂z|` with no Jacobian and no linear solve, which
also keeps it usable inside a per-pixel GPU kernel. Measured cost at a 1e-14
stage tolerance is 3.5–5.5 iterations per step, so a `gl4` step is 9–11
right-hand side evaluations against Dormand–Prince's 7 per step at 5th order.

The iteration runs to the stage tolerance rather than stopping at the local
truncation error, because an inexact stage solve perturbs the map by the
residual and reintroduces exactly the drift the method exists to avoid.

**Constant step size is a condition, not a preference.** The backward-error
argument that makes a symplectic method conserve `H` says it conserves the
Hamiltonian of a nearby modified system — and that system depends on `h`, so
changing `h` changes what is being conserved. This is a bad fit for a ray that
starts 50 M out and grazes the photon sphere. The way out is to keep the step
constant and change the *independent variable*: integrating

    K(x, p) = s(x) [ H(x, p) − H₀ ],   s(x) = (r / r_ref)^κ

at fixed step in a fictitious time `τ` is still canonical, hence still
symplectic, and on the level set `K = 0` its orbits are the orbits of `H` with
`dλ = s(x) dτ`. `κ = 1` turned out to be the right exponent for light; `κ = 2`
over-concentrates work in the strong field and loses.

## Measurements

`make bench` (about 5 s). Kerr `a = 0.9`, `rs = 2` (M = 1), double precision,
Apple M-series, clang -O2. Errors are against a converged reference — 6th-order
symplectic at a tiny step, cross-checked against explicit RK4 at a tiny step,
the two agreeing to the figure quoted per scenario. `rk4` on the canonical form
is the control: same variables as the symplectic methods, same explicitness as
Dormand–Prince, no symplectic structure, so the columns separate "what the
change of variables bought" from "what symplecticity bought".

### 1. A strongly lensed ray — one pixel of a render

Camera at r = 50, grazes to r = 3.9, escapes; λ = 0..140. Reference good to 2e-13.

| method | control | rhs evals | ms | position error | max \|ΔH\| |
| --- | --- | ---: | ---: | ---: | ---: |
| dopri54 | 1e-6 | 413 | 0.017 | 3.2e-06 | 6.0e-06 |
| dopri54 | 1e-8 *(floor)* | 749 | 0.028 | **3.8e-08** | 4.4e-08 |
| rk4 | h=0.2 | 2800 | 0.109 | 3.5e-08 | 3.9e-08 |
| gl4 | h=0.2 | 7368 | 0.210 | 8.7e-09 | 1.9e-08 |
| gl4-tt (κ=1) | dτ=0.2 | 2048 | 0.065 | **1.9e-08** | 4.6e-08 |
| gl4-tt (κ=1) | dτ=0.025 | 10682 | 0.351 | 4.7e-12 | 1.1e-11 |
| gl6 | h=0.2 | 11037 | 0.307 | 1.0e-12 | 1.5e-12 |

At the accuracy a render needs, Dormand–Prince gets there in 749 evaluations
and the best symplectic configuration needs 2048. **Dormand–Prince wins the
work–precision comparison on this ray by about 2.7×**, and nothing in the
symplectic column changes that: a photon makes one pass through the strong
field, picks up its error there, and leaves. There is no long-time drift to
prevent, which is the entire advantage being paid for.

### 2. A photon-ring ray — the case the zoom movies are made of

Tuned to the critical impact parameter to within 1e-4: four windings around the
prograde photon sphere at r = 1.566, then escapes. λ = 0..115. The winding makes
the ray exponentially sensitive to integration error. Reference good to 2.5e-12.

| method | control | rhs evals | ms | position error | max \|ΔH\| |
| --- | --- | ---: | ---: | ---: | ---: |
| dopri54 | 1e-6 | 546 | 0.021 | 9.3e-03 | 2.5e-05 |
| dopri54 | 1e-8 *(floor)* | 1141 | 0.039 | **1.1e-04** | 1.2e-07 |
| rk4 | h=0.2 | 2300 | 0.083 | 1.6e-04 | 1.2e-06 |
| rk4 | h=0.05 | 9200 | 0.336 | 5.9e-07 | 4.6e-09 |
| gl4-tt (κ=1) | dτ=0.2 | 3686 | 0.116 | 4.0e-06 | 5.7e-08 |
| gl4-tt (κ=1) | dτ=0.025 | 19042 | 0.602 | 9.7e-10 | 1.4e-11 |
| gl6 | h=0.2 | 10299 | 0.262 | **8.1e-11** | 3.6e-11 |

This is the one qualitative difference, and it is not about symplecticity — it
is about the floor. `dopri54_step` clamps its own tolerance:

```c
FP const tolerance = fmax(hole.errormax, FP(1e-8));
```

so `errormax` below 1e-8 does nothing, and this ray cannot be integrated better
than ~1e-4 by the shipped path at any cost. Four windings amplify a per-step
error of 1e-8 into 1e-4 at the exit. Every other method here has no such
ceiling: `gl6` reaches 8e-11 — six orders better — for 9× the evaluations, and
`gl4-tt` reaches 4e-6 for 3×. For photon-ring zooms, and for validating what
the renderer produces there, that is the difference between resolvable and not.

### 3. A bound timelike orbit — where symplecticity actually shows

Not a ray, but the same integrators on the same equations: an eccentric orbit
from r = 20 down to r = 8.9, 66 revolutions, λ = 0..20000. Position error here
bottoms out at 3e-10, the reference's own roundoff floor over 10⁵ steps.

| method | control | rhs evals | position error | max \|ΔH\| |
| --- | --- | ---: | ---: | ---: |
| dopri54 | 1e-8 *(floor)* | 26873 | 3.5e-04 | 4.5e-08 |
| rk4 | h=0.5 | 160000 | 5.9e-07 | 2.1e-10 |
| gl4 | h=1 | 264622 | 4.1e-07 | 1.2e-10 |
| gl6 | h=2 | 223908 | **1.1e-09** | 3.4e-13 |
| gl4-tt (κ=2) | dτ=0.02 | 553064 | 1.9e-10 | 8.2e-14 |

Growth of the `|H|` error, as the running maximum at each quarter of the run —
this is the measurement the whole exercise is about:

| method | 1/4 | 1/2 | 3/4 | 4/4 | 2nd half |
| --- | ---: | ---: | ---: | ---: | --- |
| dopri54 @1e-8 | 1.17e-08 | 2.30e-08 | 3.34e-08 | 4.50e-08 | 1.95× **growing** |
| rk4 h=1 | 1.64e-09 | 3.39e-09 | 5.03e-09 | 6.67e-09 | 1.97× **growing** |
| gl2 h=1 | 2.41e-06 | 2.41e-06 | 2.41e-06 | 2.41e-06 | 1.00× bounded |
| gl4 h=1 | 1.17e-10 | 1.17e-10 | 1.17e-10 | 1.17e-10 | 1.00× bounded |
| gl4-tt h=0.2 | 1.45e-11 | 1.45e-11 | 1.46e-11 | 1.46e-11 | 1.00× bounded |

Both explicit methods accumulate error in `H` linearly in the elapsed
parameter — the ratio is 1.95 and 1.97 where linear growth predicts 2.00. Every
symplectic method is flat to three digits, including the *second-order* one:
`gl2` holds `|ΔH|` at 2.4e-06 forever, while `rk4`, four orders more accurate
per step, is on a trajectory that will pass it given enough revolutions. That
is the textbook result, reproduced on this codebase's own equations. (`gl6` at
1e-15 is the exception, and only because it is at the roundoff floor, where
nothing is conserved exactly.)

## What to use

* **Keep Dormand–Prince as the renderer default.** It is the cheapest per digit
  on ordinary rays and it is already there. Nothing in these measurements
  argues for switching the render path.
* **Reach for `gl6` or `gl4-tt` (κ=1) for photon-ring work**, and as the
  independent check on what the renderer produces anywhere. They are the only
  options that get below the 1e-8 tolerance floor.
* **Consider raising that floor.** `fmax(hole.errormax, 1e-8)` is reasonable for
  a float32 render and is the binding constraint in double precision. It was
  not changed here — that is a render-path decision with its own testing — but
  it is where the shipped integrator's accuracy actually stops.
* `gauss_legendre_step_v` takes the renderer's `(x, v)` state directly, so
  swapping the integrator for an experiment is a one-line change; the
  conversion at each end is exact algebra.

## Caveats

* All measurements are CPU, double precision, single ray. The renderer's
  regime is GPU float32 across a million rays, where register pressure (the
  s-stage method holds `S × 8` doubles of stage state) and warp divergence in
  the iteration count matter and are not measured here.
* Dormand–Prince's evaluation count is charged at 7 per *accepted* step;
  rejected attempts are not counted, so its column is a lower bound — the
  comparison is being generous to the incumbent.
* The fixed-point iteration is not a contraction for arbitrarily large steps.
  `gauss_legendre_step_safe` retries a failed step as 2, 4, … substeps, which
  breaks the constant-step condition; it is a safety net for rays about to be
  captured, not a step-size controller.
