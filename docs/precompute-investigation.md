# Precomputing ray paths for a portable viewer — investigation

Goal: rent interruptible A100s on vast.ai, precompute geodesic data at a
resolution the desktop app could never reach live, ship that data inside an
easily-installable app that only *shades* it.

This document is step one: what is actually precomputable, what it costs, and
what has to change in the renderer before any GPU-hours are worth spending.

---

## 1. The precomputable object already exists

The codebase is already split exactly along the right seam:

| Layer | File | Cost | Depends on |
| --- | --- | --- | --- |
| Geodesic integration | `cuda_ray.h` `ray_step_T` | GPU, expensive | metric + camera + disk geometry |
| Shading / film finish | `cli_imagemaker.py` | numpy, cheap | artistic + physical choices |

`ray_step_T` writes `RAY_CHANNELS = 6` floats per pixel and nothing downstream
re-derives a trajectory. That 6-channel buffer **is** the precompute product.
Everything in `cli_imagemaker.py` — the lensed sky plate, the thermal disk
model, MRI turbulence, retarded-time animation, bloom, ACES tonemap — runs on
CPU at interactive rates today and should stay in the app. Precomputing it
would only freeze artistic choices the user should keep.

So the question is not *whether* to precompute the hit buffer, but **over which
parameters**, and **in what representation**.

---

## 2. The parameter space, and how much of it collapses

Full input surface of one trace (`ffi_bridge.cpp` + `cli_parser.h`):

```
metric:     rs, a, Q
camera:     r0, theta0, phi0, t0, Omega (axis-angle, 3 dof)
screen:     kepernyo_high, kepernyo_tav, SZELES, MAGAS
zoom:       roi_cx, roi_cy, roi_w
geometry:   sugar_ki, gyuru_sugar_kicsi, gyuru_sugar_nagy
accuracy:   errormax, de0
```

That is ~19 numbers. Most of them are not real dimensions.

### Exactly free (symmetry, not approximation)

**`phi0` and `t0`.** `christoffel()` reads only `x[1]` and `x[2]` (r and θ) —
never `x[0]` or `x[3]`. `ijk_to_vec_zoom` builds the tetrad from `r_0`,
`theta_0`, `rs`, `a`, `Q` and never touches `phi_0`. Kerr–Newman is stationary
and axisymmetric, so shifting `phi0` by Δφ shifts the whole solution by Δφ.
Recover it in the app: add Δφ to channel 2 (disk azimuth) and channel 5
(`sky_phi`); channels 0, 1, 3, 4 are unchanged. **Orbiting the hole in azimuth
costs zero storage.**

**`rs`.** The metric has an exact scaling symmetry `(r, a, Q, rs, t) → λ·(…)`
that leaves null geodesics invariant as curves. `rs = 0.05` is hardcoded in
`web_server.py` anyway. Precompute in units of M and never store `rs`.
(Current units: `rs = 2M = 0.05` ⟹ M = 0.025, so `r0 = 1` is **40 M**, disk
`0.1–0.5` is **4 M – 20 M**, horizon `0.05` is 2 M, photon sphere 3 M, ISCO 6 M.)

**`Omega` (camera orientation, 3 dof).** `ijk_to_vec_zoom` applies the rotation
to the *Minkowski* screen direction **before** the tetrad transform — it is a
pure rotation in the camera's local frame. If the precompute stores the result
as a map over **local ray direction on the full sphere** rather than over pixel
indices, then pan/tilt/roll is a rotation of the lookup and is free.
17× more rays than one 67°×37° frustum, but it kills 3 dimensions and makes the
app's camera controls instant.

**Screen params and `roi_*` (6 dof).** `kepernyo_high`, `kepernyo_tav`,
`SZELES`, `MAGAS`, and the zoom window all just select which directions out of
the same continuous ray-direction field get sampled. Against a direction map
they are resampling, not retracing. Zoom is free **up to the stored angular
resolution** — see §5, this is the hard limit of the whole approach.

**`errormax`, `de0`.** Not scene parameters; converge them once and freeze.
A100 runs fp64 at 1:2 vs fp32 (9.7 vs 19.5 TFLOP/s) — unlike the consumer 1:32
or 1:64 that `ffi_bridge.cpp`'s comment assumes. **Precompute everything in
f64.**

`OpenACC_cli/fp64_showdown.py` now measures what that buys, zooming onto the
shadow edge and scoring the f32 trace against the f64 one frame by frame. Two
f32 failures, in this order:

| | onset | what it looks like |
| --- | --- | --- |
| chaotic amplification — rays still distinct and finite, individually wrong | agreement < 90% by **~1.3e3×** | salt-and-pepper noise |
| rays going non-finite — integration fails, `sky_direction` substitutes a sentinel | **~1e5×**; 13.8% of the frame at 1e6× vs f64's 7.1% | flat patches |

**Neither precision suffers genuine ray collapse** — adjacent pixels rounding to
identical launch directions — anywhere up to 1e6×: over 24 frames, f32 collapsed
zero neighbouring pairs and f64 collapsed three out of 12,045 in one frame,
because `ijk_to_vec_mink_zoom` computes launch directions in double even for the
f32 trace. That matters for the precompute design: the ceiling on a stored
direction map is set by the map's own angular resolution, not by any ULP limit
in how the rays were launched.

Note the first failure arrives two decades before anything runs out of bits,
because the photon sphere is a chaotic scattering region that amplifies f32
rounding exponentially — so `FLOAT_ZOOM_LIMIT = 512` in `web_server.py` is a
reasonable place to switch, but it guards against amplified error, not lost
resolution.

### Halved

**`theta0`.** Kerr–Newman is symmetric under θ → π−θ. A camera below the
equator is the mirror of one above it (mirror the direction map, flip
`sky_theta → π − sky_theta`, leave azimuth and radius alone). Store θ0 ∈ (0, π/2]
only.

### Collapsed further when a = 0

Schwarzschild and Reissner–Nordström are **spherically symmetric**, so each ray
stays in the plane containing the camera, the centre, and the initial direction.
The whole map reduces from 2 direction-dimensions to **1** (the angle from
radially-inward), and `r0` folds in too because the trajectory is fixed by the
impact parameter alone. See tier B in §5 — this is where the size problem
actually gets solved.

### Genuinely irreducible

```
a, Q            metric shape          (2)
r0, theta0      camera on the sphere  (2, θ0 halved)
ray direction   local sky             (2)
```

Six dimensions, of which two are the "pixels".

---

## 3. Change the payload *before* burning GPU-hours

Rendering the current 6-channel buffer at high resolution would bake in
restrictions that make the app much less interesting. Four changes, in
descending order of value:

### 3a. Record every equatorial crossing, not "the disk hit"

`ffi_bridge.cpp` hardcodes `gyuru_sugar_kicsi = 0.1`, `gyuru_sugar_nagy = 0.5`,
and `ray_step_T` **terminates** the ray at the first crossing inside that
annulus. Precompute that and the disk's radii are frozen into the dataset
forever.

Instead: don't terminate on the disk. Record the first N equatorial-plane
crossings as `(r, φ, t)` triples. The app then decides which crossing an
annulus actually intercepts. That buys, at no extra trace cost:

- user-adjustable inner/outer disk radius,
- optically thin / semi-transparent disks (accumulate along all crossings),
- explicit control over image order n = 0, 1, 2, … — the photon-ring structure
  the project already renders GIFs of.

N = 4 covers everything visible; successive orders are separated by e^{2π} ≈ 535
in width, so n ≥ 4 is far below any storable resolution.

### 3b. Store conserved quantities instead of the redshift

Channel 1 (`disk_redshift`) is **redundant**. Look at what it computes:

```cpp
FP const p_t   = gtt * v[0] + gtphi * v[3];
FP const p_phi = gtphi * v[0] + gphiphi * v[3];
```

`∂_t` and `∂_φ` are Killing vectors, so `p_t` and `p_φ` are conserved along the
whole geodesic — identical at every point. Store them **once per ray** (really
just `b = −p_φ/p_t`, since `E` normalises out) and the app computes g in closed
form at whatever radius, for whatever emitter orbit model it likes. A retuned
disk model no longer means a re-trace.

### 3c. Push the escape sphere out

`sugar_ki = 1.01` with `r0 = 1.0`: the sky is **0.4 M beyond the camera**, at
40.4 M. Escape directions are read there, where BL coordinates are assumed
near-orthonormal (`sky_direction`'s own comment). For a live renderer that is a
sensible cost trade; for a one-off precompute it is free to push `sugar_ki` to
10³–10⁴ M and get sky directions that are genuinely asymptotic. Costs more
steps per ray (§4), buys correctness that can never be fixed later.

### 3d. Expose the frustum

`ffi_bridge.cpp` hardcodes `kepernyo_high = 0.5`, `kepernyo_tav = 0.75` —
a 67.4° × 36.9° frustum. A full-sphere direction map needs 6 cube faces at 90°,
so these must become arguments. Small change, hard prerequisite for §2's
"orientation is free".

*(Related: `gyuru_sugar_kicsi = 0.1` = 4 M is inside the a = 0 ISCO of 6 M. Fine
as a geometric annulus, but the thin-disk emission model in `cli_imagemaker.py`
assumes a no-torque inner boundary — worth deciding deliberately.)*

---

## 4. Compute is cheap. Distribution size is the whole problem.

### Compute

Rough per-ray cost: `de0 = 0.01` caps the step, so reaching a 1000 M escape
sphere from 40 M takes ≳2500 steps; each DOPRI5 step is 7 `christoffel()`
evaluations at a few hundred flops plus several transcendentals. Call it
~10 Mflop-equivalents per typical ray, and assume 20 % of A100 fp64 peak given
the branch divergence:

**≈ 1–3 × 10⁵ rays/s/GPU** (near-critical rays that wind many times cost
10–100× more, but they are a small fraction of solid angle).

A 144-state grid × 25 M rays/state ≈ 3.7 × 10⁹ rays ⟹ **~5 GPU-hours**.
At interruptible A100 pricing that is a few dollars. **Compute is not the
constraint.**

### Storage

Full-sphere cubemap at 0.044°/texel (matches a 1280-wide render across the
native 67° frustum): 6 × 2048² = 25.2 M texels.

| Payload | Bytes/texel | Per camera state |
| --- | --- | --- |
| 6 channels, f32 (today's format) | 24 | 605 MB |
| 6 channels, f16 | 12 | 302 MB |
| §3 payload (tag + 4 crossings + b), f16 | ~20 | 504 MB |
| …lossless-compressed (smooth almost everywhere) | ~2–4 | **50–100 MB** |

Against a 144-state grid: **7–14 GB**. That is not a portable app.

**This is the finding that should drive the design.** The instinct — "rent
GPUs, brute-force a grid of viewpoints" — is affordable to compute and
impossible to ship. Effort belongs in *representation*, not in GPU-hours.

---

## 5. Three tiers

### Tier A — brute-force direction maps (works now, doesn't scale)

Trace the §3 payload over a full-sphere cubemap for a handful of hand-picked
`(a, Q, r0, θ0)` states. Free camera orientation, free φ0, free FOV, free zoom
up to 0.044°/texel (≈ **1200× magnification** before texels are visible —
comparable to the current f32 live path, well short of the f64 path's 2^32).

- **Cost:** ~2 GPU-minutes and ~50–100 MB *per state*.
- **Ship:** 4–8 states, ~500 MB. A "gallery" app, not a free-flight one.
- **Value:** real, immediate, and it is the reference dataset every other tier
  gets validated against.

### Tier B — impact-parameter tables for a = 0 (the actual answer for the non-spinning case)

Spherical symmetry (§2) collapses the map to: for each impact parameter b, the
swept angle ψ(r) and coordinate time t(r) along the trajectory. The app does
spherical trig to find where the ray plane crosses the disk plane, then two
table lookups.

- **Size:** 8192 b-samples × 512 r-samples × 2 values × f32 ≈ **33 MB, total** —
  covering *every* camera radius, *every* inclination, *every* orientation,
  *every* zoom, *every* disk geometry, and Q as one extra small axis.
- **Zoom depth:** unbounded in practice, if b is sampled in `log(b − b_crit)`
  around the critical impact parameter b_crit = 3√3 M. This is what makes the
  photon-ring zoom GIFs the project already produces work at arbitrary depth
  instead of hitting a texel floor.
- **Honest caveat:** Schwarzschild null geodesics are *exactly* solvable in
  Jacobi elliptic functions. This table can be generated on a laptop in
  minutes, or replaced entirely by ~50 lines of Carlson elliptic integrals in
  the app. **Tier B needs no vast.ai at all.**

### Tier C — Carter-constant tables for a ≠ 0 (the one worth renting GPUs for)

Photons are neutral, so null geodesics in Kerr–**Newman** are separable: the
Carter constant exists, and r-motion and θ-motion decouple. Each ray is labelled
by two constants `(b = L_z/E, q = Q_C/E²)` plus turning-point signs. Precompute
per `(a, Q)`: radial and polar Mino-time integrals over the `(b, q)` plane.

- **Size:** a few tens of MB per `(a, Q)`, ~1 GB for a decent spin grid.
- **Camera position becomes a lookup offset**, not a stored dimension — the
  same collapse tier B gets from r0, now applied to (r0, θ0) both.
- **Zoom depth:** unbounded, same log-spacing trick near the critical curve.
- **Cost:** this is real physics/engineering work (Gralla & Lupsasca 2019 give
  the closed forms; a numerical table sidesteps the case analysis but not the
  turning-point bookkeeping).

### Recommendation

Do **A** first — it is a few GPU-minutes, it forces §3's payload redesign, and
it produces the ground truth everything else is checked against. Then **C**,
because it is the only one that both scales and needs the GPUs. **B** is a
weekend of maths that makes the a = 0 case perfect and ~33 MB, and it validates
the whole pipeline against a case with an exact answer.

A caveat worth stating plainly: the README frames this as *"ray tracing in
general space time"*. Tiers B and C trade that generality for compression —
they exploit Kerr–Newman's specific symmetries and would have to be rebuilt for
any other metric. Tier A does not. If metric-agnosticism is the point of the
project, tier A is the only honest precompute, and the size problem in §4 is
inherent rather than solvable.

---

## 6. vast.ai interruptible mechanics

Preemption is the design constraint, and the existing cache layer is already
close to what is needed.

- **Work unit = one cubemap face tile**, minutes not hours. `render_frame_roi_f64`
  already renders arbitrary sub-windows; note it enforces a **2:1 aspect**
  (`roi_window` returns −1 otherwise), so tiles must be 2:1.
- **Content-addressed output.** `web_server.py`'s `cache_key` (quantized params
  → hash → file) is exactly the right pattern; reuse it. A tile that already
  exists in object storage is skipped, so a preempted instance costs at most one
  tile.
- **Upload per tile, not per job.** Nothing on local disk survives preemption.
- **Idempotent driver** that lists completed tiles at startup and renders the
  complement — then "resume" and "scale to 8 GPUs" are the same code path.
- **Pin `errormax`/`de0`/`SHADING_VERSION` into the key.** A mid-run parameter
  change must invalidate, not silently mix.
- **Convergence check first:** render one tile at errormax 1e-4 / 1e-5 / 1e-6 and
  compare. Whatever is spent on the grid is wasted if the tolerance was loose.

## 7. Measure these before committing

1. **Actual f64 ray throughput on an A100**, with `sugar_ki` at 10³ M. §4's
   estimate spans 3×; everything downstream scales off it.

   *Partial answer, measured on an A100-SXM4-40GB:* a 1280×640 unzoomed frame
   (819k rays) at `errormax 1e-8` takes **~1.05 s**, i.e. **~780k rays/s**.
   But that is the **f32** path at the stock `sugar_ki = 1.01`, where rays reach
   the escape sphere almost immediately. The §4 estimate assumed f64 *and* an
   escape sphere at 10³ M — far more steps per ray — so this is an upper bound
   on the real precompute rate, not a measurement of it. Still to measure: the
   same frame in f64, and with the distant escape sphere §3c argues for.
2. **Step-count distribution.** The `idokorlat >= int(1.0/errormax)` cap couples
   the step budget to the tolerance — at HQ (`errormax = 5e-5`) that is 20 000
   steps. Find out what fraction of rays hit the cap; those are the photon-ring
   rays that matter most.
3. **Real compression ratio** on the §3 payload. §4 assumes 5–10× from
   smoothness. If it is 2×, tier A is not shippable at all.
4. **Convergence of `sugar_ki`**: how much do sky directions actually move
   between 40 M and 10³ M? Sets whether §3c matters.

---

## 8. Unrelated bug found while checking whether orientation is precomputable — **fixed**

`cuda_ray.h:768` — the Rodrigues rotation matrix's middle row had misplaced
parentheses, folding both off-diagonal `sin(phi)` terms into the `(1 − cos)`
factor (and losing a sign):

```cpp
// was:
FP x2_tmp = (u[0] * u[1] * (1 - cos(phi) + u[2] * sin(phi))) * x[1] + ...
// correct:  u[1]*u[0]*(1 - cos(phi)) + u[2]*sin(phi)   and   u[1]*u[2]*(1 - cos(phi)) - u[0]*sin(phi)
```

The result was not a rotation — at 1.2 rad about a general axis, `det = 0.444`
and `‖RRᵀ − I‖ = 0.54`, so the camera basis was sheared and the ray directions
it produced were not unit vectors. Rows 1 and 3 were already correct.

Why it survived: it degenerates to the correct answer whenever
`u[0]·sin(phi)` and `u[2]·sin(phi)` both vanish, and `web_server.py`'s camera
never left that set by any single interaction. `BASE_ORIENTATION` is a 180° flip
about the local up axis; **panning** composes to another rotation about up
(`u[0] = u[2] = 0`), and a **tilt or roll on its own** tips the axis but leaves
the angle at exactly π (`sin(phi) = 0`). Only a pan *combined* with a tilt or
roll escapes both conditions — that is where the sheared frames came from.

Fixed, with `OpenACC_cli/tests/test_camera_rotation.cpp` covering it: it drives
the real `ijk_to_vec_zoom`, inverts the (diagonal, at a = 0) tetrad back out,
and asserts the recovered vector matches textbook Rodrigues and preserves norm
to 1e-14 across four axes and seven angles. Run with `make -C OpenACC_cli test`.
Confirmed to fail loudly against a copy with the old expression restored, and to
leave the default 180° view bit-identical.
