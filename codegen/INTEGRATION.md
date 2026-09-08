# Wiring a generated metric into the OpenACC tracer

The generator and its proofs are complete and the headers are committed; this is
the remaining step, written down so it can be picked up directly.

The tracer currently hard-codes Kerr-Newman in three places. Each maps onto one
generated entry point.

| today, in `OpenACC_cli/cuda_ray.h` | replacement |
| --- | --- |
| `christoffel()` — the pasted 84-temporary block | `bh_<name>::geodesic_rhs(p, x, v, acc)` |
| the Carter tetrad written out inside `ijk_to_vec_zoom()` | host-side Gram-Schmidt from `metric()` + `observer()`, uploaded as a 4x4 |
| `gomb_be` / `gomb_ki` / `disk1` / `disk2` | `event_fields()` + `event_guards()` against `EVENT_KIND` |

### 1. Geodesic right-hand side

`RK45()` already calls `christoffel(hole, xs, vs, kv[s])` through a single seam,
so this is a one-line change plus threading the parameter vector through
`kerr_black_hole`. `tests/metric_codegen_equiv.cpp` pins the generated block to
the shipped one first, so the swap is checkable before and after.

### 2. Camera frame, on the host

All rays share one camera position, so the frame is one 4x4 per rendered frame
rather than per ray. Build it in `makeframe_T()`:

```
metric(p, x_camera, g);
observer(p, x_camera, u);          // or the ZAMO fallback, from g
gram_schmidt(g, u, SIGNATURE, e);  // 4x4, host side
assert_orthonormal(g, e, SIGNATURE);
```

then upload `e` and reduce the device path to `v[mu] = sum_b e[b][mu] * dir[b]`.
This deletes `ijk_to_vec_zoom`'s frame algebra outright, along with the class of
bug that produced the wrong frame in the first place. The runtime assertion is
cheap because it runs once per frame, not once per ray.

### 3. Events

Replace the four hard-coded predicates with a loop over `N_EVENT`, comparing
`event_fields` at the step endpoints according to `EVENT_KIND` and requiring the
event's `event_guards` slice to be positive. Two things fall out for free:

* **Root refinement.** The tracer currently reports the step endpoint as the hit
  point, so disk-edge accuracy is bounded by the step size. With a scalar event
  field in hand, bisecting on `f` between the endpoints gives the crossing to
  machine precision for a handful of extra evaluations.
* **Chart independence.** `disk1`/`disk2` test `theta` against `pi/2`, which only
  means anything in a spherical chart. The event formulation carries over to the
  Cartesian and Kerr-Schild metrics unchanged.

### 4. Selecting a metric

The generated headers are independent namespaces, so the cheapest approach is a
compile-time choice — `-DBH_METRIC=kerr_newman` with a small dispatch header —
which keeps `geodesic_rhs` inlinable into the integrator. A runtime switch would
cost the inlining, which is most of the performance on the GPU path.

`PARAM_NAME` / `PARAM_DEFAULT` / `PARAM_DOC` are emitted precisely so
`cli_parser.cpp` can build its option table from the metric instead of hard-coding
`--a`, `--Q`, `--rs`.
