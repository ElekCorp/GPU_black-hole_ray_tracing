# Physics regression tests

`geodesic_invariants` checks the two things that must hold for any correct
null-geodesic tracer in Kerr-Newman, and that were both broken before:

1. **The ray leaves the camera null**, `g_{mu nu} v^mu v^nu == 0`.  This tests
   `ijk_to_vec_zoom()` — the Rodrigues rotation and Carter's tetrad.
2. **`g(v,v)`, the energy `E` and the angular momentum `L` are conserved along
   the ray.**  `d/dt` and `d/dphi` are Killing vectors of the metric, so `E` and
   `L` are exact constants of the motion; any drift is integration error.  This
   tests `RK45()` and the step controller.

Run with `make -C tests run` (host `g++`, no GPU needed).

Measured on the code as of this commit, against the same 24 rays:

| quantity                        | before  | after   |
| ------------------------------- | ------- | ------- |
| `\|g(v,v)\|` at the camera      | 9.1e-01 | 5.6e-16 |
| `\|g(v,v)\|` after integration  | 9.1e-01 | 3.1e-11 |
| relative drift in `E`           | 3.1e-06 | 1.3e-11 |
| relative drift in `L`           | 7.7e-05 | 7.7e-11 |

---

`metric_codegen_equiv` is the bridge between the hand-written tracer and the
`codegen/` pipeline. It checks three things about the generated Kerr-Newman
header:

1. its `geodesic_rhs()` agrees with the `christoffel()` block the tracer ships;
2. its `metric()` agrees with the closed-form Boyer-Lindquist metric;
3. the two are mutually consistent, via `d/de [g(v,v)] = 0` along the geodesic
   flow with the metric derivative taken numerically.

The shipped block was independently verified against symbolically derived
Kerr-Newman Christoffel symbols, so (1) transfers that confidence to the
generated code — and any future disagreement localises immediately to the
generator rather than to the tracer.
