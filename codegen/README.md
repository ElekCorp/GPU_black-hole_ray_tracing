# Metric codegen and verification pipeline

Give it a metric; it gives you a verified C++ header the ray tracer can use.

```
codegen/metrics/<name>.py          you write this: g_{mu nu}, the camera's
                                   observer, and the scene's stopping events
        |
        |  make -C codegen
        v
OpenACC_cli/generated/bh_<name>_metric.h    generated, committed, and only
                                            written if every proof obligation
                                            that must hold does hold
```

## What is generated

All as `template <class FP>`, so one header serves both the float and the double
build:

| function | meaning |
| --- | --- |
| `metric(p, x, g)` | `g_{mu nu}` at `x` |
| `geodesic_rhs(p, x, v, acc)` | `acc^mu = -Gamma^mu_{ab} v^a v^b` |
| `observer(p, x, u)` | the camera's 4-velocity, up to scale |
| `event_fields(p, x, f)` | one scalar per stopping condition |
| `event_guards(p, x, gd)` | guard scalars for those conditions |

`p` is one flat vector holding the metric's parameters and the scene's alike, so
an event may mention both. `PARAM_NAME` / `PARAM_DEFAULT` describe it to the CLI.

## Any chart, not just (t, r, theta, phi)

Two design choices carry the "works with every given metric" requirement.

**Events are sign conditions on scalar fields.** "The ray hit the disk" is not
`theta` crossing `pi/2` — that sentence only parses in a spherical chart. It is
*a scalar field changing sign while its guard scalars stay positive*, which
parses in any chart. In Boyer-Lindquist the disk field is `theta - pi/2`; in the
Cartesian example metric it is `z`; the tracer cannot tell the difference. Same
for the horizon: Kerr-Newman declares it as `Delta`, which is the metric's own
function, so a naked singularity (`rs^2 < 4(a^2+Q^2)`, where `Delta` never
vanishes) needs no special case — the event simply never fires.

**The camera frame is built numerically on the host, not written out by hand.**
The camera occupies a single point, so its orthonormal frame is one 4x4 matrix
per rendered frame. The pipeline builds it by Lorentzian Gram-Schmidt from
`metric()` and `observer()` and uploads sixteen numbers; the device never
contains metric-specific frame algebra. That is a direct response to how the
previous frame was wrong: it was a closed-form Carter tetrad, mistranscribed. A
frame that is constructed and then checked against `e^T g e = eta` cannot be
mistranscribed.

Seeding Gram-Schmidt with Carter's observer reproduces Carter's frame exactly for
Kerr, so images are unchanged. Metrics that do not name an observer fall back to
the zero-angular-momentum observer, which stays timelike inside an ergosphere
where a static observer does not.

## Verification

Run automatically at generation time. **A header is written only if nothing is
refuted**, and the report is embedded in the header's own banner.

| stage | what it establishes | strength |
| --- | --- | --- |
| Butcher order conditions (`verify/butcher.py`) | the RK tableau has the order it claims | exact rational arithmetic |
| Cross-derivation (`verify/symbolic.py`) | the geodesic RHS from the Christoffel formula equals the one from the Euler-Lagrange equations | exact identity (see below) |
| CSE equivalence | the emitted straight-line program computes the expression it was derived from | exact identity |
| Tetrad orthonormality | `e^T g e - eta == 0` | exact identity |
| Domain safety (`verify/interval.py`) | no division by zero, no square root of a negative, anywhere in the region the tracer visits | interval arithmetic + branch and bound, rigorous over a continuum |
| Round-off (`verify/roundoff.py`) | a proven bound on `|fl(program) - program|` | Gappa proof |

Three of these are proofs of identities, so they hold for every camera position
and every parameter value at once — something no amount of sampling can give.
The interval stage is rigorous over whole boxes rather than sample points. The
Gappa stage reasons about the floating-point numbers the GPU actually computes.

### Deciding the identities (`verify/algebra.py`)

`sympy.simplify` is a heuristic. It can run for many minutes on a Kerr-sized
expression, and it can fail to reach `0` even when the identity holds — which
would surface as a spurious REFUTED. So the identity obligations use a decision
procedure instead, strongest strategy first:

1. **Rational normal form, then polynomial reduction modulo `sin^2 + cos^2 - 1`.**
   Every expression here is a rational function of the coordinates, parameters,
   velocities and the sines and cosines of the angles. Rewrite `tan`/`cot`/`sec`/
   `csc` as `sin`/`cos` quotients, expand multiple angles, replace each `sin u`
   and `cos u` by a symbol, and reduce the numerator modulo the relations
   `s^2 + c^2 = 1`. A zero remainder is a proof: polynomial reduction is exact
   and terminating.
2. **`simplify`**, for the tetrad, whose entries carry square roots and so are
   not rational functions.
3. **Exact-rational Schwartz-Zippel**, if the first two exhaust their budget.
   Reported as `probable`, never as `proved`.

A cheap Schwartz-Zippel screen runs first regardless, so a genuinely wrong
expression is refuted in milliseconds rather than after a fruitless reduction.
The sampling is subtle enough to be worth stating: substituting a random rational
for `theta` would leave `sin(1234/5678)` behind, which exact arithmetic cannot
decide, so the sines and cosines are instead given random *rational points on the
unit circle* through the Weierstrass substitution `s = 2t/(1+t^2)`,
`c = (1-t^2)/(1+t^2)`. That respects the one relation tying them together and
keeps every evaluation rational — no floating point anywhere, so a non-zero
result is a definitive refutation rather than a rounding artefact.

### Why these and not a proof assistant

Formalising pseudo-Riemannian geometry in Lean or Coq and proving generated C++
against it is a multi-person-year effort with no existing GR library to build on.
The one part where a proof assistant would genuinely pay — the Runge-Kutta order
conditions — is finite exact rational arithmetic, which `verify/butcher.py`
settles in a few dozen lines at the same rigour. SMT (dReal, Z3) would suit the
small obligations but chokes on an eighty-temporary expression block.

### Honest limits

* The identity proofs rest on SymPy's polynomial arithmetic being correct — the
  same standing assumption every computer-algebra derivation in this field makes.
  Where the procedure falls back to Schwartz-Zippel the result is labelled
  `probable`, with a failure probability of at most `(d/10^6)^24`.
* Interval branch and bound cannot terminate on a constraint surface of measure
  zero. Boxes straddling it are charged to `boundary_volume` and reported, so the
  claim is "proved on the region minus a shell of this size", never a bare
  "proved". Budget exhaustion reports `UNKNOWN`, never success.
* The search is exponential in the number of variables that actually appear.
  Pinning the metric parameters to the values being rendered, rather than leaving
  them as free ranges, is the difference between `UNKNOWN` and a proof.
* Gappa's bound covers arithmetic round-off only. It has no model of `sin` or
  `cos`, so each transcendental application becomes a fresh input variable with
  its true range (computed by the interval stage) and the exact reference chain
  uses the same variable. libm's own error is a separate, documented input.
* A `LOOSE` round-off result means Gappa returned a bound no better than
  evaluating the two chains separately. It is reported as weak, not as a pass.
* The absolute round-off bounds are judged against a magnitude that is itself an
  interval bound, and therefore an overestimate. Read the relative figure as an
  upper bound on the relative error, not as a tight estimate.

## Adding a metric

Write `codegen/metrics/<name>.py` exporting `SPEC` (a `MetricSpec`) and,
optionally, `VERIFY_BOX`, `VERIFY_CONSTRAINTS` and `ROUNDOFF_BOX` describing the
region to verify over. `minkowski_cartesian.py` is the smallest complete example
and deliberately uses a Cartesian chart. Then `make -C codegen`.

## Requirements

`make venv` creates the Python environment (sympy, mpmath, numpy). The round-off
stage additionally needs Gappa (`apt-get install gappa`); `make fast` skips it.
