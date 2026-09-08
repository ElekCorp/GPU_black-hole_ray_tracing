"""Floating-point round-off proofs for the generated code, via Gappa.

Everything else in this package reasons about real numbers.  This module reasons
about the numbers the GPU actually computes.  It hands Gappa the exact
straight-line program the emitter produces - same CSE temporaries, same
association, same integer-power expansion - and asks it to bound

    | fl(program)  -  program evaluated in exact arithmetic |

over the input box.  Gappa answers with a machine-checkable proof.

What is and is not covered
--------------------------
Covered: every addition, subtraction, multiplication, division and square root,
including the catastrophic cancellation that shows up near Delta = 0.

Not covered: the accuracy of libm itself.  Gappa has no model of sin, cos or
tan, so each transcendental application is replaced by a fresh input variable
constrained to its true range, and the exact reference chain uses the *same*
variable.  The bound returned is therefore the round-off of the arithmetic given
the transcendental results, which is the part the code generator is responsible
for; libm's own error is a separate, documented input to any end-to-end budget.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass

import sympy as sp

#: Mirrors OpenACC_cli/cuda_ray.h's pown(): small exponents are spelled out, so
#: the proof has to associate the multiplications the same way the code does or
#: it is bounding a different program.
def pown_expansion(base: str, n: int) -> str:
    if n < 0:
        return f"(1 / {pown_expansion(base, -n)})"
    if n == 0:
        return "1"
    if n == 1:
        return base
    if n == 2:
        return f"({base} * {base})"
    if n == 3:
        return f"({base} * {base} * {base})"
    if n == 4:
        return f"(({base} * {base}) * ({base} * {base}))"
    if n == 5:
        return f"((({base} * {base}) * ({base} * {base})) * {base})"
    if n == 6:
        return f"(({base} * {base} * {base}) * ({base} * {base} * {base}))"
    # generic binary exponentiation, as pown_gen() does it
    out, sq, e = None, base, n
    while e:
        if e & 1:
            out = sq if out is None else f"({out} * {sq})"
        sq = f"({sq} * {sq})"
        e >>= 1
    return out


class GappaPrinter:
    """Print a SymPy expression in Gappa's expression syntax."""

    def __init__(self, names: dict, aux: dict | None = None):
        self.names = dict(names)    # Symbol -> gappa identifier
        # Shared between the rounded and the exact printer so both chains refer
        # to the *same* libm result, isolating arithmetic round-off.
        self.aux: dict = {} if aux is None else aux
        self._n = 0

    def _fresh(self, stem: str) -> str:
        self._n = len(self.aux) + 1
        return f"{stem}{self._n}"

    def __call__(self, e) -> str:
        e = sp.sympify(e)
        if e in self.names:
            return self.names[e]
        if e.is_Integer:
            return str(int(e))
        if e.is_Rational:
            return f"({e.p} / {e.q})"
        if e.is_Float:
            return repr(float(e))
        if e is sp.pi:
            return "3.14159265358979323846"
        if e.is_Add:
            return "(" + " + ".join(self(a) for a in e.args) + ")"
        if e.is_Mul:
            num, den = e.as_numer_denom()
            if den != 1:
                return f"({self(num)} / {self(den)})"
            return "(" + " * ".join(self(a) for a in e.args) + ")"
        if e.is_Pow:
            if e.exp == sp.Rational(1, 2):
                return f"sqrt({self(e.base)})"
            if e.exp.is_Integer:
                return pown_expansion(self(e.base), int(e.exp))
            raise NotImplementedError(f"gappa: non-integer power {e}")
        if isinstance(e, (sp.sin, sp.cos, sp.tan, sp.exp, sp.log, sp.asin, sp.acos, sp.atan)):
            # Modelled as an opaque input; see the module docstring.
            if e not in self.aux:
                self.aux[e] = self._fresh("libm_")
            return self.aux[e]
        if isinstance(e, sp.Abs):
            return f"|{self(e.args[0])}|"
        raise NotImplementedError(f"gappa: {type(e).__name__}: {e}")


@dataclass
class RoundoffResult:
    """Outcome of one round-off obligation.

    ``status``
        ``TIGHT``  - the bound is small compared with the output's own magnitude,
                     i.e. Gappa genuinely tracked the cancellation.
        ``LOOSE``  - a bound was returned but it is no better than bounding the
                     two chains separately, so it says nothing useful about
                     round-off.  Usually means Gappa needs a subdivision hint.
        ``NO_BOUND``- Gappa could not bound the expression at all.
        ``ERROR``  - the tool was missing or timed out.
    """

    target: str
    status: str
    abs_error: float | None
    raw: str
    magnitude: float | None = None
    relative: float | None = None


class UnboundedLibm(ValueError):
    """A transcendental stand-in has no finite range, so no bound is provable."""


def _finite(x) -> bool:
    return x == x and abs(x) != float("inf")


_FMT = {"double": "float<ieee_64,ne>", "float": "float<ieee_32,ne>"}
_EPS = {"double": 2.0 ** -53, "float": 2.0 ** -24}


def _range_over_box(expr, input_ranges, reps=()):
    """Bound `expr` over the input box by interval arithmetic.

    `reps` is the CSE chain.  It has to be threaded through because CSE hoists
    subexpressions out of transcendental calls - an angle becomes `sin(t25)` -
    and a bare input box says nothing about `t25`.  Evaluating the chain in order
    gives every temporary an interval, so the stand-in for `sin(t25)` gets a real
    range instead of failing.

    Used for two things: the range hypothesis Gappa needs for each libm
    stand-in, and the magnitude the absolute round-off bound is judged against.

    Giving sin(x2) its true range matters.  Left at the generic [-1, 1] it
    admits sin = 0, so any 1/sin in the expression is unbounded and Gappa
    correctly refuses to prove anything at all.
    """
    from mpmath import iv
    from .interval import _iv_eval
    env = {s: iv.mpf([lo, hi]) for s, (lo, hi) in input_ranges.items()}
    try:
        for sym, val in reps:
            env[sym] = _iv_eval(sp.sympify(val), env)
        v = _iv_eval(sp.sympify(expr), env)
        return float(v.a), float(v.b)
    except Exception:
        return None


def build_script(reps, reduced, targets, input_ranges, precision="double",
                 goal_index=0, subdivide: dict | None = None) -> tuple[str, GappaPrinter]:
    """Emit a Gappa script bounding the round-off of one output.

    `reps`/`reduced` are exactly what :func:`sympy.cse` returned for the emitted
    function, so the proof tracks the shipped program rather than an idealised
    rewriting of it.
    """
    names = {s: str(s).replace("[", "_").replace("]", "") for s in input_ranges}

    # Two printers over one shared libm dictionary: the rounded chain calls the
    # temporaries t0, t1, ...; the exact chain calls them E_t0, E_t1, ...  Naming
    # them apart rather than rewriting text afterwards keeps the two programs
    # provably structurally identical.
    aux: dict = {}
    pr = GappaPrinter(names, aux)
    pr_e = GappaPrinter(names, aux)

    rounded, body_exact = [], []
    for sym, val in reps:
        rounded.append(f"{sym} rnd= {pr(val)};")
        body_exact.append(f"E_{sym} = {pr_e(val)};")
        pr.names[sym] = str(sym)
        pr_e.names[sym] = f"E_{sym}"

    tgt_r = pr(reduced[goal_index])
    tgt_e = pr_e(reduced[goal_index])

    lines = [f"@rnd = {_FMT[precision]};", ""]
    for sym, name in names.items():
        lines.append(f"{name} = rnd(dummy_{name});")
    for expr, name in pr.aux.items():
        lines.append(f"{name} = rnd(dummy_{name});")
    lines += ["", "# --- the program as executed, every operation rounded ---"]
    lines += rounded
    lines.append(f"out rnd= {tgt_r};")
    lines += ["", "# --- the same program in exact arithmetic ---"]
    lines += body_exact
    lines.append(f"E_out = {tgt_e};")

    hyp = []
    for s, (lo, hi) in input_ranges.items():
        hyp.append(f"{names[s]} in [{lo!r}, {hi!r}]")
    for expr, name in aux.items():
        rng = _range_over_box(expr, input_ranges, reps)
        if rng is None or not all(map(_finite, rng)):
            raise UnboundedLibm(
                f"{expr} is unbounded on the input box (range {rng}); "
                f"rewrite it via bhgen.normalize before emitting")
        hyp.append(f"{name} in [{rng[0]!r}, {rng[1]!r}]")
    lines += ["", "{ " + " /\\\n  ".join(hyp) + "\n  -> |out - E_out| in ? }"]
    if subdivide:
        # Gappa often cannot see through cancellation on a wide input range.
        # Splitting the dominant variable lets it bound each slice separately;
        # the result is still a proof over the whole range, just a costlier one.
        for sym, n in subdivide.items():
            lines.append(f"$ {names[sym]} in {n};")
    return "\n".join(lines) + "\n", pr


# Gappa reports either a plain number ("[0, 0]") or its own binary format with a
# decimal gloss in braces ("[0, 3b-52 {6.66134e-16, 2^(-50.415)}]").  Gappa also
# rewrites the goal's left-hand side when the two chains coincide, so match on
# the "- E_out|" tail rather than on a fixed name.
_GOAL = re.compile(r"\|\S+ - E_out\| in \[([^,]+),\s*([^\]]+)\]")


def _parse_bound(text: str):
    m = _GOAL.search(text)
    if not m:
        return None
    hi = m.group(2).strip()
    brace = re.search(r"\{([0-9.eE+-]+)", hi)
    if brace:
        return float(brace.group(1))
    try:
        return float(hi)
    except ValueError:
        return None


def classify(res: "RoundoffResult", magnitude: float | None, precision: str) -> "RoundoffResult":
    """Judge an absolute bound against the output's magnitude on the same box."""
    if res.status != "PROVED" or magnitude is None or res.abs_error is None:
        return res
    if res.abs_error == 0.0:
        return RoundoffResult(res.target, "TIGHT", 0.0, res.raw, magnitude, 0.0)
    rel = res.abs_error / magnitude if magnitude > 0 else float("inf")
    # A meaningful round-off bound has to be far below the value itself.  Anything
    # near the magnitude is just |out| + |E_out| in disguise.
    status = "TIGHT" if rel < 1e-6 else "LOOSE"
    return RoundoffResult(res.target, status, res.abs_error, res.raw, magnitude, rel)


def run(script: str, target: str = "out", timeout: int = 600) -> RoundoffResult:
    if shutil.which("gappa") is None:
        return RoundoffResult(target, "ERROR", None, "gappa not installed")
    with tempfile.NamedTemporaryFile("w", suffix=".gappa", delete=False) as fh:
        fh.write(script)
        path = fh.name
    try:
        p = subprocess.run(["gappa", path], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return RoundoffResult(target, "ERROR", None, "gappa timed out")
    out = p.stdout + p.stderr
    bound = _parse_bound(out)
    if bound is not None:
        return RoundoffResult(target, "PROVED", bound, out)
    return RoundoffResult(target, "NO_BOUND", None, out)
