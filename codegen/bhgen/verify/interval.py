"""Rigorous domain-safety proofs by interval arithmetic with branch and bound.

Sampling cannot prove that a generated expression never divides by zero or takes
the square root of a negative number: those failures live on thin sets that
random points miss.  Interval arithmetic evaluates over whole boxes at once, so a
result that excludes zero is a statement about *every* point in the box, not
about the finitely many that were tried.

This is the stage that would have caught the shipped
``delta = r^2 - 4*rs*r + a^2 + Q^2`` bug directly: its range over the camera's
own parameter box contains negative values, and it sat under a sqrt.

Each obligation is discharged by recursive bisection: split the widest variable
until every leaf box is decided, a box is pruned as being outside the region of
interest, or the budget runs out.  Running out of budget is reported as UNKNOWN,
never as success.
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass
from typing import Sequence

import sympy as sp
from mpmath import iv


# --------------------------------------------------------------------------- #
#  evaluating a SymPy expression in interval arithmetic
# --------------------------------------------------------------------------- #

def _iv_eval(expr, env):
    """Evaluate `expr` with every free symbol bound to an mpmath interval."""
    if expr.is_Symbol:
        return env[expr]
    if expr.is_Integer:
        return iv.mpf(int(expr))
    if expr.is_Rational:
        return iv.mpf(int(expr.p)) / iv.mpf(int(expr.q))
    if expr.is_Float:
        return iv.mpf(str(expr))
    if expr is sp.pi:
        return iv.pi
    if expr.is_Add:
        acc = iv.mpf(0)
        for a in expr.args:
            acc = acc + _iv_eval(a, env)
        return acc
    if expr.is_Mul:
        acc = iv.mpf(1)
        for a in expr.args:
            acc = acc * _iv_eval(a, env)
        return acc
    if expr.is_Pow:
        base = _iv_eval(expr.base, env)
        if expr.exp.is_Integer:
            return base ** int(expr.exp)
        if expr.exp == sp.Rational(1, 2):
            return iv.sqrt(base)
        return iv.exp(_iv_eval(expr.exp, env) * iv.log(base))
    fn = {sp.sin: iv.sin, sp.cos: iv.cos, sp.tan: iv.tan,
          sp.exp: iv.exp, sp.log: iv.log, sp.sqrt: iv.sqrt}
    for key, f in fn.items():
        if isinstance(expr, key):
            return f(_iv_eval(expr.args[0], env))
    if isinstance(expr, sp.Abs):
        return abs(_iv_eval(expr.args[0], env))
    raise NotImplementedError(f"interval evaluation of {type(expr).__name__}: {expr}")


# --------------------------------------------------------------------------- #
#  proof obligations extracted from an expression graph
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class Obligation:
    """`expr` must satisfy `kind` everywhere on the region of interest."""
    kind: str          # "nonzero" or "nonneg"
    expr: sp.Expr
    origin: str        # human-readable: which subexpression demanded it


def collect_obligations(exprs: Sequence[sp.Expr]) -> list[Obligation]:
    """Every division and even root in `exprs` becomes a proof obligation."""
    seen, out = set(), []
    def walk(e, where):
        if e in seen or e.is_Atom:
            return
        seen.add(e)
        if e.is_Pow:
            if e.exp.is_negative:
                out.append(Obligation("nonzero", e.base, f"{where}: 1/({e.base})"))
            if e.exp.is_Rational and not e.exp.is_Integer and e.exp.q % 2 == 0:
                out.append(Obligation("nonneg", e.base, f"{where}: sqrt({e.base})"))
        if isinstance(e, sp.sqrt.__class__) or (e.func is sp.Pow and e.exp == sp.Rational(1, 2)):
            out.append(Obligation("nonneg", e.base, f"{where}: sqrt({e.base})"))
        for a in e.args:
            walk(a, where)
    for i, e in enumerate(exprs):
        walk(sp.sympify(e), f"expr[{i}]")
    # de-duplicate on (kind, expr)
    uniq = {}
    for o in out:
        uniq.setdefault((o.kind, o.expr), o)
    return list(uniq.values())


# --------------------------------------------------------------------------- #
#  branch and bound
# --------------------------------------------------------------------------- #

@dataclass
class Result:
    """Outcome of one obligation.

    ``PROVED``
        Holds everywhere in the region of interest.
    ``PROVED_MODULO_BOUNDARY``
        Holds everywhere except within `boundary_volume` of the constraint
        surface itself.  Exact branch and bound cannot terminate on a boundary of
        measure zero - boxes straddling it never resolve - so the honest
        statement is "proved on the region minus an epsilon shell", with the
        shell's volume reported so you can judge whether it matters.
    ``REFUTED``
        A box lying *strictly inside* the region was found on which the property
        provably fails.  `witness` is that box: a genuine counterexample, not a
        failure to decide.
    ``UNKNOWN``
        The budget ran out, or a box strictly inside the region stayed undecided
        down to `min_width`.  Never treat this as success.
    """

    status: str
    boxes_examined: int
    witness: dict | None = None
    margin: float | None = None
    boundary_volume: float = 0.0
    region_volume: float = 0.0


def _decide(kind, val):
    """Decide an obligation on one box, or return None for "split further"."""
    lo, hi = val.a, val.b
    if kind == "nonzero":
        if lo > 0 or hi < 0:
            return True
        return False if (lo == 0 and hi == 0) else None
    if kind == "nonneg":
        if lo >= 0:
            return True
        return False if hi < 0 else None
    raise ValueError(kind)


def _volume(cur):
    v = 1.0
    for lo, hi in cur:
        v *= float(hi - lo)
    return v


def discharged_by_constraints(obligation: Obligation, constraints) -> str | None:
    """Discharge an obligation from the constraints alone, with no search.

    If some constraint is ``expr + d`` for a numeric ``d <= 0``, then
    ``constraint > 0`` already gives ``expr > -d >= 0``, which settles both a
    "nonneg" and a "nonzero" obligation everywhere in the region of interest.

    This matters because the obligation and the constraint are frequently the
    same quantity: "Delta must be positive under the square root" versus "the
    camera is outside the horizon, where Delta is positive".  Without the
    shortcut, branch and bound would chew through the entire boundary shell to
    rediscover an implication that is true by inspection.
    """
    for c in constraints:
        for expr, sign in ((obligation.expr, "+"), (-obligation.expr, "-")):
            # For a "nonzero" obligation the sign is irrelevant: expr != 0 iff
            # -expr != 0, and the constraint may have been written either way
            # round (r - rs > 0 discharges 1/(rs - r) just as well as 1/(r - rs)).
            if sign == "-" and obligation.kind != "nonzero":
                continue
            d = sp.simplify(sp.sympify(c) - expr)
            if d.is_number and d.is_real and d <= 0:
                return f"implied by constraint ({c}) > 0"
    return None


def _ivt_refutes(expr, cur, syms):
    """Rigorously exhibit a zero of `expr` inside `cur` by the intermediate value theorem.

    A "nonzero" obligation can essentially never be refuted by interval evaluation
    alone: that would require an interval that is exactly [0, 0].  But if two
    points of the box give values of strictly opposite sign - each established by
    evaluating a degenerate interval, so each sign is certain - then continuity
    puts a genuine zero between them.  That is a proof, not a heuristic.
    """
    corners = []
    for pick in (0, 1):
        env = {s: iv.mpf([c[pick], c[pick]]) for s, c in zip(syms, cur)}
        try:
            corners.append(_iv_eval(expr, env))
        except (ZeroDivisionError, ValueError):
            return False
    (lo_v, hi_v) = corners
    return (lo_v.b < 0 and hi_v.a > 0) or (lo_v.a > 0 and hi_v.b < 0)


def prove(obligation: Obligation, box: dict, constraints: Sequence[sp.Expr] = (),
          max_boxes: int = 20000, min_width: float = 1e-9,
          boundary_width: float = 1e-3) -> Result:
    """Try to prove `obligation` on `box`, subject to every constraint being > 0.

    `constraints` carve the region the ray tracer actually visits out of the
    parameter box - for a black hole, "outside the horizon".  A sub-box is
    classified three ways:

    * every constraint provably non-positive somewhere -> outside, pruned;
    * every constraint provably positive              -> strictly inside, so a
      failure here is a real counterexample;
    * otherwise                                        -> straddles the boundary,
      and is split until it is narrow enough to charge to `boundary_volume`.

    Search is best-first, ordered so that the most promising box for a
    *refutation* is examined first: a real bug is found in a handful of boxes
    rather than after exhausting the tree.

    Two width thresholds, because the two kinds of undecided box mean different
    things.  A box straddling the constraint surface is charged to
    `boundary_volume` once it is narrower than `boundary_width` - splitting it
    further only ever resolves a measure-zero surface and would not terminate.  A
    box strictly *inside* the region is pursued all the way down to `min_width`
    before giving up, because there an undecided result is a real gap in the
    proof rather than an artefact of the region's edge.
    """
    # Only split on variables the obligation or a constraint actually mentions.
    # Bisecting a coordinate the expression does not depend on - t and phi are
    # usually absent from a stationary, axisymmetric metric - doubles the search
    # for nothing.
    relevant = set(obligation.expr.free_symbols)
    for c in constraints:
        relevant |= sp.sympify(c).free_symbols
    syms = [s for s in box if s in relevant] or list(box)
    start = tuple(box[s] for s in syms)
    constraints = [sp.sympify(c) for c in constraints]

    why = discharged_by_constraints(obligation, constraints)
    if why is not None:
        return Result("PROVED", 0, margin=None)

    counter = 0
    heap = [(0.0, counter, start)]
    examined = 0
    boundary_volume = 0.0
    interior_volume = 0.0
    worst = None

    while heap:
        if examined >= max_boxes:
            _, _, cur = heap[0]
            return Result("UNKNOWN", examined,
                          witness=dict(zip((str(s) for s in syms), cur)),
                          boundary_volume=boundary_volume, region_volume=interior_volume)
        _, _, cur = heapq.heappop(heap)
        examined += 1
        env = {s: iv.mpf([lo, hi]) for s, (lo, hi) in zip(syms, cur)}

        outside = inside = False
        try:
            cvals = [_iv_eval(c, env) for c in constraints]
            outside = any(cv.b <= 0 for cv in cvals)
            inside = all(cv.a > 0 for cv in cvals)
        except (ZeroDivisionError, ValueError):
            cvals, inside = [], False
        if outside:
            continue

        try:
            val = _iv_eval(obligation.expr, env)
        except (ZeroDivisionError, ValueError):
            val = None
        verdict = _decide(obligation.kind, val) if val is not None else None

        if verdict is True:
            m = float(min(abs(val.a), abs(val.b)) if obligation.kind == "nonzero" else val.a)
            worst = m if worst is None else min(worst, m)
            if inside:
                interior_volume += _volume(cur)
            continue
        if verdict is False and inside:
            return Result("REFUTED", examined,
                          witness=dict(zip((str(s) for s in syms), cur)),
                          boundary_volume=boundary_volume, region_volume=interior_volume)
        if (verdict is None and inside and obligation.kind == "nonzero"
                and _ivt_refutes(obligation.expr, cur, syms)):
            return Result("REFUTED", examined,
                          witness=dict(zip((str(s) for s in syms), cur)),
                          boundary_volume=boundary_volume, region_volume=interior_volume)

        widths = [hi - lo for lo, hi in cur]
        k = max(range(len(cur)), key=lambda i: widths[i])
        limit = min_width if inside else boundary_width
        if widths[k] <= limit:
            if inside:
                return Result("UNKNOWN", examined,
                              witness=dict(zip((str(s) for s in syms), cur)),
                              boundary_volume=boundary_volume, region_volume=interior_volume)
            boundary_volume += _volume(cur)
            continue

        lo, hi = cur[k]
        mid = (lo + hi) / 2
        for half in ((lo, mid), (mid, hi)):
            child = cur[:k] + (half,) + cur[k + 1:]
            # Priority: chase the most negative lower bound, i.e. the box most
            # likely to contain a violation.
            key = float(val.a) if val is not None else -1e300
            counter += 1
            heapq.heappush(heap, (key, counter, child))

    status = "PROVED" if boundary_volume == 0.0 else "PROVED_MODULO_BOUNDARY"
    return Result(status, examined, margin=worst,
                  boundary_volume=boundary_volume, region_volume=interior_volume)
