"""SymPy expressions -> a CSE'd, precision-generic C++ header.

The printer here differs from a plain `sympy.ccode` call in three ways that
matter for this project:

* Every numeric literal is wrapped as ``FP(...)``.  Emitting a bare ``1.0/8.0``
  inside a ``template <class FP>`` silently promotes the whole subexpression to
  double, which costs roughly half the throughput on the float path.
* Integer powers become ``pown(x, n)``, the repo's binary-exponentiation helper,
  rather than a libm ``pow`` call.
* All the expressions of one function are passed to ``cse`` together, so the
  four geodesic components share their common subexpressions the way the
  hand-generated block did.
"""

from __future__ import annotations

import textwrap
from typing import Iterable, Sequence

import sympy as sp
from sympy.printing.c import C99CodePrinter

from .normalize import normalize
from .spec import DIM, EventKind, MetricSpec

#: CSE temporaries are named t0, t1, ...  Coordinates and velocities are printed
#: as x[i] / v[i], so there is no way for the two to collide.
_TEMP = sp.numbered_symbols("t")


class FPPrinter(C99CodePrinter):
    """C99 printer that keeps every literal in the template's own precision."""

    def _fp(self, text: str) -> str:
        return f"FP({text})"

    def _print_Float(self, expr):
        return self._fp(super()._print_Float(expr))

    def _print_Rational(self, expr):
        return self._fp(f"{expr.p}.0 / {expr.q}.0")

    def _print_Integer(self, expr):
        # Integers are exact in every format and promote cleanly, so leave them
        # bare; wrapping them would only add noise.
        return str(expr.p)

    def _print_Pi(self, expr):
        return self._fp("3.14159265358979323846")

    def _print_Exp1(self, expr):
        return self._fp("2.71828182845904523536")

    def _print_Half(self, expr):
        return self._fp("1.0 / 2.0")

    def _print_Pow(self, expr):
        base, exp = expr.base, expr.exp
        if exp == sp.Rational(1, 2):
            return f"sqrt({self._print(base)})"
        if exp == sp.Rational(-1, 2):
            return f"(FP(1) / sqrt({self._print(base)}))"
        if exp.is_Integer:
            n = int(exp)
            # pown() handles negative exponents, but only through its generic
            # binary-exponentiation loop.  A reciprocal is both clearer and
            # cheaper, and -1 is by far the most common case in these expressions.
            if n == -1:
                return f"(FP(1) / {self._print(base)})"
            if n < 0:
                return f"(FP(1) / pown({self._print(base)}, {-n}))"
            return f"pown({self._print(base)}, {n})"
        return f"pow({self._print(base)}, {self._print(exp)})"


def _printer() -> FPPrinter:
    return FPPrinter()


def _subs_map(spec: MetricSpec, with_velocity: bool) -> dict:
    """Coordinates, parameters and velocities -> C array element symbols.

    The replacement symbols are named ``x[0]``, ``p[3]`` and so on.  SymPy treats
    a symbol name as opaque, and the printer emits it verbatim, so this produces
    correct C without any string post-processing.
    """
    out = {c: sp.Symbol(f"x[{i}]") for i, c in enumerate(spec.coords)}
    out.update({p.symbol: sp.Symbol(f"p[{i}]") for i, p in enumerate(spec.params)})
    if with_velocity:
        out.update({sp.Symbol(f"v{i}", real=True): sp.Symbol(f"v[{i}]") for i in range(DIM)})
    return out


def render_body(exprs: Sequence[sp.Expr], targets: Sequence[str],
                subs: dict, indent: str = "    ") -> str:
    """CSE a group of expressions and render them as a C++ function body."""
    exprs = [normalize(sp.sympify(e)).xreplace(subs) for e in exprs]
    reps, reduced = sp.cse(exprs, symbols=sp.numbered_symbols("t"),
                           optimizations="basic")
    pr = _printer()
    lines = []
    for sym, val in reps:
        lines.append(f"{indent}FP const {sym} = {pr.doprint(val)};")
    if reps:
        lines.append("")
    for tgt, val in zip(targets, reduced):
        lines.append(f"{indent}{tgt} = {pr.doprint(val)};")
    return "\n".join(lines), len(reps)


def _fn(signature: str, doc: str, body: str) -> str:
    doc = textwrap.indent(textwrap.dedent(doc).strip("\n"), " * ")
    return f"/**\n{doc}\n */\ntemplate <class FP>\ninline void {signature}\n{{\n{body}\n}}\n"


def generate_header(spec: MetricSpec, acc: Sequence[sp.Expr],
                    frame: Sequence[Sequence[sp.Expr]] | None = None,
                    provenance: str = "") -> str:
    """Render the complete generated header for `spec`.

    `acc` is the geodesic acceleration (passed in rather than recomputed so the
    caller can reuse the copy the verification stage already proved correct).
    """
    ns = f"bh_{spec.name}"
    sub_x = _subs_map(spec, with_velocity=False)
    sub_xv = _subs_map(spec, with_velocity=True)
    stats = {}

    # ---- metric -----------------------------------------------------------
    gexprs, gtargets = [], []
    for i in range(DIM):
        for j in range(i, DIM):
            gexprs.append(spec.g[i, j])
            gtargets.append(f"g[{i}][{j}]")
    body, n = render_body(gexprs, gtargets, sub_x)
    body += "\n\n" + "\n".join(
        f"    g[{j}][{i}] = g[{i}][{j}];" for i in range(DIM) for j in range(i + 1, DIM))
    stats["metric"] = n
    f_metric = _fn(
        "metric(FP const* p, FP const* x, FP g[4][4])",
        """
        The metric g_{mu nu} at x.  Symmetric; only the upper triangle is
        computed and the rest is mirrored.
        """, body)

    # ---- geodesic RHS -----------------------------------------------------
    body, n = render_body(acc, [f"acc[{i}]" for i in range(DIM)], sub_xv)
    stats["geodesic_rhs"] = n
    f_geo = _fn(
        "geodesic_rhs(FP const* p, FP const* x, FP const* v, FP* acc)",
        """
        Geodesic acceleration acc^mu = -Gamma^mu_{alpha beta} v^alpha v^beta,
        i.e. the right-hand side of dv^mu / de = acc^mu.

        The Christoffel symbols are never formed: the full double contraction is
        simplified and common-subexpression eliminated as one block.
        """, body)

    # ---- observer ---------------------------------------------------------
    f_obs = ""
    if spec.observer is not None:
        body, n = render_body(list(spec.observer), [f"u[{i}]" for i in range(DIM)], sub_x)
        stats["observer"] = n
        f_obs = _fn(
            "observer(FP const* p, FP const* x, FP* u)",
            """
            The camera's 4-velocity at x, up to normalisation - the time leg the
            frame builder is seeded with.  Only its direction matters.
            """, body)

    # ---- events -----------------------------------------------------------
    ev_exprs = [e.field for e in spec.events]
    body, n = render_body(ev_exprs, [f"f[{i}]" for i in range(len(ev_exprs))], sub_x)
    stats["event_fields"] = n
    f_ev = _fn(
        "event_fields(FP const* p, FP const* x, FP* f)",
        """
        The scalar field of every event, evaluated at x.  An event fires on the
        sign or sign change of its scalar; see EVENT_KIND.
        """, body)

    guard_exprs, guard_off, k = [], [], 0
    for e in spec.events:
        guard_off.append(k)
        guard_exprs.extend(e.guards)
        k += len(e.guards)
    if guard_exprs:
        body, n = render_body(guard_exprs, [f"gd[{i}]" for i in range(len(guard_exprs))], sub_x)
    else:
        body, n = "    (void)p; (void)x; (void)gd;", 0
    stats["event_guards"] = n
    f_gd = _fn(
        "event_guards(FP const* p, FP const* x, FP* gd)",
        """
        Guard scalars, concatenated over all events in order.  Event i owns the
        slice [EVENT_GUARD_OFFSET[i], EVENT_GUARD_OFFSET[i] + EVENT_N_GUARD[i]);
        it can only fire while every one of them is positive.
        """, body)

    # ---- metadata ---------------------------------------------------------
    kind_of = {EventKind.REGION_NEG: "EVENT_REGION_NEG",
               EventKind.CROSS_POS_TO_NEG: "EVENT_CROSS_POS_TO_NEG",
               EventKind.CROSS_NEG_TO_POS: "EVENT_CROSS_NEG_TO_POS"}
    lit = lambda xs: ", ".join(xs)
    meta = f"""\
constexpr int N_PARAM = {len(spec.params)};
constexpr int N_EVENT = {len(spec.events)};
constexpr int N_GUARD = {len(guard_exprs)};

inline constexpr char const* PARAM_NAME[N_PARAM] = {{ {lit(f'"{p.name}"' for p in spec.params)} }};
inline constexpr double PARAM_DEFAULT[N_PARAM] = {{ {lit(repr(float(p.default)) for p in spec.params)} }};
inline constexpr char const* PARAM_DOC[N_PARAM] = {{ {lit(f'"{p.doc}"' for p in spec.params)} }};

inline constexpr char const* EVENT_NAME[N_EVENT] = {{ {lit(f'"{e.name}"' for e in spec.events)} }};
inline constexpr bh_event_kind EVENT_KIND[N_EVENT] = {{ {lit(kind_of[e.kind] for e in spec.events)} }};
inline constexpr int EVENT_N_GUARD[N_EVENT] = {{ {lit(str(len(e.guards)) for e in spec.events)} }};
inline constexpr int EVENT_GUARD_OFFSET[N_EVENT] = {{ {lit(str(o) for o in guard_off)} }};

//! Signature convention: +1 means (+,-,-,-), -1 means (-,+,+,+).
constexpr int SIGNATURE = {spec.signature};
"""
    for i, e in enumerate(spec.events):
        meta += f"constexpr int EV_{e.name} = {i};\n"

    stat_lines = "\n".join(f"//   {k:<14} {v:>4} CSE temporaries" for k, v in stats.items())
    header = f"""\
// ===========================================================================
//  GENERATED FILE - DO NOT EDIT.
//
//  Metric  : {spec.name}
//  Chart   : {spec.chart_doc}
//
//  Regenerate with:  make -C codegen
//  Source of truth:  codegen/metrics/{spec.name}.py
//
{textwrap.indent(provenance.strip(), '//  ') if provenance else '//'}
//
//  Code size:
{stat_lines}
// ===========================================================================
#pragma once

#include "bh_metric_api.h"

namespace {ns}
{{

{meta}
{f_metric}
{f_geo}
{f_obs}{f_ev}
{f_gd}
}} // namespace {ns}
"""
    return header, stats
