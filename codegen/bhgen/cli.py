"""Driver: metric module in, verified C++ header out.

Nothing is written unless every proof obligation that is supposed to hold does
hold.  The report is printed and also embedded in the generated header, so a
checked-in header carries the provenance of its own verification.
"""

from __future__ import annotations

import argparse
import importlib
import json
import pathlib
import sys
import time

import sympy as sp

from . import emit, geodesic, tetrad
from .cache import cached
from .normalize import normalize
from .spec import DIM
from .verify import interval as IV, roundoff as RO, symbolic as SYM


def _t():
    return time.time()


def run(spec, outdir: pathlib.Path, do_roundoff=True, roundoff_box=None,
        interval_box=None, constraints=(), verbose=True):
    log = (lambda *a: print(*a, flush=True)) if verbose else (lambda *a: None)
    report = {"metric": spec.name, "obligations": []}

    def record(name, status, detail=""):
        report["obligations"].append({"name": name, "status": status, "detail": detail})
        mark = {"PROVED": "proved", "TIGHT": "proved", "PROVED_MODULO_BOUNDARY": "proved*",
                "PROBABLE": "probable",
                "LOOSE": "weak", "UNKNOWN": "unknown", "NO_BOUND": "unknown",
                "REFUTED": "REFUTED", "ERROR": "error"}.get(status, status)
        log(f"    [{mark:>8}] {name}{('  -- ' + detail) if detail else ''}")

    v = geodesic.velocity_symbols()

    log(f"\n=== {spec.name} ===\n  {spec.chart_doc}")
    log("  deriving geodesic right-hand side ...")
    t0 = _t()
    acc, hit = cached(spec, "acc",
                      lambda: [normalize(e) for e in geodesic.acceleration_from_christoffel(spec, v)])
    log(f"    {'from cache' if hit else 'done'} in {_t()-t0:.1f}s, "
        f"ops = {[sp.count_ops(e) for e in acc]}")

    log("  symbolic obligations:")
    t0 = _t()
    ok, how, bad = SYM.is_zero(SYM.check_cross_derivation(spec, acc=acc, v=v))
    record("geodesic RHS: Christoffel formula == Euler-Lagrange",
           ("PROBABLE" if "Schwartz" in how else "PROVED") if ok else "REFUTED",
           f"{how}, {_t()-t0:.1f}s" if ok else f"residual {bad}")

    t0 = _t()
    frame, _ = cached(spec, "frame", lambda: tetrad.spec_frame(spec))
    ok, how, bad = SYM.is_zero(SYM.check_tetrad(spec, frame))
    record("camera tetrad: e^T g e == eta",
           ("PROBABLE" if "Schwartz" in how else "PROVED") if ok else "REFUTED",
           f"{how}, {_t()-t0:.1f}s" if ok else f"residual {bad}")

    sub = emit._subs_map(spec, with_velocity=True)
    ex = [e.xreplace(sub) for e in acc]
    reps, red = sp.cse(ex, symbols=sp.numbered_symbols("t"), optimizations="basic")
    t0 = _t()
    ok, how, bad = SYM.is_zero(SYM.check_cse_equivalence(ex, reps, red))
    record(f"CSE ({len(reps)} temporaries) preserves the expression",
           ("PROBABLE" if "Schwartz" in how else "PROVED") if ok else "REFUTED",
           f"{how}, {_t()-t0:.1f}s" if ok else f"residual {bad}")

    if interval_box:
        log("  domain-safety obligations (interval arithmetic):")
        obls = IV.collect_obligations(acc)
        for o in obls:
            t0 = _t()
            r = IV.prove(o, interval_box, constraints=constraints, max_boxes=40000)
            detail = f"{r.boxes_examined} boxes, {_t()-t0:.1f}s"
            if r.witness:
                detail += "  witness " + " ".join(f"{k}=[{a:.4g},{b:.4g}]" for k, (a, b) in r.witness.items())
            elif r.boundary_volume:
                detail += f", unproved boundary shell {r.boundary_volume:.3g}"
            record(f"{o.kind}: {o.origin[:70]}", r.status, detail)

    if do_roundoff and roundoff_box:
        log("  floating-point round-off obligations (Gappa):")
        for prec in ("double", "float"):
            for i in range(DIM):
                rng = RO._range_over_box(red[i], roundoff_box, reps)
                mag = max(abs(rng[0]), abs(rng[1])) if rng else None
                try:
                    script, _ = RO.build_script(reps, red, [f"acc{i}"], roundoff_box,
                                                precision=prec, goal_index=i)
                    r = RO.classify(RO.run(script, timeout=300), mag, prec)
                    if r.status == "LOOSE":
                        sym0 = list(roundoff_box)[0]
                        script, _ = RO.build_script(reps, red, [f"acc{i}"], roundoff_box,
                                                    precision=prec, goal_index=i,
                                                    subdivide={sym0: 40})
                        r = RO.classify(RO.run(script, timeout=600), mag, prec)
                    detail = (f"|fl-exact| <= {r.abs_error:.3e}, relative {r.relative:.2e}"
                              if r.abs_error is not None else r.raw.strip().splitlines()[-1][:70])
                except RO.UnboundedLibm as e:
                    r, detail = RO.RoundoffResult(f"acc{i}", "ERROR", None, str(e)), str(e)[:70]
                record(f"round-off {prec} acc[{i}]", r.status, detail)

    refuted = [o for o in report["obligations"] if o["status"] == "REFUTED"]
    if refuted:
        log(f"\n  REFUTED: {len(refuted)} obligation(s) failed; nothing written.")
        return report, None

    prov = "Verification at generation time:\n" + "\n".join(
        f"  [{o['status']:<22}] {o['name']}" for o in report["obligations"])
    header, stats = emit.generate_header(spec, acc, frame, provenance=prov)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"bh_{spec.name}_metric.h"
    path.write_text(header)
    log(f"\n  wrote {path}  ({len(header.splitlines())} lines, {stats})")
    return report, path


def main(argv=None):
    ap = argparse.ArgumentParser(description="Generate and verify a metric header.")
    ap.add_argument("metrics", nargs="+", help="module names under codegen/metrics/")
    ap.add_argument("-o", "--outdir", default="../OpenACC_cli/generated")
    ap.add_argument("--no-roundoff", action="store_true")
    ap.add_argument("--report", default=None, help="write the JSON report here")
    args = ap.parse_args(argv)

    reports = []
    for name in args.metrics:
        mod = importlib.import_module(f"metrics.{name}")
        spec = mod.SPEC
        box = getattr(mod, "VERIFY_BOX", None)
        cons = getattr(mod, "VERIFY_CONSTRAINTS", ())
        rbox = getattr(mod, "ROUNDOFF_BOX", None)
        rep, path = run(spec, pathlib.Path(args.outdir), do_roundoff=not args.no_roundoff,
                        roundoff_box=rbox, interval_box=box, constraints=cons)
        reports.append(rep)
        if path is None:
            return 1
    if args.report:
        pathlib.Path(args.report).write_text(json.dumps(reports, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
