#!/usr/bin/env python3
"""Bisect the shadow edge at cx=0.5 using the just-built libblackhole.so, so
run_bench.sh's "zoomed" region sits on the photon-ring edge (the critical-
impact-parameter ring) instead of blindly zooming into the shadow's (empty)
centre.

Run from the repo root (OpenACC_cli/), after `make` has produced
libblackhole.so for whichever scene parameters you pass below. Prints the
resulting cy to stdout; paste it into run_bench.sh's EDGE_CY if you change
a/Q/rs/r0/theta0 from the defaults baked in here.
"""
import ctypes
import sys
from pathlib import Path

LIB = str(Path(__file__).resolve().parent.parent / "libblackhole.so")

lib = ctypes.CDLL(LIB)
scene_args = [ctypes.c_double] * 14 + [ctypes.c_uint64, ctypes.c_uint64]
lib.render_frame_roi_f64.argtypes = scene_args + [ctypes.POINTER(ctypes.c_double)]
lib.render_frame_roi_f64.restype = ctypes.c_int
lib.render_frame_channels.restype = ctypes.c_int
channels = lib.render_frame_channels()

# Matches run_bench.sh's scene: a=0.02 (nonzero spin, so christoffel() always
# takes the christoffel_general branch), Q=0, rs=0.05/r0=1.0/theta0=1.63
# (cli_parser.h defaults).
R0, THETA0, PHI0 = 1.0, 1.57 + 0.06, 0.0
A, Q, RS = 0.02, 0.0, 0.05
ERRORMAX, DE0 = 0.0001, 0.01
OMEGA = (0.0, 3.14159265358979, 0.0)


def captured_at(cx, cy, w=1e-12):
    buf = (ctypes.c_double * (channels * 2 * 1))()
    rc = lib.render_frame_roi_f64(R0, THETA0, PHI0, A, Q, RS, ERRORMAX, DE0,
                                   OMEGA[0], OMEGA[1], OMEGA[2],
                                   cx, cy, w, 2, 1, buf)
    assert rc == 0, rc
    return buf[0] == 0.0


def find_edge(samples=64, bisection_steps=60):
    rows = [(i + 0.5) / samples for i in range(samples)]
    flags = [captured_at(0.5, cy) for cy in rows]
    print("scan:", "".join("#" if f else "." for f in flags), file=sys.stderr)

    bracket = None
    for i in range(samples - 1):
        if flags[i] != flags[i + 1]:
            bracket = (rows[i], rows[i + 1], flags[i + 1])
            break
    if bracket is None:
        raise SystemExit("no shadow edge found in scan - widen the scan or check scene params")

    upper, lower, captured_below = bracket
    sky_cy, shadow_cy = (upper, lower) if captured_below else (lower, upper)
    for _ in range(bisection_steps):
        mid = (sky_cy + shadow_cy) / 2.0
        if captured_at(0.5, mid):
            shadow_cy = mid
        else:
            sky_cy = mid
    return (sky_cy + shadow_cy) / 2.0


if __name__ == "__main__":
    print(f"{find_edge()!r}")
