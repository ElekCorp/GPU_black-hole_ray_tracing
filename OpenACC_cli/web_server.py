#!/usr/bin/env python3
"""Interactive web UI: drag the mouse to orbit, scroll to zoom.

Renders through libblackhole.so (see ffi_bridge.cpp / render_worker.py) - two
persistent worker processes that dlopen the ray tracer once and keep their
CUDA context warm across many frames, instead of spawning ./main fresh (with
a full CUDA context init) on every single request. This file tracks camera
state and turns mouse/slider input into renders through those workers.

Three accuracy tiers trade render time for image quality:
  - DRAG:   while the mouse button/wheel/slider is actively moving (fast, coarse)
  - SETTLE: one frame right after the camera stops (medium)
  - HQ:     kicked off in the background a moment after that, and swapped in
            once done, if the camera hasn't moved again in the meantime

HQ renders through a second worker (its own process/CUDA context) so a slow
background HQ pass never blocks interactive frames.

Every rendered frame is cached to disk under cache/orbit/, keyed by a
quantized (r0, theta0, phi0, a, Q, tier) so revisiting a viewpoint -
including the default start view on every page load - is instant.
"""

import hashlib
import math
import shutil
import threading
import time
from pathlib import Path
from threading import Lock
from typing import Optional

import numpy as np
from fastapi import FastAPI
from fastapi.responses import FileResponse, HTMLResponse
from PIL import Image
from pydantic import BaseModel

import cli_imagemaker
import render_worker

STATIC_DIR = Path("static")
IMAGE_PATH = Path("web_images") / "blackhole_cli.png"
CACHE_DIR = Path("cache") / "orbit"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
IMAGE_PATH.parent.mkdir(parents=True, exist_ok=True)

# All three keep the 2:1 aspect ratio implied by main.cpp's default
# kepernyoSZELES/kepernyoMAGAS (10240:5120); ffi_bridge.cpp rejects a mismatch.
DRAG_RESOLUTION = (160, 80)      # 1/16 the pixels of settle - the fps knob
SETTLE_RESOLUTION = (640, 320)   # original baseline resolution
HQ_RESOLUTION = (1280, 640)      # 4x the pixels of settle - background-only, so slow is fine

DRAG_ACCURACY = {"errormax": 0.01, "de0": 0.05}
SETTLE_ACCURACY = {"errormax": 0.001, "de0": 0.01}
HQ_ACCURACY = {"errormax": 0.00005, "de0": 0.0005}

POLE_EPS = 0.02  # keep theta0 off the axis; cuda_ray.h floors sin(theta) but stay clear of it
MAX_R0 = 6.0
MIN_R0_PAD = 1.3   # multiplier over the horizon radius the camera must stay outside of
MIN_R0_FLOOR = 0.03  # absolute floor regardless of horizon size

RS = 0.05  # matches cli_parser.h's default; not exposed as adjustable here
SPIN_CHARGE_MARGIN = 0.98  # stay this fraction below extremal (a^2+Q^2 = (rs/2)^2) to keep a horizon

QUANT_ANGLE = 0.02
QUANT_R = 0.05
QUANT_AQ = 0.001
QUANT_OMEGA = 0.01

# ijk_to_vec_zoom's camera-orientation axes, in the local (pre-rotation) screen
# frame ijk_to_vec_mink_zoom builds: +x is forward (along kepernyo_tav), +y is
# the screen's vertical axis, +z is the screen's horizontal axis.
LOCAL_FORWARD = np.array([1.0, 0.0, 0.0])
LOCAL_UP = np.array([0.0, 1.0, 0.0])
LOCAL_RIGHT = np.array([0.0, 0.0, 1.0])

HQ_SETTLE_DEBOUNCE_S = 0.35  # let a quick flick-through cancel before burning a slow HQ render


# Orientation is tracked as a quaternion (w, x, y, z), not a matrix/raw
# axis-angle vector. The default orientation is a full 180-degree rotation
# (see BASE_ORIENTATION below), and axis-angle extraction from a matrix is
# numerically unstable there (the standard R_32-R_23-style formula divides by
# sin(angle), which is ~0 right at 180 degrees - confirmed directly: two
# matrices that differed by ~1e-16 extracted wildly different axes). A
# quaternion's axis-angle extraction instead divides by sin(angle/2), which
# is perfectly well-conditioned at 180 degrees (only degenerates near 0
# degrees, i.e. "no rotation", which is comparatively harmless).
def axis_angle_to_quat(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    s = math.sin(angle / 2)
    return np.array([math.cos(angle / 2), axis[0] * s, axis[1] * s, axis[2] * s])


def quat_multiply(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ])


def quat_rotate_vector(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    qv = np.array([0.0, v[0], v[1], v[2]])
    q_conj = q * np.array([1.0, -1.0, -1.0, -1.0])
    return quat_multiply(quat_multiply(q, qv), q_conj)[1:]


def quat_to_axis_angle(q: np.ndarray) -> tuple[np.ndarray, float]:
    w = min(1.0, max(-1.0, q[0]))
    angle = 2 * math.acos(w)
    s = math.sqrt(max(0.0, 1 - w * w))  # sin(angle/2)
    if s < 1e-8:
        return np.array([0.0, 1.0, 0.0]), 0.0
    return q[1:] / s, angle


def apply_look_delta(orientation: np.ndarray, dpan: float, dtilt: float, droll: float) -> np.ndarray:
    """Pan/tilt/roll the camera's FACING direction, independent of its position on the
    orbit sphere - i.e. rotate ijk_to_vec_zoom's Omega, not r0/theta0/phi0. Axes are
    taken from the camera's CURRENT orientation (not fixed world axes), so this composes
    the way first-person look controls are expected to: panning always turns relative to
    whatever the camera is presently facing."""
    if dpan == 0.0 and dtilt == 0.0 and droll == 0.0:
        return orientation
    forward = quat_rotate_vector(orientation, LOCAL_FORWARD)
    up = quat_rotate_vector(orientation, LOCAL_UP)
    right = quat_rotate_vector(orientation, LOCAL_RIGHT)
    delta = quat_multiply(
        quat_multiply(axis_angle_to_quat(up, dpan), axis_angle_to_quat(right, dtilt)),
        axis_angle_to_quat(forward, droll),
    )
    result = quat_multiply(delta, orientation)
    return result / np.linalg.norm(result)  # renormalize against drift from repeated compositions


def current_omega(orientation: np.ndarray) -> tuple[float, float, float]:
    axis, angle = quat_to_axis_angle(orientation)
    vector = axis * angle
    return float(vector[0]), float(vector[1]), float(vector[2])


# main.cpp/ffi_bridge.cpp's previous hardcoded default: a fixed 180-degree
# flip around the local up axis, which is what makes the camera face the hole
# by default. User pan/tilt/roll composes on top of this, not instead of it.
BASE_ORIENTATION = axis_angle_to_quat(LOCAL_UP, math.pi)

app = FastAPI()

state_lock = Lock()
state = {
    "r0": 1.0,
    "theta0": 1.57 + 0.06,
    "phi0": 0.0,
    "a": 0.0,
    "Q": 0.0,
    "orientation": BASE_ORIENTATION.copy(),
}
generation = 0
hq_ready = {"generation": -1}


def horizon_radius(a: float, Q: float) -> float:
    """Outer Kerr-Newman horizon r_+ for the fixed RS; see the sugar_be calculation in
    ray_step's collision setup, which uses the same discriminant."""
    disc = RS * RS - 4 * (a * a + Q * Q)
    if disc <= 0:
        return RS / 2  # extremal edge; clamp_spin_charge keeps us just short of this
    return (RS + math.sqrt(disc)) / 2


def min_r0_for(a: float, Q: float) -> float:
    return max(horizon_radius(a, Q) * MIN_R0_PAD, MIN_R0_FLOOR)


def clamp_spin_charge(a: float, Q: float) -> tuple[float, float]:
    """Keep a^2+Q^2 safely below the extremal (rs/2)^2 bound so a horizon always exists,
    scaling both down proportionally rather than favoring whichever was set last."""
    limit = (RS / 2) * SPIN_CHARGE_MARGIN
    magnitude = math.hypot(a, Q)
    if magnitude > limit and magnitude > 0:
        scale = limit / magnitude
        a *= scale
        Q *= scale
    return a, Q


def cache_key(r0: float, theta0: float, phi0: float, a: float, Q: float, omega: tuple, tier: str, resolution: tuple[int, int]) -> str:
    qr = round(r0 / QUANT_R) * QUANT_R
    qt = round(theta0 / QUANT_ANGLE) * QUANT_ANGLE
    qp = round(phi0 / QUANT_ANGLE) * QUANT_ANGLE
    qa = round(a / QUANT_AQ) * QUANT_AQ
    qq = round(Q / QUANT_AQ) * QUANT_AQ
    qo = tuple(round(v / QUANT_OMEGA) * QUANT_OMEGA for v in omega)
    raw = f"{tier}:{resolution[0]}x{resolution[1]}:{qr:.3f}:{qt:.4f}:{qp:.4f}:{qa:.4f}:{qq:.4f}:{qo[0]:.3f}:{qo[1]:.3f}:{qo[2]:.3f}"
    return hashlib.sha256(raw.encode()).hexdigest()[:24]


def render_frame(worker: "render_worker.RenderWorker", r0: float, theta0: float, phi0: float, a: float, Q: float,
                  omega: tuple, accuracy: dict, resolution: tuple[int, int]) -> Image.Image:
    szeles, magas = resolution
    buffer = worker.render(r0, theta0, phi0, a, Q, RS, accuracy["errormax"], accuracy["de0"], omega, szeles, magas)
    hits, redshift, phi = cli_imagemaker.split_channels(szeles, magas, buffer.reshape(-1))
    return cli_imagemaker.cinematic_image(hits, redshift, phi, seed=23071969, peak_temperature=12_000.0)


def render_cached(worker: "render_worker.RenderWorker", r0: float, theta0: float, phi0: float, a: float, Q: float,
                   omega: tuple, tier: str, accuracy: dict, resolution: tuple[int, int]) -> tuple[float, Path]:
    """Render (or reuse a cached render of) one frame. Returns (elapsed_seconds, png_path)."""
    key = cache_key(r0, theta0, phi0, a, Q, omega, tier, resolution)
    cached = CACHE_DIR / f"{key}.png"
    if cached.exists():
        return 0.0, cached

    start = time.time()
    image = render_frame(worker, r0, theta0, phi0, a, Q, omega, accuracy, resolution)
    elapsed = time.time() - start
    image.save(cached)
    return elapsed, cached


def hq_background_task(gen: int, r0: float, theta0: float, phi0: float, a: float, Q: float, omega: tuple) -> None:
    time.sleep(HQ_SETTLE_DEBOUNCE_S)
    with state_lock:
        if generation != gen:
            return  # camera moved again before we even started; skip the slow render entirely

    try:
        _, cached_path = render_cached(render_worker.HQ_WORKER, r0, theta0, phi0, a, Q, omega, "hq", HQ_ACCURACY, HQ_RESOLUTION)
    except RuntimeError:
        return

    with state_lock:
        if generation == gen:  # still the same viewpoint -> promote to the live image
            shutil.copy(cached_path, IMAGE_PATH)
            hq_ready["generation"] = gen


class OrbitRequest(BaseModel):
    dtheta: float = 0.0
    dphi: float = 0.0
    dzoom: float = 1.0  # multiplicative factor applied to r0
    dpan: float = 0.0    # camera orientation deltas (independent of r0/theta0/phi0):
    dtilt: float = 0.0   # pan turns around the camera's current up axis,
    droll: float = 0.0   # tilt around its current right axis, roll around its current forward axis
    a: Optional[float] = None  # absolute spin; omitted/None means "leave unchanged"
    Q: Optional[float] = None  # absolute charge; omitted/None means "leave unchanged"
    dragging: bool = True


@app.on_event("shutdown")
def _shutdown_workers() -> None:
    render_worker.shutdown_all()


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return (STATIC_DIR / "index.html").read_text()


@app.post("/orbit")
def orbit(req: OrbitRequest):
    global generation
    with state_lock:
        state["theta0"] = min(max(state["theta0"] + req.dtheta, POLE_EPS), math.pi - POLE_EPS)
        state["phi0"] = (state["phi0"] + req.dphi) % (2 * math.pi)
        state["orientation"] = apply_look_delta(state["orientation"], req.dpan, req.dtilt, req.droll)

        if req.a is not None or req.Q is not None:
            new_a = req.a if req.a is not None else state["a"]
            new_q = req.Q if req.Q is not None else state["Q"]
            state["a"], state["Q"] = clamp_spin_charge(new_a, new_q)

        min_r0 = min_r0_for(state["a"], state["Q"])
        zoom = min(max(req.dzoom, 0.5), 2.0)
        state["r0"] = min(max(state["r0"] * zoom, min_r0), MAX_R0)
        generation += 1
        gen = generation
        omega = current_omega(state["orientation"])

        tier = "drag" if req.dragging else "settle"
        accuracy = DRAG_ACCURACY if req.dragging else SETTLE_ACCURACY
        resolution = DRAG_RESOLUTION if req.dragging else SETTLE_RESOLUTION
        elapsed, cached_path = render_cached(
            render_worker.INTERACTIVE_WORKER, state["r0"], state["theta0"], state["phi0"], state["a"], state["Q"],
            omega, tier, accuracy, resolution,
        )
        shutil.copy(cached_path, IMAGE_PATH)

        if not req.dragging:
            threading.Thread(
                target=hq_background_task,
                args=(gen, state["r0"], state["theta0"], state["phi0"], state["a"], state["Q"], omega),
                daemon=True,
            ).start()

        return {
            "elapsed": elapsed,
            "r0": state["r0"],
            "theta0": state["theta0"],
            "phi0": state["phi0"],
            "a": state["a"],
            "Q": state["Q"],
            "min_r0": min_r0,
            "generation": gen,
            "image": f"/image?t={time.time()}",
        }


@app.get("/status")
def status():
    with state_lock:
        return {"generation": generation, "hq_ready": hq_ready["generation"] == generation}


@app.get("/image")
def image():
    return FileResponse(IMAGE_PATH, headers={"Cache-Control": "no-store"})


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
