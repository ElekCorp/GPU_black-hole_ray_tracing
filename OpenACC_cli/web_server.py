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

HQ_SETTLE_DEBOUNCE_S = 0.35  # let a quick flick-through cancel before burning a slow HQ render

app = FastAPI()

state_lock = Lock()
state = {
    "r0": 1.0,
    "theta0": 1.57 + 0.06,
    "phi0": 0.0,
    "a": 0.0,
    "Q": 0.0,
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


def cache_key(r0: float, theta0: float, phi0: float, a: float, Q: float, tier: str, resolution: tuple[int, int]) -> str:
    qr = round(r0 / QUANT_R) * QUANT_R
    qt = round(theta0 / QUANT_ANGLE) * QUANT_ANGLE
    qp = round(phi0 / QUANT_ANGLE) * QUANT_ANGLE
    qa = round(a / QUANT_AQ) * QUANT_AQ
    qq = round(Q / QUANT_AQ) * QUANT_AQ
    raw = f"{tier}:{resolution[0]}x{resolution[1]}:{qr:.3f}:{qt:.4f}:{qp:.4f}:{qa:.4f}:{qq:.4f}"
    return hashlib.sha256(raw.encode()).hexdigest()[:24]


def render_frame(worker: "render_worker.RenderWorker", r0: float, theta0: float, phi0: float, a: float, Q: float,
                  accuracy: dict, resolution: tuple[int, int]) -> Image.Image:
    szeles, magas = resolution
    buffer = worker.render(r0, theta0, phi0, a, Q, RS, accuracy["errormax"], accuracy["de0"], szeles, magas)
    hits, redshift, phi = cli_imagemaker.split_channels(szeles, magas, buffer.reshape(-1))
    return cli_imagemaker.cinematic_image(hits, redshift, phi, seed=23071969, peak_temperature=12_000.0)


def render_cached(worker: "render_worker.RenderWorker", r0: float, theta0: float, phi0: float, a: float, Q: float,
                   tier: str, accuracy: dict, resolution: tuple[int, int]) -> tuple[float, Path]:
    """Render (or reuse a cached render of) one frame. Returns (elapsed_seconds, png_path)."""
    key = cache_key(r0, theta0, phi0, a, Q, tier, resolution)
    cached = CACHE_DIR / f"{key}.png"
    if cached.exists():
        return 0.0, cached

    start = time.time()
    image = render_frame(worker, r0, theta0, phi0, a, Q, accuracy, resolution)
    elapsed = time.time() - start
    image.save(cached)
    return elapsed, cached


def hq_background_task(gen: int, r0: float, theta0: float, phi0: float, a: float, Q: float) -> None:
    time.sleep(HQ_SETTLE_DEBOUNCE_S)
    with state_lock:
        if generation != gen:
            return  # camera moved again before we even started; skip the slow render entirely

    try:
        _, cached_path = render_cached(render_worker.HQ_WORKER, r0, theta0, phi0, a, Q, "hq", HQ_ACCURACY, HQ_RESOLUTION)
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

        if req.a is not None or req.Q is not None:
            new_a = req.a if req.a is not None else state["a"]
            new_q = req.Q if req.Q is not None else state["Q"]
            state["a"], state["Q"] = clamp_spin_charge(new_a, new_q)

        min_r0 = min_r0_for(state["a"], state["Q"])
        zoom = min(max(req.dzoom, 0.5), 2.0)
        state["r0"] = min(max(state["r0"] * zoom, min_r0), MAX_R0)
        generation += 1
        gen = generation

        tier = "drag" if req.dragging else "settle"
        accuracy = DRAG_ACCURACY if req.dragging else SETTLE_ACCURACY
        resolution = DRAG_RESOLUTION if req.dragging else SETTLE_RESOLUTION
        elapsed, cached_path = render_cached(
            render_worker.INTERACTIVE_WORKER, state["r0"], state["theta0"], state["phi0"], state["a"], state["Q"],
            tier, accuracy, resolution,
        )
        shutil.copy(cached_path, IMAGE_PATH)

        if not req.dragging:
            threading.Thread(
                target=hq_background_task,
                args=(gen, state["r0"], state["theta0"], state["phi0"], state["a"], state["Q"]),
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
