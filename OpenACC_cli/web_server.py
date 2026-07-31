#!/usr/bin/env python3
"""Interactive web UI: drag the mouse to orbit, scroll to zoom.

Shells out to the existing, unmodified ./main CLI + cli_imagemaker.py
pipeline for every frame; this file only tracks camera state and turns mouse
input into --r0/--theta0/--phi0 renders.

Three accuracy tiers trade render time for image quality:
  - DRAG:   while the mouse button/wheel is actively moving (fast, coarse)
  - SETTLE: one frame right after the camera stops (medium)
  - HQ:     kicked off in the background a moment after that, and swapped in
            once done, if the camera hasn't moved again in the meantime

Every rendered frame is cached to disk under cache/orbit/, keyed by a
quantized (r0, theta0, phi0, tier) so revisiting a viewpoint - including the
default start view on every page load - is instant.
"""

import hashlib
import math
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path
from threading import Lock

from fastapi import FastAPI
from fastapi.responses import FileResponse, HTMLResponse
from pydantic import BaseModel

STATIC_DIR = Path("static")
IMAGE_PATH = Path("web_images") / "blackhole_cli.png"
CACHE_DIR = Path("cache") / "orbit"
HQ_SCRATCH_DIR = Path("web_images") / "hq_scratch"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Resolved once at startup (before any cwd assumptions matter) so renders
# into a different workdir - the HQ scratch dir - can still find them.
MAIN_BIN = Path("main").resolve()
IMAGEMAKER = Path("cli_imagemaker.py").resolve()

SZELES = 640
MAGAS = 320

DRAG_ACCURACY = {"--errormax": "0.01", "--de0": "0.05"}
SETTLE_ACCURACY = {"--errormax": "0.001", "--de0": "0.01"}
HQ_ACCURACY = {"--errormax": "0.00005", "--de0": "0.0005"}

POLE_EPS = 0.02  # keep theta0 off the axis; cuda_ray.h floors sin(theta) but stay clear of it
MIN_R0 = 0.6     # stay outside the disk's outer edge (gyuru_sugar_nagy default 0.5)
MAX_R0 = 6.0

QUANT_ANGLE = 0.02
QUANT_R = 0.05

HQ_SETTLE_DEBOUNCE_S = 0.35  # let a quick flick-through cancel before burning a slow HQ render

app = FastAPI()

state_lock = Lock()
state = {
    "r0": 1.0,
    "theta0": 1.57 + 0.06,
    "phi0": 0.0,
}
generation = 0
hq_ready = {"generation": -1}


def cache_key(r0: float, theta0: float, phi0: float, tier: str) -> str:
    qr = round(r0 / QUANT_R) * QUANT_R
    qt = round(theta0 / QUANT_ANGLE) * QUANT_ANGLE
    qp = round(phi0 / QUANT_ANGLE) * QUANT_ANGLE
    raw = f"{tier}:{SZELES}x{MAGAS}:{qr:.3f}:{qt:.4f}:{qp:.4f}"
    return hashlib.sha256(raw.encode()).hexdigest()[:24]


def render_frame(r0: float, theta0: float, phi0: float, accuracy: dict, workdir: Path) -> Path:
    """Run ./main + cli_imagemaker.py for one frame into workdir, return the PNG path."""
    workdir.mkdir(parents=True, exist_ok=True)
    args = [str(MAIN_BIN), "--r0", str(r0), "--theta0", str(theta0), "--phi0", str(phi0),
            "--SZELES", str(SZELES), "--MAGAS", str(MAGAS)]
    for flag, value in accuracy.items():
        args += [flag, value]
    subprocess.run(args, check=True, capture_output=True, cwd=workdir)
    subprocess.run([sys.executable, str(IMAGEMAKER)], check=True, capture_output=True, cwd=workdir)
    return workdir / "web_images" / "blackhole_cli.png"


def render_cached(r0: float, theta0: float, phi0: float, tier: str, accuracy: dict, workdir: Path) -> tuple[float, Path]:
    """Render (or reuse a cached render of) one frame. Returns (elapsed_seconds, png_path)."""
    key = cache_key(r0, theta0, phi0, tier)
    cached = CACHE_DIR / f"{key}.png"
    if cached.exists():
        return 0.0, cached

    start = time.time()
    img_path = render_frame(r0, theta0, phi0, accuracy, workdir)
    elapsed = time.time() - start
    shutil.copy(img_path, cached)
    return elapsed, cached


def hq_worker(gen: int, r0: float, theta0: float, phi0: float) -> None:
    time.sleep(HQ_SETTLE_DEBOUNCE_S)
    with state_lock:
        if generation != gen:
            return  # camera moved again before we even started; skip the slow render entirely

    try:
        _, cached_path = render_cached(r0, theta0, phi0, "hq", HQ_ACCURACY, HQ_SCRATCH_DIR)
    except subprocess.CalledProcessError:
        return

    with state_lock:
        if generation == gen:  # still the same viewpoint -> promote to the live image
            shutil.copy(cached_path, IMAGE_PATH)
            hq_ready["generation"] = gen


class OrbitRequest(BaseModel):
    dtheta: float = 0.0
    dphi: float = 0.0
    dzoom: float = 1.0  # multiplicative factor applied to r0
    dragging: bool = True


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return (STATIC_DIR / "index.html").read_text()


@app.post("/orbit")
def orbit(req: OrbitRequest):
    global generation
    with state_lock:
        state["theta0"] = min(max(state["theta0"] + req.dtheta, POLE_EPS), math.pi - POLE_EPS)
        state["phi0"] = (state["phi0"] + req.dphi) % (2 * math.pi)
        zoom = min(max(req.dzoom, 0.5), 2.0)
        state["r0"] = min(max(state["r0"] * zoom, MIN_R0), MAX_R0)
        generation += 1
        gen = generation

        tier = "drag" if req.dragging else "settle"
        accuracy = DRAG_ACCURACY if req.dragging else SETTLE_ACCURACY
        elapsed, cached_path = render_cached(state["r0"], state["theta0"], state["phi0"], tier, accuracy, Path("."))
        if cached_path != IMAGE_PATH:
            shutil.copy(cached_path, IMAGE_PATH)

        if not req.dragging:
            threading.Thread(target=hq_worker, args=(gen, state["r0"], state["theta0"], state["phi0"]), daemon=True).start()

        return {
            "elapsed": elapsed,
            "r0": state["r0"],
            "theta0": state["theta0"],
            "phi0": state["phi0"],
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
