#!/usr/bin/env python3
"""Interactive web UI: drag the mouse to orbit the camera around the hole.

Shells out to the existing, unmodified ./main CLI + cli_imagemaker.py
pipeline for every frame; this file only tracks camera state and turns mouse
deltas into --theta0/--phi0 renders. While the mouse is held down we render
with relaxed --errormax/--de0 (fewer adaptive steps, lower accuracy, higher
fps); on mouseup we re-render once at normal accuracy.
"""

import math
import subprocess
import sys
import time
from pathlib import Path
from threading import Lock

from fastapi import FastAPI
from fastapi.responses import FileResponse, HTMLResponse
from pydantic import BaseModel

STATIC_DIR = Path("static")
IMAGE_PATH = Path("web_images") / "blackhole_cli.png"

SZELES = 640
MAGAS = 320

# Loose tolerances while dragging trade accuracy for fps; tighter ones on
# release give a crisp resting frame. Same knobs cli_wrapper_server.py
# already exposes on ./main, just two fixed presets instead of a UI control.
DRAG_ACCURACY = {"--errormax": "0.01", "--de0": "0.05"}
SETTLE_ACCURACY = {"--errormax": "0.001", "--de0": "0.01"}

POLE_EPS = 0.02  # keep theta0 off the axis; cuda_ray.h floors sin(theta) but stay clear of it

app = FastAPI()

state_lock = Lock()
state = {
    "r0": 1.0,
    "theta0": 1.57 + 0.06,
    "phi0": 0.0,
}


def render(accuracy: dict) -> float:
    args = ["./main", "--r0", str(state["r0"]), "--theta0", str(state["theta0"]),
            "--phi0", str(state["phi0"]), "--SZELES", str(SZELES), "--MAGAS", str(MAGAS)]
    for flag, value in accuracy.items():
        args += [flag, value]

    start = time.time()
    subprocess.run(args, check=True, capture_output=True)
    subprocess.run([sys.executable, "cli_imagemaker.py"], check=True, capture_output=True)
    return time.time() - start


class OrbitRequest(BaseModel):
    dtheta: float = 0.0
    dphi: float = 0.0
    dragging: bool = True


@app.get("/", response_class=HTMLResponse)
def index() -> str:
    return (STATIC_DIR / "index.html").read_text()


@app.post("/orbit")
def orbit(req: OrbitRequest):
    with state_lock:
        state["theta0"] = min(max(state["theta0"] + req.dtheta, POLE_EPS), math.pi - POLE_EPS)
        state["phi0"] = (state["phi0"] + req.dphi) % (2 * math.pi)
        elapsed = render(DRAG_ACCURACY if req.dragging else SETTLE_ACCURACY)
        return {
            "elapsed": elapsed,
            "theta0": state["theta0"],
            "phi0": state["phi0"],
            "image": f"/image?t={time.time()}",
        }


@app.get("/image")
def image():
    return FileResponse(IMAGE_PATH, headers={"Cache-Control": "no-store"})


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
