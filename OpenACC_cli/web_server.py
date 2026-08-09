#!/usr/bin/env python3
"""Interactive web UI: drag the mouse to orbit, scroll to zoom.

Two different things are called "zoom" here, and they are not the same:
  - TRUE ZOOM (the scroll wheel, the box-zoom, the roi_* fields): narrows the
    camera's field of view onto a sub-rectangle of the frame and retraces the
    full pixel budget across it, so magnified detail is really traced rather
    than interpolated. The camera does not move. See render_frame_roi_f32.
  - DOLLY (dzoom): moves the camera in or out along r0. That changes the
    physics of what is visible - the lensing, the redshift, the disk's
    apparent geometry - not just how much of the sky the frame covers.

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
import os
import shutil
import tempfile
import threading
import time
from collections import OrderedDict
from pathlib import Path
from threading import Lock
from typing import Optional

import numpy as np
from fastapi import FastAPI
from fastapi.responses import HTMLResponse, Response
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

RENDER_SEED = 23071969  # fixed star-field/grain seed, shared by every render/recolor call
PEAK_TEMPERATURE = 12_000.0  # peak local disk temperature in kelvin

POLE_EPS = 0.02  # keep theta0 off the axis; cuda_ray.h floors sin(theta) but stay clear of it
MAX_R0 = 6.0
MIN_R0_PAD = 1.3   # multiplier over the horizon radius the camera must stay outside of
MIN_R0_FLOOR = 0.03  # absolute floor regardless of horizon size

RS = 0.05  # matches cli_parser.h's default; the rs slider's starting point
SPIN_CHARGE_MARGIN = 0.98  # stay this fraction below extremal (a^2+Q^2 = (rs/2)^2) to keep a horizon

# Bounds for the rs slider. rs is a pure length scale - rescaling r, a, Q, rs
# and t together leaves null geodesics unchanged as curves - so on its own it
# would do nothing. What makes it a real control here is that the geometry it is
# measured against does NOT scale with it: ffi_bridge.cpp hardcodes the disk
# annulus at 0.1 .. 0.5 and the outer integration sphere at 1.01, and the camera
# sits at r0 <= ~1. Moving rs therefore resizes the hole against a fixed disk and
# a fixed camera, which is exactly the knob it looks like.
#
# Both ends are set by that fixed geometry, measured at r0 = 1, theta0 = 1.63:
#
#   rs     shadow   note
#   0.005    0%     critical impact parameter 0.013 against a disk starting at
#   0.010    0%     0.1 - the hole is a dot, there is nothing to look at
#   0.050    4%     the current default
#   0.120   22%     the classic view: disk lensed up and over a large shadow
#   0.200   58%     hole dominates the frame, disk squeezed to the edges
#   0.350   96%     almost nothing but shadow left
#   0.500  100%     b_crit = 1.299*rs now exceeds r0 itself, so the camera is
#                   inside the capture region and every ray falls in: solid black
#
# rs = 0 is allowed as a distinct case rather than as "very small": the metric
# collapses to Minkowski, so the same disk and star field are traced with no
# lensing, no shadow and no redshift at all. That is the reference image every
# lensed frame should be read against, which is worth a slider stop of its own.
# Nothing downstream needs a special case for it. christoffel_static (a = 0
# takes that branch) reduces to the flat-space spherical-coordinate connection -
# still nonzero coefficients, but zero curvature, and its denominators are
# powers of r and of r^2-rs*r+Q^2 = r^2, neither of which vanishes outside the
# origin. ray_step leaves its capture radius at 0 because the discriminant
# rs^2-4(a^2+Q^2) is not positive, so no ray is captured and there is no shadow.
# disk_redshift returns an unshifted 1 once rs*r/2 - Q^2 goes non-positive. And
# the a and Q sliders collapse along with it: clamp_spin_charge forces both to 0
# at rs = 0, since a hole with no mass cannot carry spin or charge.
MIN_RS = 0.0
MAX_RS = 0.20

# What actually limits the true zoom is the precision the geodesics are traced
# in, so the renderer switches paths rather than stopping at the float ceiling.
# Both numbers below are measured, not assumed: rendering a smooth patch of
# lensed sky at rising magnification and asking (a) what fraction of adjacent
# pixels still resolve distinct rays, and (b) how much of the row-to-row
# variation is integration hash rather than gradient.
#
#   f32: clean to ~512x; by 1024x the gradient is half noise, by 4096x one pixel
#        in ten has collapsed onto its neighbour's geodesic, and past ~1e6x the
#        frame is flat blocks. Tightening errormax does not move this - what
#        runs out is the float representation of the ray state itself.
#   f64: distinct 1.000 and no measurable noise all the way to 2^32, first
#        showing any at 2^36.
#
# So: trace in f32 while the zoom is shallow enough for it (it is much cheaper
# on GPUs with a reduced fp64 rate), and switch to f64 beyond that.
FLOAT_ZOOM_LIMIT = 512.0
MAX_ROI_ZOOM = float(2**32)

# "auto" applies FLOAT_ZOOM_LIMIT above; "f32"/"f64" pin the trace regardless of
# zoom. Pinning is worth having in both directions and not just as an escape
# hatch: forcing f64 at low zoom is how the float path's error can be measured
# against a reference, and forcing f32 at high zoom is how its breakdown can be
# seen rather than taken on trust.
PRECISION_MODES = ("auto", "f32", "f64")

# The integrator is tightened as the field of view narrows, so that what
# neighbouring rays differ by stays geometry rather than integration error.
# Capped only to bound the worst case; in practice this costs very little,
# because de0 caps the step size anyway and the tolerance rarely binds - and
# dopri54_step floors the tolerance at 1e-8 regardless of what is passed.
# de0 itself is deliberately not touched: it is the step *ceiling*, and
# shrinking it forces small steps everywhere rather than only where the
# geometry needs them.
ZOOM_ACCURACY_CAP = 1024.0

# Bumped whenever the renderer's output or the shading model changes, so a
# cache full of frames drawn by the previous version is bypassed rather than
# served alongside new ones. Old entries under cache/orbit/ are then dead
# weight and can simply be deleted.
SHADING_VERSION = 5

QUANT_ANGLE = 0.02
QUANT_R = 0.05
QUANT_AQ = 0.001
QUANT_RS = 0.001  # finer than the rs slider's own step, so no two slider stops share a key
QUANT_OMEGA = 0.01
# The zoom level is quantized geometrically (in log2 steps) and the window's
# centre in units of the window's own width - so the cache's tolerance shrinks
# with the window instead of staying a fixed slice of the unzoomed frame, which
# at high magnification would collapse wildly different views onto one key.
QUANT_ROI_ZOOM_STEP = 1.0 / 16.0
QUANT_ROI_CENTRE = 0.002

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
    "rs": RS,
    "orientation": BASE_ORIENTATION.copy(),
    # The true-zoom window, as a fraction of the camera's full field of view:
    # centre (cx, cy) and width w, with the height fraction equal to w (window
    # and frame share the 2:1 aspect). w == 1 is unzoomed. See apply_roi_delta.
    "roi": {"cx": 0.5, "cy": 0.5, "w": 1.0},
    "precision_mode": "auto",  # one of PRECISION_MODES
}

# The opening view, snapshotted so the UI's "reset frame" can return to it after
# the camera has been flown somewhere unrecoverable.
#
# Deliberately the CAMERA only - where it sits, which way it faces, how narrow
# its field of view is. rs, a and Q are the scene's own parameters rather than
# part of the framing, and they have their own sliders showing their current
# values; a reset that silently pulled those back too would discard a setup the
# user can see they chose, to fix a problem they did not describe.
INITIAL_VIEW = {
    "r0": state["r0"],
    "theta0": state["theta0"],
    "phi0": state["phi0"],
    "orientation": state["orientation"].copy(),
    "roi": dict(state["roi"]),
}

generation = 0
hq_ready = {"generation": -1}

# The accretion-disk flow animation: the ray-traced geometry for a fixed camera
# pose doesn't change over time, only how brightly the disk material each pixel
# looks at is shining at that instant. So once a settle/HQ render supplies a
# shading cache, animate_flow below just re-colorizes the disk pixels on a
# steady clock - no GPU re-render per animation frame.
ANIM_FPS = 15
ANIM_INTERVAL_S = 1.0 / ANIM_FPS
# Instants of the disk averaged into each published frame. A camera integrates
# over the whole time its shutter is open instead of sampling one instant, and
# at this frame rate the inner disk sweeps several pixels between frames -
# without it the material strobes from position to position rather than flowing.
ANIM_SHUTTER_SAMPLES = 2
FLOW_START = time.time()
anim_lock = Lock()
anim_state = {"cache": None, "physics": None, "generation": -1, "frame": 0}

# Frames are cached to disk as PNGs, but a PNG is a finished instant: reopening
# a cached viewpoint left the disk frozen, because nothing on disk can be
# re-colorized at a new time. Holding the last few shading caches in memory
# lets a revisited viewpoint resume flowing immediately instead. Deliberately
# small and bounded - an HQ cache runs to tens of megabytes.
SHADE_CACHE_LIMIT = 6
shade_cache_lock = Lock()
shade_cache: "OrderedDict[str, tuple]" = OrderedDict()


def disk_physics(a: float, Q: float, rs: float = RS) -> cli_imagemaker.DiskPhysics:
    """The scene parameters the shading model needs, as actually traced.

    rs, a and Q all have to be threaded through rather than left at their
    defaults: the disk's orbital angular velocity depends on all three, and that
    velocity is what the flow animation advects at - it has to stay the same one
    the renderer used when it computed every pixel's Doppler shift.

    rs is keyword-with-default rather than positional so the movie scripts, which
    call this as disk_physics(a, Q) against the fixed default, keep working.
    """
    return cli_imagemaker.DiskPhysics(rs=rs, a=a, Q=Q)


def flow_time() -> float:
    return time.time() - FLOW_START


def stash_anim_buffers(cache: cli_imagemaker.ShadingCache, physics: cli_imagemaker.DiskPhysics, gen: int) -> None:
    """Hand a freshly traced frame's shading cache to animate_flow.

    The cache holds the lensed sky and the disk's traced geometry - everything
    that depends on where the camera is rather than on when. The observer never
    moves during flow animation, so this is exactly the work the per-tick path
    must not repeat."""
    with anim_lock:
        anim_state.update(cache=cache, physics=physics, generation=gen)


# rs defaults on all three of these so the movie scripts, which work at the fixed
# default, are unaffected by rs becoming adjustable in the UI.
def horizon_radius(a: float, Q: float, rs: float = RS) -> float:
    """Outer Kerr-Newman horizon r_+; see the sugar_be calculation in ray_step's
    collision setup, which uses the same discriminant."""
    disc = rs * rs - 4 * (a * a + Q * Q)
    if disc <= 0:
        return rs / 2  # extremal edge; clamp_spin_charge keeps us just short of this
    return (rs + math.sqrt(disc)) / 2


def min_r0_for(a: float, Q: float, rs: float = RS) -> float:
    return max(horizon_radius(a, Q, rs) * MIN_R0_PAD, MIN_R0_FLOOR)


def spin_charge_limit(rs: float = RS) -> float:
    """The largest hypot(a, Q) that still leaves a horizon, with margin.

    Reported to the client because the a and Q sliders' own ranges depend on it,
    and it moves whenever rs does.
    """
    return (rs / 2) * SPIN_CHARGE_MARGIN


def clamp_spin_charge(a: float, Q: float, rs: float = RS) -> tuple[float, float]:
    """Keep a^2+Q^2 safely below the extremal (rs/2)^2 bound so a horizon always exists,
    scaling both down proportionally rather than favoring whichever was set last.

    This has to be re-applied when rs *shrinks*, not only when a or Q are moved:
    the bound is proportional to rs, so a spin that was legal at rs = 0.2 is
    superextremal at rs = 0.01 and would leave the hole with no horizon at all.
    """
    limit = spin_charge_limit(rs)
    magnitude = math.hypot(a, Q)
    if magnitude > limit and magnitude > 0:
        scale = limit / magnitude
        a *= scale
        Q *= scale
    return a, Q


def clamp_roi(cx: float, cy: float, w: float) -> dict:
    """Hold the zoom window inside the camera's full field of view.

    A window is never allowed to hang off the edge: there is nothing outside the
    frame to zoom into, the renderer would only trace rays the camera's screen
    does not cover, and letting it drift out is how a pan-while-zoomed ends up
    looking at nothing at all.
    """
    w = min(max(w, 1.0 / MAX_ROI_ZOOM), 1.0)
    half = w / 2.0
    return {"cx": min(max(cx, half), 1.0 - half), "cy": min(max(cy, half), 1.0 - half), "w": w}


def apply_roi_delta(roi: dict, req: "OrbitRequest") -> dict:
    """Fold one request's true-zoom input into the zoom window.

    Everything arriving here is expressed relative to what the client is
    currently *looking at* (the anchor and the box-zoom selection are in
    current-view coordinates, the pan is in units of the current window's
    width), never in absolute frame coordinates. That is deliberate: requests
    are coalesced and can land out of order, so a client that has not yet seen
    the latest window must still be able to say "zoom in on this spot" without
    knowing where that spot sits in the unzoomed frame.
    """
    cx, cy, w = roi["cx"], roi["cy"], roi["w"]
    if req.roi_reset:
        return clamp_roi(0.5, 0.5, 1.0)

    if req.roi_select is not None:
        # Box zoom: the drawn rectangle becomes the whole frame. Its width alone
        # sets the new window, since the window's aspect is fixed - the client
        # is responsible for having snapped the box to that aspect.
        selection = req.roi_select
        left = (cx - w / 2.0) + selection.x * w
        top = (cy - w / 2.0) + selection.y * w
        w = min(max(w * selection.w, 1.0 / MAX_ROI_ZOOM), 1.0)
        cx, cy = left + w / 2.0, top + w / 2.0

    if req.zoom != 1.0 and req.zoom > 0.0:
        # Zoom about the anchor: whatever is under the cursor stays under it, the
        # way an image viewer's wheel zoom does. Clamp the new width first, so a
        # request that runs into the zoom limit still holds the anchor exactly
        # instead of drifting by however much the clamp took off.
        new_w = min(max(w / req.zoom, 1.0 / MAX_ROI_ZOOM), 1.0)
        anchor_x = min(max(req.anchor_x, 0.0), 1.0)
        anchor_y = min(max(req.anchor_y, 0.0), 1.0)
        point_x = (cx - w / 2.0) + anchor_x * w
        point_y = (cy - w / 2.0) + anchor_y * w
        cx = point_x + new_w * (0.5 - anchor_x)
        cy = point_y + new_w * (0.5 - anchor_y)
        w = new_w

    # Pan in units of the current window, so dragging by a given fraction of the
    # view moves by that same fraction of the view at every zoom level.
    cx += req.roi_pan_x * w
    cy += req.roi_pan_y * w
    return clamp_roi(cx, cy, w)


def zoom_accuracy(accuracy: dict, roi: dict) -> dict:
    """The tier's accuracy, tightened for how far the field of view is narrowed."""
    scale = min(1.0 / roi["w"], ZOOM_ACCURACY_CAP)
    return {"errormax": accuracy["errormax"] / scale, "de0": accuracy["de0"]}


def zoom_precision(roi: dict) -> str:
    """Which trace the zoom level needs - see FLOAT_ZOOM_LIMIT for the measurements."""
    return "f32" if 1.0 / roi["w"] <= FLOAT_ZOOM_LIMIT else "f64"


def effective_precision(roi: dict, mode: str) -> str:
    """The trace that will actually run: the zoom's own requirement, or an override."""
    return zoom_precision(roi) if mode == "auto" else mode


def quantize_roi(roi: dict) -> tuple[float, float, float]:
    steps = round(math.log2(1.0 / roi["w"]) / QUANT_ROI_ZOOM_STEP)
    w = 2.0 ** (-steps * QUANT_ROI_ZOOM_STEP)
    centre_step = w * QUANT_ROI_CENTRE
    return round(roi["cx"] / centre_step) * centre_step, round(roi["cy"] / centre_step) * centre_step, w


def cache_key(r0: float, theta0: float, phi0: float, a: float, Q: float, rs: float, omega: tuple,
              roi: dict, precision: str, tier: str, resolution: tuple[int, int]) -> str:
    qr = round(r0 / QUANT_R) * QUANT_R
    qt = round(theta0 / QUANT_ANGLE) * QUANT_ANGLE
    qp = round(phi0 / QUANT_ANGLE) * QUANT_ANGLE
    qa = round(a / QUANT_AQ) * QUANT_AQ
    qq = round(Q / QUANT_AQ) * QUANT_AQ
    qs = round(rs / QUANT_RS) * QUANT_RS
    qo = tuple(round(v / QUANT_OMEGA) * QUANT_OMEGA for v in omega)
    rx, ry, rw = quantize_roi(roi)
    # The precision belongs in the key, not just the ROI: the same view can be
    # traced by either path - at FLOAT_ZOOM_LIMIT, or anywhere at all once the
    # mode is pinned - and those are not the same frame.
    #
    # The ROI is formatted with full round-trip precision, not a fixed number of
    # decimals. Fixed decimals are fine for the camera angles, whose scale never
    # changes, but the zoom window's does: at a billion times magnification the
    # window is ~1e-9 of the frame wide and its centre has to be pinned far
    # tighter than that, so anything less would quietly collapse completely
    # different views onto one cache entry. quantize_roi has already done the
    # deliberate rounding; this must not add its own.
    # rs belongs in the key for the same reason a and Q do: it changes every
    # geodesic in the frame. Adding it changes the string for every view, so the
    # entries already under cache/orbit/ are orphaned rather than mis-served -
    # they are simply dead weight and can be deleted.
    raw = (
        f"v{SHADING_VERSION}:{tier}:{precision}:{resolution[0]}x{resolution[1]}"
        f":{qr:.3f}:{qt:.4f}:{qp:.4f}:{qa:.4f}:{qq:.4f}:{qs:.4f}"
        f":{qo[0]:.3f}:{qo[1]:.3f}:{qo[2]:.3f}"
        f":{rx!r}:{ry!r}:{rw!r}"
    )
    return hashlib.sha256(raw.encode()).hexdigest()[:24]


def publish_image(source: "Image.Image | Path") -> None:
    """Write to IMAGE_PATH via a temp-file-plus-rename so a concurrent GET /image never
    observes a truncated/empty file mid-write. Necessary now that animate_flow publishes
    several times a second on top of the interactive/HQ paths - Image.save and shutil.copy
    both write in place, so without this a poll can land exactly mid-write."""
    fd, tmp_name = tempfile.mkstemp(dir=IMAGE_PATH.parent, suffix=".png")
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        if isinstance(source, Path):
            shutil.copy(source, tmp)
        else:
            source.save(tmp)
        os.replace(tmp, IMAGE_PATH)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def render_frame(worker: "render_worker.RenderWorker", r0: float, theta0: float, phi0: float, a: float, Q: float,
                  rs: float, omega: tuple, roi: dict, precision: str, accuracy: dict,
                  resolution: tuple[int, int], physics: cli_imagemaker.DiskPhysics
                  ) -> tuple[Image.Image, cli_imagemaker.ShadingCache]:
    szeles, magas = resolution
    accuracy = zoom_accuracy(accuracy, roi)
    buffer = worker.render(r0, theta0, phi0, a, Q, rs, accuracy["errormax"], accuracy["de0"], omega,
                           (roi["cx"], roi["cy"], roi["w"]), precision, szeles, magas)
    frame = cli_imagemaker.split_channels(szeles, magas, buffer)
    # The shading cache is built here rather than thrown away and rebuilt for
    # the flow loop: this first image is just its first tick.
    cache = cli_imagemaker.base_scene(frame, RENDER_SEED, physics)
    scene = cli_imagemaker.composite(cache, PEAK_TEMPERATURE, RENDER_SEED, flow_time(), physics)
    return cli_imagemaker.finish_frame(scene, RENDER_SEED), cache


def remember_shading(key: str, cache: cli_imagemaker.ShadingCache, physics: cli_imagemaker.DiskPhysics) -> None:
    with shade_cache_lock:
        shade_cache[key] = (cache, physics)
        shade_cache.move_to_end(key)
        while len(shade_cache) > SHADE_CACHE_LIMIT:
            shade_cache.popitem(last=False)


def render_cached(worker: "render_worker.RenderWorker", r0: float, theta0: float, phi0: float, a: float, Q: float,
                   rs: float, omega: tuple, roi: dict, precision: str, tier: str, accuracy: dict,
                   resolution: tuple[int, int], physics: cli_imagemaker.DiskPhysics
                   ) -> tuple[float, Path, Optional[tuple[cli_imagemaker.ShadingCache, cli_imagemaker.DiskPhysics]]]:
    """Render (or reuse a cached render of) one frame. Returns (elapsed_seconds, png_path, shading).

    shading is what the flow animation needs to keep re-colorizing this viewpoint. A disk cache
    hit still supplies it when the viewpoint is recent enough to be in shade_cache; otherwise it
    is None and the disk holds still on the cached still until a fresh trace supplies one.
    """
    key = cache_key(r0, theta0, phi0, a, Q, rs, omega, roi, precision, tier, resolution)
    cached = CACHE_DIR / f"{key}.png"
    if cached.exists():
        with shade_cache_lock:
            remembered = shade_cache.get(key)
            if remembered is not None:
                shade_cache.move_to_end(key)
        return 0.0, cached, remembered

    start = time.time()
    image, cache = render_frame(worker, r0, theta0, phi0, a, Q, rs, omega, roi, precision, accuracy,
                                resolution, physics)
    elapsed = time.time() - start
    image.save(cached)
    remember_shading(key, cache, physics)
    return elapsed, cached, (cache, physics)


def hq_background_task(gen: int, r0: float, theta0: float, phi0: float, a: float, Q: float, rs: float,
                        omega: tuple, roi: dict, precision: str) -> None:
    time.sleep(HQ_SETTLE_DEBOUNCE_S)
    with state_lock:
        if generation != gen:
            return  # camera moved again before we even started; skip the slow render entirely

    try:
        _, cached_path, shading = render_cached(
            render_worker.HQ_WORKER, r0, theta0, phi0, a, Q, rs, omega, roi, precision, "hq", HQ_ACCURACY,
            HQ_RESOLUTION, disk_physics(a, Q, rs),
        )
    except RuntimeError:
        return

    with state_lock:
        if generation == gen:  # still the same viewpoint -> promote to the live image
            publish_image(cached_path)
            hq_ready["generation"] = gen

    if shading is not None:
        stash_anim_buffers(*shading, gen)


def animate_flow() -> None:
    """Background loop: re-colorize only the disk of the latest settled frame on a steady
    clock, so the accretion disk reads as flowing material instead of a still. The observer's
    view is never touched here - no camera state changes and no GPU ray-tracing call happens
    in this loop, only a CPU recolor of the disk pixels against the already-built shading
    cache (see stash_anim_buffers/base_scene). Skips silently (leaving the last frame on
    screen) whenever the camera has moved past the generation that cache was traced for, or
    nothing has been captured yet.

    The shutter is set to the interval actually achieved rather than the one aimed for, so
    the motion blur matches the rate frames are really being published at: if a big HQ frame
    makes the loop run at half speed, each frame integrates over twice as long, exactly as a
    camera shooting at that rate would."""
    interval = ANIM_INTERVAL_S
    while True:
        started = time.time()
        with anim_lock:
            cache, physics, gen = anim_state["cache"], anim_state["physics"], anim_state["generation"]
            index = anim_state["frame"]
            anim_state["frame"] = index + 1
        if cache is not None:
            with state_lock:
                current = generation
            if current == gen:  # otherwise the camera moved on; wait for the next settle
                scene = cli_imagemaker.composite(
                    cache, PEAK_TEMPERATURE, RENDER_SEED, flow_time(), physics,
                    shutter=interval, shutter_samples=ANIM_SHUTTER_SAMPLES,
                )
                image = cli_imagemaker.finish_frame(scene, RENDER_SEED, index)
                with state_lock:
                    if generation == gen:  # still current; publish
                        publish_image(image)
        spent = time.time() - started
        interval = max(spent, ANIM_INTERVAL_S)
        time.sleep(max(0.0, ANIM_INTERVAL_S - spent))


threading.Thread(target=animate_flow, daemon=True).start()


class RoiSelect(BaseModel):
    """A box-zoom rectangle, in fractions of what the client is currently viewing.

    Only the width is taken: the window's height follows from the fixed 2:1
    aspect, so the client snaps the drawn box to that aspect before sending it.
    """
    x: float
    y: float
    w: float


class OrbitRequest(BaseModel):
    dtheta: float = 0.0
    dphi: float = 0.0
    dzoom: float = 1.0  # DOLLY: multiplicative factor applied to r0, i.e. the camera moves
    dpan: float = 0.0    # camera orientation deltas (independent of r0/theta0/phi0):
    dtilt: float = 0.0   # pan turns around the camera's current up axis,
    droll: float = 0.0   # tilt around its current right axis, roll around its current forward axis
    # TRUE ZOOM: narrows the field of view instead of moving the camera. All of
    # these are relative to the client's current view, never to the unzoomed
    # frame - see apply_roi_delta.
    zoom: float = 1.0        # >1 zooms in, <1 out
    anchor_x: float = 0.5    # the point `zoom` holds fixed (0..1 across the current view)
    anchor_y: float = 0.5
    roi_select: Optional[RoiSelect] = None  # box zoom: make this rectangle the whole frame
    roi_pan_x: float = 0.0   # slide the zoom window, in units of its own width/height
    roi_pan_y: float = 0.0
    roi_reset: bool = False  # back to the full field of view
    # Back to the opening camera: position, orientation and zoom together. rs, a
    # and Q are left alone - see INITIAL_VIEW.
    reset_view: bool = False
    a: Optional[float] = None  # absolute spin; omitted/None means "leave unchanged"
    Q: Optional[float] = None  # absolute charge; omitted/None means "leave unchanged"
    # Absolute Schwarzschild radius, same convention. Moving it rescales the hole
    # against the disk and camera, and drags the a/Q bound with it - see
    # MIN_RS/MAX_RS and clamp_spin_charge.
    rs: Optional[float] = None
    # "auto", "f32" or "f64"; omitted/None means "leave unchanged". Anything
    # else is ignored rather than rejected - a stale client should not be able
    # to break the render loop over a field it merely got wrong.
    precision_mode: Optional[str] = None
    dragging: bool = True
    hq: bool = False  # opt-in background HQ pass on settle; off by default since it can take ~15s


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
        # First, so that a request combining the reset with movement deltas
        # applies them to the restored view rather than having them overwritten.
        # The r0 clamp further down re-applies afterwards, which matters: a large
        # rs can push min_r0 past where the opening view sat.
        if req.reset_view:
            state["r0"] = INITIAL_VIEW["r0"]
            state["theta0"] = INITIAL_VIEW["theta0"]
            state["phi0"] = INITIAL_VIEW["phi0"]
            state["orientation"] = INITIAL_VIEW["orientation"].copy()
            state["roi"] = dict(INITIAL_VIEW["roi"])

        state["theta0"] = min(max(state["theta0"] + req.dtheta, POLE_EPS), math.pi - POLE_EPS)
        state["phi0"] = (state["phi0"] + req.dphi) % (2 * math.pi)
        state["orientation"] = apply_look_delta(state["orientation"], req.dpan, req.dtilt, req.droll)

        if req.rs is not None:
            state["rs"] = min(max(req.rs, MIN_RS), MAX_RS)

        # Re-clamped when rs moves as well as when a or Q do: the extremal bound
        # is proportional to rs, so shrinking rs alone can leave a spin that was
        # legal a moment ago superextremal, and the hole with no horizon.
        if req.a is not None or req.Q is not None or req.rs is not None:
            new_a = req.a if req.a is not None else state["a"]
            new_q = req.Q if req.Q is not None else state["Q"]
            state["a"], state["Q"] = clamp_spin_charge(new_a, new_q, state["rs"])

        # Likewise the horizon, and so the closest the camera may sit: growing rs
        # can swallow an r0 that was outside it.
        min_r0 = min_r0_for(state["a"], state["Q"], state["rs"])
        dolly = min(max(req.dzoom, 0.5), 2.0)
        state["r0"] = min(max(state["r0"] * dolly, min_r0), MAX_R0)
        state["roi"] = apply_roi_delta(state["roi"], req)
        roi = dict(state["roi"])
        if req.precision_mode in PRECISION_MODES:
            state["precision_mode"] = req.precision_mode
        precision = effective_precision(roi, state["precision_mode"])
        generation += 1
        gen = generation
        omega = current_omega(state["orientation"])

        tier = "drag" if req.dragging else "settle"
        accuracy = DRAG_ACCURACY if req.dragging else SETTLE_ACCURACY
        resolution = DRAG_RESOLUTION if req.dragging else SETTLE_RESOLUTION
        elapsed, cached_path, shading = render_cached(
            render_worker.INTERACTIVE_WORKER, state["r0"], state["theta0"], state["phi0"], state["a"], state["Q"],
            state["rs"], omega, roi, precision, tier, accuracy, resolution,
            disk_physics(state["a"], state["Q"], state["rs"]),
        )
        publish_image(cached_path)

        if not req.dragging and req.hq:
            threading.Thread(
                target=hq_background_task,
                args=(gen, state["r0"], state["theta0"], state["phi0"], state["a"], state["Q"],
                      state["rs"], omega, roi, precision),
                daemon=True,
            ).start()

        result = {
            "elapsed": elapsed,
            "r0": state["r0"],
            "theta0": state["theta0"],
            "phi0": state["phi0"],
            "a": state["a"],
            "Q": state["Q"],
            "rs": state["rs"],
            # The a and Q sliders' own range moves with rs, so the client is told
            # what it is rather than left to recompute the extremal bound.
            "spin_charge_limit": spin_charge_limit(state["rs"]),
            "min_rs": MIN_RS,
            "max_rs": MAX_RS,
            "min_r0": min_r0,
            "roi": roi,
            "zoom": 1.0 / roi["w"],
            "precision": precision,                       # what actually traced this frame
            "precision_mode": state["precision_mode"],    # and whether that was chosen or forced
            "float_zoom_limit": FLOAT_ZOOM_LIMIT,
            "generation": gen,
            "image": f"/image?t={time.time()}",
        }

    # Outside state_lock: nothing here needs the camera-state lock, it only reads the
    # shading cache/gen already captured above. Drag-tier frames are too transient and
    # low-res to bother animating; only settle frames seed the flow loop (HQ renders also
    # seed it, at higher quality, once ready).
    if not req.dragging and shading is not None:
        stash_anim_buffers(*shading, gen)

    return result


@app.get("/status")
def status():
    with state_lock:
        return {"generation": generation, "hq_ready": hq_ready["generation"] == generation}


@app.get("/image")
def image():
    # Not FileResponse: it stats the file for Content-Length and opens/reads it as two
    # separate steps, which can race publish_image's atomic replace (a differently-sized
    # frame swapped in between the stat and the read -> "Too little data for declared
    # Content-Length"). Reading the bytes in one call instead means the Content-Length is
    # derived from what was actually read - os.replace() is atomic at the directory-entry
    # level, so an in-flight read always sees one complete file's contents, never a mix.
    return Response(IMAGE_PATH.read_bytes(), media_type="image/png", headers={"Cache-Control": "no-store"})


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)
