"""Shared plumbing for the movie scripts: tracing, shading, captions, GIF output.

Everything scene-related is taken from web_server rather than restated, for the
same reason photon_ring_zoom.py does it: a movie that hardcodes its own copy of
RS, the accuracy tiers or the zoom/precision policy silently stops showing the
same hole as the interactive view the moment either is touched.  Importing that
module also starts its flow-animation thread, which idles harmlessly in a script
because nothing ever hands it a frame to re-colorize.

The helpers here are the ones photon_ring_zoom.py grew and the newer movies need
too - font fitting, captioning, GIF encoding, the trace/shade pair.  That script
predates this module and still carries its own copies; it is left alone rather
than migrated, since it works and rewriting a finished movie script to prove a
point is not worth the risk of changing what it renders.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

import cli_imagemaker
import web_server as scene


FONT_CANDIDATES = (
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/System/Library/Fonts/Helvetica.ttc",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/TTF/DejaVuSans.ttf",
)

# ffi_bridge.cpp hardcodes the outer integration sphere at sugar_ki = 1.01, and
# ray_step_T ends any ray already outside it.  A camera placed beyond that
# radius therefore has every one of its rays terminate on the first step, and
# the frame comes back 100% "escaped" - a blank starfield with no hole in it.
# Measured directly: r0 = 1.01 renders normally, r0 = 1.02 renders empty.
#
# web_server's MAX_R0 is 6.0 and does not know this, so the interactive UI can
# dolly out to a blank screen too; that is worth fixing there, but a movie
# script should at minimum refuse to spend an hour rendering nothing.
SUGAR_KI = 1.01
MAX_USEFUL_R0 = 1.0  # a hair inside it, so the camera is never on the boundary


# --------------------------------------------------------------------------
# Tracing and shading
# --------------------------------------------------------------------------


def trace(worker, camera: dict, roi: dict, accuracy: dict, resolution: tuple,
          omega: tuple = None, precision: str = None) -> cli_imagemaker.RayFrame:
    """One frame of raw ray data: the 6 channels cuda_ray.h's ray_step_T writes.

    `omega` defaults to the same 180-degree flip web_server opens with, so a
    caller that is not moving the camera's gaze does not have to think about
    orientation at all.  `precision` defaults to web_server's own zoom-driven
    policy (f32 until the window is narrower than FLOAT_ZOOM_LIMIT, f64 after).
    """
    szeles, magas = resolution
    if omega is None:
        omega = scene.current_omega(scene.BASE_ORIENTATION)
    if precision is None:
        precision = scene.effective_precision(roi, "auto")
    buffer = worker.render(
        camera["r0"], camera["theta0"], camera["phi0"],
        camera["a"], camera["Q"], scene.RS,
        accuracy["errormax"], accuracy["de0"],
        omega, (roi["cx"], roi["cy"], roi["w"]), precision, szeles, magas,
    )
    return cli_imagemaker.split_channels(szeles, magas, buffer)


def shade(frame: cli_imagemaker.RayFrame, physics, flow_time: float = 0.0,
          frame_index: int = 0) -> Image.Image:
    """The same compositing path the interactive view uses, at one instant."""
    cache = cli_imagemaker.base_scene(frame, scene.RENDER_SEED, physics)
    lit = cli_imagemaker.composite(cache, scene.PEAK_TEMPERATURE, scene.RENDER_SEED,
                                   flow_time, physics)
    return cli_imagemaker.finish_frame(lit, scene.RENDER_SEED, frame_index)


# --------------------------------------------------------------------------
# Finding the shadow's edge
# --------------------------------------------------------------------------
#
# Any movie that zooms has to converge on its target far more precisely than the
# final frame is wide, or the subject drifts out of shot.  The shadow's edge is
# not taken from theory here: it is bisected for, by tracing single rays and
# asking which ones the renderer captures, so the target is the edge this
# renderer actually draws.  (photon_ring_zoom.py carries its own copy of this
# from before this module existed; the two agree.)

# A probe window narrow enough that both of its two pixels are the same ray to
# any relevant precision - this is how one point of the image plane gets traced
# through an interface that only renders whole frames.
PROBE_WINDOW = 1e-12
PROBE_RESOLUTION = (2, 1)  # the smallest frame keeping the required 2:1 aspect


def captured_at(worker, camera: dict, cx: float, cy: float, accuracy: dict,
                precision: str = "f64") -> bool:
    """Does the ray through this point of the image plane fall into the hole?

    Probed in f64 by default. The edge is bisected far past what f32 can
    resolve, and an f32 probe would stop discriminating partway down and hand
    back whichever side its last resolvable step landed on.
    """
    roi = {"cx": cx, "cy": cy, "w": PROBE_WINDOW}
    frame = trace(worker, camera, roi, accuracy, PROBE_RESOLUTION, precision=precision)
    return bool(frame.radius.flat[0] == 0.0)


def scan_column(worker, camera: dict, cx: float, accuracy: dict, samples: int,
                precision: str = "f64") -> list:
    """Every place the column at cx crosses between sky and shadow.

    Seeds for bisection are found rather than assumed: spin drags the shadow
    sideways, and a column away from centre meets it at different rows or not at
    all, so a guessed bracket straddles nothing or the wrong edge.

    Returns (upper_cy, lower_cy, captured_below) brackets, ordered top to bottom.
    """
    rows = [(index + 0.5) / samples for index in range(samples)]
    flags = [captured_at(worker, camera, cx, cy, accuracy, precision) for cy in rows]
    return [(rows[i], rows[i + 1], flags[i + 1])
            for i in range(samples - 1) if flags[i] != flags[i + 1]]


def find_shadow_edge(worker, camera: dict, cx: float, bracket: tuple, accuracy: dict,
                     steps: int = 48, precision: str = "f64") -> float:
    """Bisect a bracket known to straddle the shadow's edge down to the grid's limit."""
    upper, lower, captured_below = bracket
    sky_cy, shadow_cy = (upper, lower) if captured_below else (lower, upper)
    for _ in range(steps):
        middle = (sky_cy + shadow_cy) / 2.0
        if captured_at(worker, camera, cx, middle, accuracy, precision):
            shadow_cy = middle
        else:
            sky_cy = middle
    return (sky_cy + shadow_cy) / 2.0


# --------------------------------------------------------------------------
# Captions
# --------------------------------------------------------------------------


def load_font(size: int) -> ImageFont.ImageFont:
    for path in FONT_CANDIDATES:
        if Path(path).exists():
            try:
                return ImageFont.truetype(path, size)
            except OSError:
                continue
    return ImageFont.load_default()


def fit_font(draw: ImageDraw.ImageDraw, lines: list, width: int, start_size: int):
    """Largest candidate size whose longest line fits, or None if none do.

    Worth doing rather than picking a size from the frame height alone: at a
    fixed size a narrow frame silently truncates the caption mid-word, which is
    how a frame ends up labelled "traced in fp3".
    """
    for size in range(start_size, 6, -1):
        font = load_font(size)
        if max(draw.textlength(line, font=font) for line in lines) <= width - 12:
            return font
    return None


def caption_line_height(size: int) -> int:
    """Vertical advance draw_caption uses for one line at this font size."""
    return size + 4


def draw_caption(canvas: Image.Image, lines: list, origin: tuple = (6, 6),
                 width: int = None, colour: tuple = (235, 235, 235),
                 max_size: int = None) -> Image.Image:
    """Drop-shadowed text, sized to fit, drawn in place.  Returns the canvas.

    `max_size` caps the font, for callers that reserved a strip of a known
    height for the text: without it the size scales with the canvas, and a
    caller that sized its strip from a guess gets its last line clipped once
    the render is big enough.
    """
    if not lines:
        return canvas
    draw = ImageDraw.Draw(canvas)
    width = canvas.width if width is None else width
    start_size = max(11, canvas.height // 22)
    if max_size is not None:
        start_size = min(start_size, max_size)
    font = fit_font(draw, lines, width, start_size) or load_font(8)
    x, y = origin
    for line in lines:
        # Plain drop shadow, so text stays legible over both the bright ring and
        # the black shadow it sits against.
        draw.text((x + 1, y + 1), line, font=font, fill=(0, 0, 0))
        draw.text((x, y), line, font=font, fill=colour)
        y += caption_line_height(getattr(font, "size", 10))
    return canvas


def format_zoom(zoom: float) -> str:
    if zoom < 10_000:
        return f"{zoom:,.0f}x"
    exponent = int(math.floor(math.log10(zoom)))
    return f"{zoom / 10**exponent:.1f}e{exponent}x"


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------


def write_gif(images: list, output: Path, ms_per_frame: int, pingpong: bool = False) -> None:
    """Encode a GIF with per-frame adaptive palettes.

    Per-frame rather than one global palette because these frames are nearly all
    near-black: a single 256-colour palette spends most of its entries on the
    background and bands the disk and ring badly.
    """
    ordered = images + images[-2:0:-1] if pingpong else images
    quantized = [im.convert("RGB").quantize(colors=256, method=Image.Quantize.FASTOCTREE)
                 for im in ordered]
    output.parent.mkdir(parents=True, exist_ok=True)
    quantized[0].save(output, save_all=True, append_images=quantized[1:],
                      duration=ms_per_frame, loop=0, optimize=False, disposal=2)
    print(f"wrote {output} ({len(quantized)} frames, {output.stat().st_size / 1e6:.1f} MB)")


# --------------------------------------------------------------------------
# Shared argument surface
# --------------------------------------------------------------------------


def add_scene_args(parser: argparse.ArgumentParser) -> None:
    """The camera/metric group, identical across the movie scripts."""
    group = parser.add_argument_group(
        "scene", "Defaults match web_server's opening view. a and Q are in the same code "
                 "units as the web UI's sliders, where the extremal bound is rs/2 = "
                 f"{scene.RS / 2:g}; they are scaled down together if they exceed it.")
    group.add_argument("--a", type=float, default=0.0, help="black hole spin")
    group.add_argument("--Q", type=float, default=0.0, help="black hole charge")
    group.add_argument("--r0", type=float, default=scene.state["r0"],
                       help=f"camera radius (default {scene.state['r0']:g}; must stay inside the "
                            f"hardcoded outer sphere sugar_ki={SUGAR_KI:g} or every ray escapes "
                            "immediately and the frame comes back empty)")
    group.add_argument("--theta0", type=float, default=scene.state["theta0"],
                       help="camera polar angle")
    group.add_argument("--phi0", type=float, default=scene.state["phi0"], help="camera azimuth")


def resolve_scene(args) -> tuple:
    """Clamp the scene arguments the way the interactive UI does, and report.

    Returns (camera dict, physics).  Raises SystemExit on an r0 that would
    render nothing, rather than letting the movie run to completion empty.
    """
    if args.r0 >= SUGAR_KI:
        raise SystemExit(
            f"--r0 {args.r0:g} is outside the renderer's outer integration sphere "
            f"(sugar_ki = {SUGAR_KI:g}, hardcoded in ffi_bridge.cpp). Every ray would terminate "
            f"on its first step and every frame would be an empty starfield. Use r0 <= "
            f"{MAX_USEFUL_R0:g}.")
    # Clamp first, then check r0 against the horizon the clamped values actually
    # produce - checking against the raw ones would test a hole that is not the
    # one about to be rendered.
    spin, charge = scene.clamp_spin_charge(args.a, args.Q)
    if (spin, charge) != (args.a, args.Q):
        print(f"a and Q scaled to the extremal bound: a={spin:.6g}, Q={charge:.6g}")

    horizon_floor = scene.min_r0_for(spin, charge)
    if args.r0 <= horizon_floor:
        raise SystemExit(f"--r0 {args.r0:g} is at or inside the horizon "
                         f"(needs > {horizon_floor:g} for a={spin:g}, Q={charge:g})")
    camera = {"r0": args.r0, "theta0": args.theta0, "phi0": args.phi0, "a": spin, "Q": charge}
    mass = scene.RS / 2.0
    print(f"scene: a={spin:.6g} (a/M={spin / mass:.3f})  Q={charge:.6g} (Q/M={charge / mass:.3f})  "
          f"camera r0={camera['r0']:g} theta0={camera['theta0']:g} phi0={camera['phi0']:g}")
    return camera, scene.disk_physics(spin, charge)


def clamp_polar(theta: float) -> float:
    """Keep theta off the coordinate axis.

    cuda_ray.h floors |sin(theta)| in the phi Christoffel term, but the tetrad in
    ijk_to_vec_zoom divides by sin(theta_0) outright, so the camera itself has to
    stay clear of the pole regardless.
    """
    return min(max(theta, scene.POLE_EPS), math.pi - scene.POLE_EPS)


def ease(t: float) -> float:
    """Smoothstep, so a sweep starts and ends at rest instead of jerking."""
    return t * t * (3.0 - 2.0 * t)


def hit_fractions(frame: cli_imagemaker.RayFrame) -> tuple:
    """(escaped, disk, captured) as fractions of the frame - a cheap sanity read."""
    return (float((frame.radius < 0).mean()),
            float((frame.radius > 0).mean()),
            float((frame.radius == 0).mean()))
