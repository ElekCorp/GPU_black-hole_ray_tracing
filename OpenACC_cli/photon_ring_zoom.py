#!/usr/bin/env python3
"""Render a GIF that zooms into the photon ring, one photon orbit at a time.

WHAT IS ACTUALLY ON SCREEN

Light that passes close enough to the hole loops around it before escaping to
the camera. How many times it loops depends on how close: rays arriving just
outside the critical impact parameter loop once, rays a little closer loop
twice, and so on without limit - the loop count diverges as the critical value
is approached from outside. Every one of those families is drawn in the image
as its own thin arc just outside the shadow's edge, each one nested inside the
last. They are the photon ring's subrings, and for a Schwarzschild hole each is
thinner than the one before it by a factor of e^(2*pi) ~ 535.

So they are unreachable by rendering a bigger image: the second subring is
already finer than one pixel of any frame that shows the whole hole, and no
amount of resolution helps because the detail is not in the frame's sampling,
it is outside its field of view budget. What reaches them is narrowing the
field of view onto the shadow's edge and retracing - which is what this script
does, geometrically, over the frames of a GIF. Each factor of ~535 in
magnification buys one more orbit.

HOW THE TARGET IS FOUND

The zoom has to converge on the shadow's edge to a precision far finer than the
final frame is wide, or the ring drifts out of shot. The edge is not assumed
from theory: it is bisected for, by tracing individual rays and asking which
ones the renderer captures, so the target is exactly the edge this renderer
draws rather than an idealisation of it. The camera never moves - the edge sits
at one fixed place in the image plane - so this is done once and every frame
zooms toward the same point.

One honesty note about that edge. cuda_ray.h gives each ray a step budget of
int(1/errormax), and a ray that exhausts it is recorded as having escaped. Rays
asymptotically close to the critical value orbit indefinitely, so some of them
run out of budget and are drawn as sky rather than shadow. The boundary found
here is therefore the renderer's shadow edge at the accuracy being used, which
approaches the true critical curve as errormax is tightened. That is the right
target anyway: it is the edge of what the frames actually show.

The per-frame loop count is measured, not predicted. The renderer records each
ray's elapsed coordinate time, so the orbit count is read off as how far that
exceeds a direct ray's, in units of the coordinate period of the circular
photon orbit.

HOW DEEP IT IS WORTH GOING

Not as deep as the renderer will zoom. dopri54_step clamps its step tolerance
to fmax(errormax, 1e-8), and each further orbit needs the impact parameter
resolved ~535 times more finely, which leaves about ln(1e8)/(2*pi) ~ 2.9 orbits
of headroom past the direct image - four to five in total, which is what the
frames report. Past ~1e6 magnification the edge being zoomed into is no longer
the critical curve but the locus where rays exhaust their step budget, and the
lensed disk leaves the frame: measured, it covers 70% of the frame at 1e6x and
0% at 1e8x. The default stops before that. The bottleneck is the integrator's
tolerance floor, not the zoom, the ROI grid, or fp64.

Usage:
    python photon_ring_zoom.py                       # 48 frames to 1e6x
    python photon_ring_zoom.py --frames 96 --width 960 --pingpong
    python photon_ring_zoom.py --max-zoom 1e8 -o deep.gif    # past the useful range
"""

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

import cli_imagemaker
import render_worker
# The renderer's constants and its zoom/precision policy are taken from the web
# UI rather than restated, so the movie is guaranteed to be showing the same
# hole, shaded the same way, as the interactive view. Importing it also starts
# that module's flow-animation thread, which idles harmlessly here because
# nothing ever hands it a frame to re-colorize.
import web_server as scene

# A probe window narrow enough that both of its two pixels are, to any relevant
# precision, the same ray - this is how a single point in the image plane gets
# traced through an interface that only renders whole frames.
PROBE_WINDOW = 1e-12
PROBE_RESOLUTION = (2, 1)  # the smallest frame keeping the required 2:1 aspect

# dopri54_step clamps the step tolerance it actually uses to fmax(errormax,
# 1e-8), so nothing below this changes a trajectory - it only buys step budget,
# which is int(1/errormax) and is not clamped. Measured directly: errormax 1e-8
# and 1e-9 produce byte-identical frames and an identical shadow edge.
TOLERANCE_FLOOR = 1e-8

# How far the zoom can usefully go, and why it stops there. Each extra orbit
# requires resolving an impact parameter e^(2*pi) ~ 535 times more finely, so
# from a tolerance of 1e-8 there are only about ln(1e8)/(2*pi) ~ 2.9 orbits of
# headroom past the direct image. Beyond roughly 1e6 magnification the frame
# stops containing lensed disk at all - measured: the disk covers 70% of the
# frame at 1e6x and 0% at 1e8x - because the shadow edge being zoomed into is
# by then the locus where rays exhaust their budget rather than the true
# critical curve, and it sits outside what is left of the cascade. This is a
# limit of the integrator, not of the zoom machinery or of fp64.
USEFUL_MAX_ZOOM = 1e6

FONT_CANDIDATES = (
    "/System/Library/Fonts/Supplemental/Arial.ttf",
    "/System/Library/Fonts/Helvetica.ttc",
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    "/usr/share/fonts/TTF/DejaVuSans.ttf",
)


def photon_orbit_period(rs: float) -> float:
    """Coordinate time for one turn of the circular photon orbit at r = 3M.

    For Schwarzschild, a null circular orbit has d(phi)/dt = sqrt(M/r^3), so at
    r = 3M the period is 2*pi*3^(3/2)*M. This is the unit the ring's subrings
    are spaced in, and the unit an orbit count is reported in below.
    """
    mass = rs / 2.0
    return 2.0 * math.pi * 3.0**1.5 * mass


def trace(worker, camera: dict, roi: dict, accuracy: dict, resolution: tuple) -> cli_imagemaker.RayFrame:
    """One frame of raw ray data for a zoom window, at the precision it needs.

    Note this does NOT apply web_server's zoom_accuracy scaling. That exists to
    tighten an interactive frame as the view narrows, from a loose starting
    point and under a cap; here the accuracy is already pinned at the tightest
    setting the integrator responds to, for every frame including the widest,
    so scaling it could only add confusion. Crucially the *same* accuracy is
    used for the frames and for the edge bisection: the shadow edge moves with
    accuracy (measured: by 0.0016 of the frame between errormax 9.8e-7 and
    1e-8, which is a thousand frame-widths at the zooms involved here), so a
    mismatch would centre every frame on the wrong point.
    """
    szeles, magas = resolution
    buffer = worker.render(
        camera["r0"], camera["theta0"], camera["phi0"],
        camera["a"], camera["Q"], scene.RS,
        accuracy["errormax"], accuracy["de0"],
        scene.current_omega(scene.BASE_ORIENTATION),
        (roi["cx"], roi["cy"], roi["w"]),
        scene.effective_precision(roi, "auto"), szeles, magas,
    )
    return cli_imagemaker.split_channels(szeles, magas, buffer)


def captured_at(worker, camera: dict, cx: float, cy: float, accuracy: dict) -> bool:
    """Does the ray through this point of the image plane fall into the hole?"""
    roi = {"cx": cx, "cy": cy, "w": PROBE_WINDOW}
    frame = trace(worker, camera, roi, accuracy, PROBE_RESOLUTION)
    return bool(frame.radius.flat[0] == 0.0)


def scan_column(worker, camera: dict, cx: float, accuracy: dict, samples: int) -> list:
    """Every place the column at cx crosses between sky and shadow.

    The seeds for bisection are found rather than assumed, because where the
    shadow sits is exactly what changes in the cases this is most wanted for. A
    column away from the frame's centre meets the shadow at different rows, or
    not at all; spin drags the shadow sideways and flattens it on the prograde
    side, so its edge is no longer where it was at a = 0. Guessing a bracket
    then either straddles nothing or straddles the wrong edge.

    Returns (upper_cy, lower_cy, captured_below) brackets, ordered top to
    bottom, each known to contain exactly one crossing.
    """
    rows = [(index + 0.5) / samples for index in range(samples)]
    flags = [captured_at(worker, camera, cx, cy, accuracy) for cy in rows]
    return [(rows[i], rows[i + 1], flags[i + 1])
            for i in range(samples - 1) if flags[i] != flags[i + 1]]


def columns_crossing_shadow(worker, camera: dict, accuracy: dict, samples: int = 24) -> tuple:
    """The span of columns that meet the shadow at all, for when the chosen one misses.

    Easy to misjudge by eye, which is why it is worth reporting rather than
    leaving the user to guess: cx and cy are fractions of their own axis and the
    frame is 2:1, so a round shadow covers only half as much of the horizontal
    range as it does of the vertical one. A column that looks comfortably inside
    the silhouette on screen can still miss it in these coordinates.
    """
    hits = [cx for cx in ((index + 0.5) / samples for index in range(samples))
            if scan_column(worker, camera, cx, accuracy, 32)]
    return (min(hits), max(hits)) if hits else (None, None)


def find_shadow_edge(worker, camera: dict, cx: float, bracket: tuple,
                     accuracy: dict, steps: int) -> float:
    """Bisect a bracket known to straddle the shadow's edge down to the grid's limit."""
    upper, lower, captured_below = bracket
    sky_cy, shadow_cy = (upper, lower) if captured_below else (lower, upper)
    for _ in range(steps):
        middle = (sky_cy + shadow_cy) / 2.0
        if captured_at(worker, camera, cx, middle, accuracy):
            shadow_cy = middle
        else:
            sky_cy = middle
    return (sky_cy + shadow_cy) / 2.0


def elapsed_times(frame: cli_imagemaker.RayFrame) -> np.ndarray:
    """Coordinate time elapsed along each ray that reached somewhere, non-finite dropped.

    Rays that fell through the horizon are excluded: they stop at the hole
    rather than at the camera, so their elapsed time is not a path to anywhere.

    The non-finite filter is not defensive habit. cuda_ray.h has a bail-out for
    a ray whose state has gone non-finite, and it records the last good time -
    which is itself non-finite if the ray went bad immediately. Those turn up
    exactly where this movie spends its time, a hair from the critical curve,
    and one of them is enough to poison a percentile; worse, max(0.0, nan)
    silently returns 0.0 in Python, so an unfiltered frame reports no orbits at
    all rather than reporting a problem.
    """
    elapsed = np.abs(frame.travel_time[frame.radius != 0.0])
    return elapsed[np.isfinite(elapsed)]


def orbit_count(frame: cli_imagemaker.RayFrame, direct_time: float, period: float) -> float:
    """The most turns any light in this frame made, from its elapsed coordinate time.

    The most-wound ray is the point of the frame, so this is a high percentile
    rather than a median - a median just reports the ordinary light that happens
    to share the shot. It is a percentile rather than the outright maximum
    because right at the critical curve there is always a ray or two that ran
    away with the step budget, and one of those should not caption the frame.
    """
    elapsed = elapsed_times(frame)
    if elapsed.size == 0:
        return 0.0
    return max(0.0, float(np.percentile(elapsed, 99) - direct_time) / period)


def shade(frame: cli_imagemaker.RayFrame, physics, flow_time: float) -> Image.Image:
    """The same compositing the interactive view uses, at a fixed instant.

    The disk's flow is held at one time across the whole movie on purpose: the
    subject here is the zoom, and letting the material advect between frames
    would put motion in the disk competing with it.
    """
    cache = cli_imagemaker.base_scene(frame, scene.RENDER_SEED, physics)
    lit = cli_imagemaker.composite(cache, scene.PEAK_TEMPERATURE, scene.RENDER_SEED, flow_time, physics)
    return cli_imagemaker.finish_frame(lit, scene.RENDER_SEED)


def load_font(size: int) -> ImageFont.ImageFont:
    for path in FONT_CANDIDATES:
        if Path(path).exists():
            try:
                return ImageFont.truetype(path, size)
            except OSError:
                continue
    return ImageFont.load_default()


def format_zoom(zoom: float) -> str:
    if zoom < 10_000:
        return f"{zoom:,.0f}x"
    exponent = int(math.floor(math.log10(zoom)))
    return f"{zoom / 10**exponent:.1f}e{exponent}x"


def caption(zoom: float, orbits: float, precision: str, compact: bool) -> list:
    if compact:
        lines = [format_zoom(zoom), f"~{orbits:.1f} orbits"]
        return lines + ["fp64"] if precision == "f64" else lines
    lines = [
        f"field of view narrowed {format_zoom(zoom)}",
        f"light here has looped up to ~{orbits:.1f}x around the hole",
    ]
    return lines + ["traced in fp64"] if precision == "f64" else lines


def fit_font(draw: ImageDraw.ImageDraw, lines: list, width: int, start_size: int):
    """Largest of our candidate sizes whose longest line fits, or None if none do.

    Worth doing rather than picking a size from the height alone: the caption is
    the only thing saying how deep the zoom is, and at a fixed size a narrow
    frame silently truncates it mid-sentence - which is how a frame ends up
    captioned "looped up to ~4.6x around the hol".
    """
    for size in range(start_size, 7, -1):
        font = load_font(size)
        if max(draw.textlength(line, font=font) for line in lines) <= width - 12:
            return font
    return None


def annotate(image: Image.Image, zoom: float, orbits: float, precision: str) -> Image.Image:
    canvas = image.convert("RGB")
    draw = ImageDraw.Draw(canvas)
    start_size = max(11, canvas.height // 22)
    # Below roughly 200px wide the sentence does not fit at any legible size, so
    # drop to the bare numbers rather than shrinking until it is unreadable.
    lines = caption(zoom, orbits, precision, compact=False)
    font = fit_font(draw, lines, canvas.width, start_size)
    if font is None:
        lines = caption(zoom, orbits, precision, compact=True)
        font = fit_font(draw, lines, canvas.width, start_size) or load_font(8)
    y = 6
    for line in lines:
        # A plain drop shadow, so the text stays readable over both the bright
        # ring and the black shadow it sits against.
        draw.text((7, y + 1), line, font=font, fill=(0, 0, 0))
        draw.text((6, y), line, font=font, fill=(235, 235, 235))
        y += getattr(font, "size", 10) + 4
    return canvas


def render_movie(worker, camera: dict, args, target_cy: float, accuracy: dict,
                  resolution: tuple, physics, period: float) -> int:
    """Trace the zoom, frame by frame, and write the GIF."""
    # A direct ray's elapsed time, taken from the widest frame, is the baseline
    # every later frame's orbit count is measured against.
    wide = trace(worker, camera, {"cx": 0.5, "cy": 0.5, "w": 1.0}, accuracy, resolution)
    wide_elapsed = elapsed_times(wide)
    if wide_elapsed.size == 0:
        raise SystemExit("no ray in the unzoomed frame reached anywhere - cannot set a baseline")
    direct_time = float(wide_elapsed.min())
    print(f"  direct-ray coordinate time {direct_time:.3f}, photon-orbit period {period:.3f}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    frames = []
    for index in range(args.frames):
        fraction = index / (args.frames - 1)
        zoom = args.max_zoom**fraction  # geometric, so each frame gains the same factor
        roi = scene.clamp_roi(args.cx, target_cy, 1.0 / zoom)
        precision = scene.effective_precision(roi, "auto")

        started = time.time()
        frame = trace(worker, camera, roi, accuracy, resolution)
        orbits = orbit_count(frame, direct_time, period)
        frames.append(annotate(shade(frame, physics, 0.0), zoom, orbits, precision))
        print(f"  frame {index + 1:3d}/{args.frames}  {format_zoom(zoom):>10}  "
              f"{precision}  ~{orbits:5.2f} orbits  {time.time() - started:5.1f}s")

    ordered = frames + frames[-2:0:-1] if args.pingpong else frames
    # Per-frame adaptive palettes: these frames are nearly all near-black, and a
    # single global 256-colour palette bands the ring badly.
    quantized = [image.quantize(colors=256, method=Image.Quantize.FASTOCTREE) for image in ordered]
    quantized[0].save(args.output, save_all=True, append_images=quantized[1:],
                      duration=args.ms_per_frame, loop=0, optimize=False, disposal=2)
    print(f"wrote {args.output} ({len(quantized)} frames, {args.output.stat().st_size / 1e6:.1f} MB)")
    return 0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Zoom into the photon ring and write it out as a GIF.")
    parser.add_argument("-o", "--output", type=Path, default=Path("web_images/photon_ring_zoom.gif"))
    parser.add_argument("--frames", type=int, default=48, help="frames in the zoom (default 48)")
    parser.add_argument("--max-zoom", type=float, default=USEFUL_MAX_ZOOM,
                        help=f"magnification of the final frame (default {USEFUL_MAX_ZOOM:g}; past this "
                             "the integrator's tolerance floor, not the zoom, is the limit)")
    parser.add_argument("--width", type=int, default=640, help="frame width; height is half (default 640)")
    parser.add_argument("--errormax", type=float, default=TOLERANCE_FLOOR,
                        help=f"integrator tolerance (default {TOLERANCE_FLOOR:g}, the tightest "
                             "dopri54_step responds to); also sets the step budget, int(1/errormax)")
    parser.add_argument("--de0", type=float, default=scene.SETTLE_ACCURACY["de0"],
                        help="integrator step ceiling (default matches web_server's settle tier)")
    parser.add_argument("--cx", type=float, default=0.5,
                        help="column of the image plane to zoom down, 0..1 (default 0.5)")
    parser.add_argument("--edge-index", type=int, default=0,
                        help="which crossing of the shadow's edge down that column to zoom "
                             "into, top first (default 0). A column through the hole meets it "
                             "twice; --list-edges shows what was found")
    parser.add_argument("--list-edges", action="store_true",
                        help="scan the column, report every shadow edge on it, and stop")
    parser.add_argument("--target-cy", type=float, default=None,
                        help="zoom into this row instead of onto a shadow edge, skipping the "
                             "search entirely - for any other feature worth magnifying. Note the "
                             "subring cascade lives at the shadow's edge, so elsewhere the zoom "
                             "runs out of new detail much sooner")
    parser.add_argument("--scan-samples", type=int, default=96,
                        help="rows sampled when locating the edges (default 96)")
    camera = parser.add_argument_group(
        "scene", "Defaults match web_server's opening view. a and Q are in the same code "
                 "units as the web UI's sliders, where the extremal bound is rs/2 = "
                 f"{scene.RS / 2:g}; they are scaled down together if they exceed it.")
    camera.add_argument("--a", type=float, default=0.0, help="black hole spin")
    camera.add_argument("--Q", type=float, default=0.0, help="black hole charge")
    camera.add_argument("--r0", type=float, default=scene.state["r0"], help="camera radius")
    camera.add_argument("--theta0", type=float, default=scene.state["theta0"], help="camera polar angle")
    camera.add_argument("--phi0", type=float, default=scene.state["phi0"], help="camera azimuth")
    # ~45 halvings take a 0.3-wide bracket below the ROI grid's own quantum
    # (1/2^47 of the frame height), past which the answer is quantized away
    # anyway; the extra few are slack, and each one costs a single ray.
    parser.add_argument("--bisection-steps", type=int, default=48)
    parser.add_argument("--ms-per-frame", type=int, default=110)
    parser.add_argument("--pingpong", action="store_true",
                        help="play the zoom back out again, so the loop is seamless")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.max_zoom > scene.MAX_ROI_ZOOM:
        print(f"--max-zoom clamped to the renderer's limit, {scene.MAX_ROI_ZOOM:.0f}x")
        args.max_zoom = scene.MAX_ROI_ZOOM
    if args.frames < 2:
        raise SystemExit("--frames must be at least 2")

    resolution = (args.width, args.width // 2)
    accuracy = {"errormax": args.errormax, "de0": args.de0}
    # Same clamp the interactive UI applies, so a and Q can never be pushed past
    # the extremal bound where the horizon stops existing.
    spin, charge = scene.clamp_spin_charge(args.a, args.Q)
    camera = {"r0": args.r0, "theta0": args.theta0, "phi0": args.phi0, "a": spin, "Q": charge}
    physics = scene.disk_physics(spin, charge)
    period = photon_orbit_period(scene.RS)
    worker = render_worker.INTERACTIVE_WORKER

    mass = scene.RS / 2.0
    if (spin, charge) != (args.a, args.Q):
        print(f"a and Q scaled to the extremal bound: a={spin:.6g}, Q={charge:.6g}")
    print(f"scene: a={spin:.6g} (a/M={spin / mass:.3f})  Q={charge:.6g} (Q/M={charge / mass:.3f})  "
          f"camera r0={camera['r0']:g} theta0={camera['theta0']:g} phi0={camera['phi0']:g}")

    # cuda_ray.h gives every ray a budget of int(1/errormax) steps and records
    # one that runs out as having escaped rather than captured. So errormax
    # governs how many orbits can be shown at all, not merely how accurately.
    print(f"errormax {args.errormax:g} -> {int(1.0 / args.errormax):,} steps per ray, "
          f"step tolerance {max(args.errormax, TOLERANCE_FLOOR):g}")
    if args.errormax > TOLERANCE_FLOOR * 10:
        print("  note: loose enough that rays will be cut off mid-orbit and drawn as sky;")
        print(f"        the subrings will be faint or missing. {TOLERANCE_FLOOR:g} is the useful floor.")
    if args.max_zoom > USEFUL_MAX_ZOOM:
        print(f"  note: past ~{USEFUL_MAX_ZOOM:g}x the shadow edge being zoomed into is where rays run")
        print("        out of budget, not the true critical curve, and the lensed disk leaves the")
        print("        frame entirely. Frames beyond that point show shadow and sky only.")

    started = time.time()
    if args.target_cy is not None:
        edge_cy = args.target_cy
        print(f"zooming into the given point cx={args.cx}, cy={edge_cy!r} (no edge search)")
        return render_movie(worker, camera, args, edge_cy, accuracy, resolution, physics, period)

    print(f"scanning column cx={args.cx} for the shadow's edge ...")
    brackets = scan_column(worker, camera, args.cx, accuracy, args.scan_samples)
    for index, (upper, lower, captured_below) in enumerate(brackets):
        entering = "into shadow" if captured_below else "out of shadow"
        print(f"  edge {index}: between cy {upper:.4f} and {lower:.4f} ({entering})")
    if not brackets:
        print(f"  none - column cx={args.cx} misses the shadow. Looking for ones that do not ...")
        low, high = columns_crossing_shadow(worker, camera, accuracy)
        if low is None:
            raise SystemExit("no column meets the shadow at all: this camera is not looking at it.")
        raise SystemExit(
            f"no shadow edge on column cx={args.cx}. Columns roughly {low:.3f} to {high:.3f} do "
            f"cross it - note that is narrower than the shadow looks, because cx is a fraction of "
            f"a frame twice as wide as it is tall.")
    if args.list_edges:
        return 0
    if not 0 <= args.edge_index < len(brackets):
        raise SystemExit(f"--edge-index {args.edge_index} but this column has "
                         f"{len(brackets)} edge(s), numbered 0..{len(brackets) - 1}")

    edge_cy = find_shadow_edge(worker, camera, args.cx, brackets[args.edge_index],
                               accuracy, args.bisection_steps)
    print(f"  zooming into edge {args.edge_index} at cy = {edge_cy!r}  "
          f"({time.time() - started:.1f}s, {args.scan_samples} scan + "
          f"{args.bisection_steps} bisection rays)")

    return render_movie(worker, camera, args, edge_cy, accuracy, resolution, physics, period)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    finally:
        render_worker.shutdown_all()
