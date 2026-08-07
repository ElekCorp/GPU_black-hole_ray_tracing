#!/usr/bin/env python3
"""Render a looping GIF of the camera moving around the hole, optionally pushing in.

FOUR MOVES, AND WHY THEY ARE NOT EQUALLY INTERESTING

  inclination  the camera swings from looking down on the disk to seeing it
               edge-on. This is the one worth rendering. Face-on, the disk is a
               broad annulus with the shadow punched through it; as the camera
               drops into the disk's plane the far side of the disk is lensed up
               and over the shadow and the whole thing closes into the familiar
               band with a bright ring around the silhouette. Nothing about the
               disk changes - only which geodesics reach the camera.

  azimuth      the camera orbits the spin axis. Be aware of what this does and
               does not show: Kerr-Newman is axisymmetric and the disk is an
               axisymmetric annulus, so the geometry of the image is *invariant*
               under this move. The shadow, the ring and the lensed disk are
               pixel-for-pixel the same at every phi0. What rotates is the disk's
               turbulence texture, which is a function of the traced azimuth.
               So this is a pretty loop, not a new view - and it is exactly the
               symmetry that makes phi0 free to precompute (see
               docs/precompute-investigation.md).

  look         the camera stays put and swings its gaze in a cone, so the hole
               travels around the frame. This composes a pan and a tilt at once,
               which is the case cuda_ray.h's Rodrigues rotation used to get
               wrong (see tests/test_camera_rotation.cpp) - a pan alone, or a
               tilt alone, both happened to land in the degenerate set where the
               old expression was accidentally correct.

  dolly        the camera moves radially. Only inwards is available: ffi_bridge.cpp
               fixes the outer integration sphere at sugar_ki = 1.01, and a camera
               outside it has every ray terminate on its first step, so the frame
               comes back as an empty starfield. Measured: r0 = 1.01 renders,
               r0 = 1.02 is blank.

ZOOM ON TOP OF THE MOVE

--zoom-to narrows the field of view over the course of the movie while the move
runs, which is a true optical zoom rather than a crop - the same pixel count is
retraced across a narrower solid angle. The integrator tolerance is tightened as
the window narrows, following the same policy the interactive UI uses, and the
trace switches to fp64 once the window is narrower than f32 can resolve.

The window centre (--cx, --cy) is FIXED for the whole movie, and there is no
tracking. That is fine for a zoom on its own, and fine alongside `look`, whose
whole job is to move the frame. It is a trap alongside `inclination` and
`dolly`, which move the subject: a window aimed at the ring in the first frame
is aimed at empty disk, or at the inside of the shadow, by the last one. Past
about 30x the two fight each other badly enough that the end of the movie is a
flat single-colour frame. The script reports how many frames came back with no
structure in them, so this is visible from the log rather than only from the
GIF; if you want both moves, zoom gently or render them as separate shots.

LOOPING

Frozen disk flow is the default so every move loops seamlessly: the cyclic moves
(azimuth, look) return to their starting frame exactly, and the sweeps
(inclination, dolly) are played back out again. Turning the flow on with
--flow-seconds animates the accretion disk, at the cost of a visible seam where
the loop closes, because the turbulence field does not have the same period as
the camera move.

Usage:
    python orbit_flyaround.py                                   # inclination sweep
    python orbit_flyaround.py --mode look --frames 60
    python orbit_flyaround.py --mode inclination --zoom-to 200 -o push.gif
    python orbit_flyaround.py --a 0.02 --mode azimuth --flow-seconds 8
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

from tqdm import tqdm

import movie_common as mc
import render_worker
import web_server as scene


MODES = ("inclination", "azimuth", "look", "dolly")

# Cyclic moves come back to their own first frame, so replaying them backwards
# would just render the loop twice. Sweeps do not, so they ping-pong.
CYCLIC = ("azimuth", "look")

# How far the gaze swings in `look` mode. Kept modest on purpose: the frustum is
# only 67 degrees wide, so a larger cone walks the hole clean out of shot.
LOOK_CONE_RADIANS = 0.22


def move_state(mode: str, fraction: float, base_camera: dict, args) -> tuple:
    """The camera and orientation at `fraction` (0..1) through the move.

    Returns (camera, omega, label) where label is the caption's description of
    where in the move this frame sits.
    """
    camera = dict(base_camera)
    omega = scene.current_omega(scene.BASE_ORIENTATION)

    if mode == "inclination":
        # Eased, so the sweep starts and ends at rest rather than jerking into
        # the turn at the ends of the ping-pong.
        theta = args.theta_from + (args.theta_to - args.theta_from) * mc.ease(fraction)
        camera["theta0"] = mc.clamp_polar(theta)
        label = f"inclination {math.degrees(camera['theta0']):.0f} deg"

    elif mode == "azimuth":
        camera["phi0"] = base_camera["phi0"] + 2.0 * math.pi * fraction
        label = f"azimuth {math.degrees(2.0 * math.pi * fraction):.0f} deg"

    elif mode == "look":
        # Pan and tilt together, tracing a circle: this is the compound rotation
        # the fixed Rodrigues matrix is needed for.
        angle = 2.0 * math.pi * fraction
        orientation = scene.apply_look_delta(
            scene.BASE_ORIENTATION,
            LOOK_CONE_RADIANS * math.cos(angle),
            LOOK_CONE_RADIANS * math.sin(angle),
            0.0,
        )
        omega = scene.current_omega(orientation)
        label = f"gaze {math.degrees(angle):.0f} deg around a {math.degrees(LOOK_CONE_RADIANS):.0f} deg cone"

    elif mode == "dolly":
        camera["r0"] = args.r0_from + (args.r0_to - args.r0_from) * mc.ease(fraction)
        label = f"r0 = {camera['r0']:.3f}"

    else:
        raise SystemExit(f"unknown --mode {mode!r}; choose from {', '.join(MODES)}")

    return camera, omega, label


def render_movie(worker, base_camera: dict, physics, args) -> int:
    accuracy = {"errormax": args.errormax, "de0": args.de0}
    resolution = (args.width, args.width // 2)
    pingpong = args.pingpong if args.pingpong is not None else args.mode not in CYCLIC

    # A cyclic move must not render both endpoints of the cycle, or the loop
    # stutters on a duplicated frame; a ping-ponged sweep must render both, or it
    # never reaches the far end.
    divisor = args.frames if args.mode in CYCLIC else max(args.frames - 1, 1)

    images = []
    empty_frames = 0
    for index in tqdm(range(args.frames), desc=args.mode, unit="frame"):
        fraction = index / divisor
        camera, omega, label = move_state(args.mode, fraction, base_camera, args)

        # Geometric in the zoom, so each frame gains the same factor - a linear
        # ramp would spend most of the movie almost unzoomed and then lurch.
        zoom = args.zoom_to**fraction if args.zoom_to > 1.0 else 1.0
        roi = scene.clamp_roi(args.cx, args.cy, 1.0 / zoom)
        frame_accuracy = scene.zoom_accuracy(accuracy, roi)
        precision = scene.effective_precision(roi, args.precision)

        started = time.time()
        frame = mc.trace(worker, camera, roi, frame_accuracy, resolution, omega, precision)
        escaped, disk, captured = mc.hit_fractions(frame)
        # A frame where one outcome covers essentially everything has no
        # structure in it: all sky (camera looking away, or outside sugar_ki),
        # all shadow, or all disk. That is the usual result of a fixed zoom
        # window plus a move that carries the subject out of it.
        if max(escaped, disk, captured) > 0.995:
            empty_frames += 1

        flow_time = args.flow_seconds * fraction
        image = mc.shade(frame, physics, flow_time, index)
        if not args.no_caption:
            lines = [label]
            if zoom > 1.0:
                lines.append(f"field of view narrowed {mc.format_zoom(zoom)}"
                             + (" (fp64)" if precision == "f64" else ""))
            mc.draw_caption(image, lines)
        images.append(image)

        tqdm.write(f"  {index + 1:3d}/{args.frames}  {label:<34}  {precision}  "
                   f"sky {escaped:.2f} disk {disk:.2f} shadow {captured:.2f}  "
                   f"{time.time() - started:5.1f}s")

    if empty_frames:
        print(f"warning: {empty_frames} of {args.frames} frame(s) came back featureless - one "
              "outcome covered >99.5% of the frame.")
        if args.zoom_to > 1.0:
            print("         With --zoom-to this is usually the fixed zoom window losing the "
                  "subject as the move carries it away;")
            print("         zoom less far, or retarget with --cx/--cy, or render the zoom as its "
                  "own shot.")

    mc.write_gif(images, args.output, args.ms_per_frame, pingpong)
    return 0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Move the camera around the hole and write it out as a looping GIF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--output", type=Path, default=Path("web_images/flyaround.gif"))
    parser.add_argument("--mode", choices=MODES, default="inclination",
                        help="which camera move to render; see the module docstring for what "
                             "each one does and does not show")
    parser.add_argument("--frames", type=int, default=48)
    parser.add_argument("--width", type=int, default=480, help="frame width; height is half")
    parser.add_argument("--ms-per-frame", type=int, default=80)
    parser.add_argument("--pingpong", action="store_true", default=None,
                        help="play the move back out again; None means decide per mode - on for "
                             "inclination and dolly, off for the cyclic moves which already "
                             "close their own loop")
    parser.add_argument("--no-pingpong", dest="pingpong", action="store_false")
    parser.add_argument("--no-caption", action="store_true")

    move = parser.add_argument_group("move")
    move.add_argument("--theta-from", type=float, default=0.32,
                      help="inclination mode: starting polar angle, radians (near face-on)")
    move.add_argument("--theta-to", type=float, default=math.pi / 2,
                      help="inclination mode: ending polar angle (pi/2 is exactly edge-on)")
    move.add_argument("--r0-from", type=float, default=1.0,
                      help=f"dolly mode: starting radius (must be < sugar_ki = {mc.SUGAR_KI:g})")
    move.add_argument("--r0-to", type=float, default=0.45, help="dolly mode: ending radius")
    move.add_argument("--flow-seconds", type=float, default=0.0,
                      help="seconds of accretion-disk flow across the movie; 0 freezes it, "
                           "which is what keeps the loop seamless")

    zoom = parser.add_argument_group("zoom")
    zoom.add_argument("--zoom-to", type=float, default=1.0,
                      help="narrow the field of view to this magnification by the last frame "
                           "(1 = no zoom). A true optical zoom: the same rays are retraced "
                           "across a narrower solid angle")
    zoom.add_argument("--cx", type=float, default=0.5, help="zoom window centre, column 0..1")
    zoom.add_argument("--cy", type=float, default=0.5, help="zoom window centre, row 0..1")
    zoom.add_argument("--precision", choices=scene.PRECISION_MODES, default="auto",
                      help="'auto' follows the zoom (f32 until it can no longer resolve "
                           "neighbouring rays, then f64); f32/f64 pin it")

    accuracy = parser.add_argument_group("accuracy")
    accuracy.add_argument("--errormax", type=float, default=scene.SETTLE_ACCURACY["errormax"],
                          help="integrator tolerance; also sets the per-ray step budget, "
                               "int(1/errormax)")
    accuracy.add_argument("--de0", type=float, default=scene.SETTLE_ACCURACY["de0"],
                          help="integrator step ceiling")

    mc.add_scene_args(parser)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.frames < 2:
        raise SystemExit("--frames must be at least 2")
    if args.zoom_to > scene.MAX_ROI_ZOOM:
        print(f"--zoom-to clamped to the renderer's limit, {scene.MAX_ROI_ZOOM:.0f}x")
        args.zoom_to = scene.MAX_ROI_ZOOM

    # dolly overrides r0, so validate the radii it will actually visit rather
    # than only the static --r0.
    if args.mode == "dolly":
        for name, value in (("--r0-from", args.r0_from), ("--r0-to", args.r0_to)):
            if value >= mc.SUGAR_KI:
                raise SystemExit(
                    f"{name} {value:g} is outside the outer integration sphere "
                    f"(sugar_ki = {mc.SUGAR_KI:g}); those frames would render empty.")
        args.r0 = args.r0_from

    camera, physics = mc.resolve_scene(args)

    if args.mode == "azimuth":
        print("note: the metric and the disk are both axisymmetric, so an azimuth orbit leaves "
              "the image geometry unchanged -")
        print("      only the disk's turbulence texture rotates. Use --mode inclination for a "
              "genuinely new view.")

    # A cyclic move closes its own loop, but a zoom ramped across it does not:
    # the last frame is near full magnification and the first is unzoomed, so
    # the GIF cuts rather than loops.
    if args.zoom_to > 1.0 and args.mode in CYCLIC and not args.pingpong:
        print(f"note: --zoom-to with --mode {args.mode} will not loop cleanly - the move returns "
              "to its start but")
        print("      the zoom does not. Add --pingpong to ramp the zoom back out again.")

    started = time.time()
    result = render_movie(worker=render_worker.INTERACTIVE_WORKER,
                          base_camera=camera, physics=physics, args=args)
    print(f"total {time.time() - started:.1f}s")
    return result


if __name__ == "__main__":
    try:
        sys.exit(main())
    finally:
        render_worker.shutdown_all()
