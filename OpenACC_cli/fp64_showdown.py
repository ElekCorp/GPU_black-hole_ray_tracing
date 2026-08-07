#!/usr/bin/env python3
"""Render fp32 and fp64 side by side, zooming in until the fp32 trace falls apart.

WHAT THIS IS SHOWING

The two panels are the same camera, the same zoom window, the same integrator
tolerance and the same shading. The only difference is the working precision of
the geodesic integration: render_frame_roi_f32 on the left,
render_frame_roi_f64 on the right (ffi_bridge.cpp). At the start they are the
same picture. They stop being the same picture, and the point of the movie is
where that happens and what it looks like.

Why it happens is not the integrator's tolerance - tightening errormax does not
move it. There are two separate failures, they look different on screen, and the
movie reaches them in this order:

  1. CHAOTIC AMPLIFICATION, from around a thousand x, and this is the main
     event. A ray that skims the photon sphere is exponentially sensitive to
     where it started: that is the same instability that produces the subring
     cascade in the first place. It amplifies fp32's rounding just as
     efficiently as it amplifies a real difference in aim, so a ray ends up on
     the far side of the hole from where it belonged. The rays are all still
     distinct and still finite - this is not a resolution failure - they are
     individually wrong, and the frame goes to salt-and-pepper noise. Measured
     at 240px panels: agreement with fp64 falls to 84% at 1,350x, 64% at 2,462x
     and 40% at 8,185x, with zero collapsed pairs throughout.

  2. RAYS GOING NON-FINITE, from around 1e5 x. Deeper in, integrations start
     failing outright: the state goes NaN or inf, ray_step_T takes its bail-out,
     and sky_direction - which cannot name a direction for a state that is not
     finite - substitutes one fixed direction, (asin(1), 0). Those pixels are
     not traced results at all, and because they all carry the same sentinel
     they read on screen as flat patches. fp32 reaches 13.8% of the frame at
     1e6x against fp64's 7.1%.

WHAT IS NOT HAPPENING, THOUGH IT IS THE OBVIOUS GUESS

The textbook failure - two adjacent pixels' launch directions rounding to the
same float, so one geodesic gets traced where the frame asked for two - does not
meaningfully occur anywhere in this range, in either precision. Measured
directly, excluding the non-finite sentinel above: over 24 frames from 1x to
1e6x at 240px panels, fp32 collapsed exactly zero neighbouring pairs, and fp64
collapsed three (out of 12,045) in a single frame at 1.65e5x. Both are nothing.

This is worth stating because it is easy to measure wrong. Counting the sentinel
pixels as collapse makes fp32 look like it is losing resolution - the number
climbs smoothly to 26% at 1e6x - and it is not; those are broken rays. The
script therefore reports the two separately, and ijk_to_vec_mink_zoom's decision
to compute the launch direction in double even for the f32 trace is a large part
of why genuine collapse stays away.

FOUR MEASURED NUMBERS, NOT FOUR ASSERTED ONES

Each frame is captioned with quantities read off the two traces:

  agree      fraction of pixels where the two precisions put the ray in the same
             place: through the horizon, onto the disk, or out to the sky. This
             is the one that tracks what you can see, and it catches both
             failures above.

  non-finite fraction of the frame that came back as sky_direction's fallback
             sentinel, i.e. rays whose integration blew up. Failure 2, isolated.

  collapse   fraction of horizontally adjacent pixel pairs that came back
             bit-identical, among pairs that landed on the same kind of target
             and excluding the sentinel above. Genuine ray collapse, isolated -
             and it stays at zero, which is the point.

  pixels     mean absolute difference between the two finished 8-bit images.
             Both panels are shaded with the same seed and the same frame index,
             so their film grain is identical and this is signal, not noise.

The zoom target is the shadow's edge, bisected for in fp64 rather than assumed,
because that is where the subring cascade lives and therefore where the
precision runs out first.

HOW FAR TO GO

Past roughly 1e6 magnification fp64 has its own limit, but a different one:
dopri54_step clamps its step tolerance to fmax(errormax, 1e-8), so the *step*
budget rather than the number format becomes the binding constraint, and the
frame stops containing lensed disk at all. So the honest claim this movie
supports is bounded: fp64 buys roughly two more decades of usable zoom, after
which the integrator, not the representation, is what needs fixing. The default
--max-zoom stops inside that range.

Usage:
    python fp64_showdown.py                              # 36 frames to 1e6x
    python fp64_showdown.py --frames 60 --width 400 -o deep.gif
    python fp64_showdown.py --max-zoom 1e8 --csv stats.csv
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np
from PIL import Image
from tqdm import tqdm

import movie_common as mc
import render_worker
import web_server as scene


# dopri54_step clamps the tolerance it actually uses to fmax(errormax, 1e-8), so
# nothing below this changes a trajectory; it only buys step budget, which is
# int(1/errormax) and is not clamped.
TOLERANCE_FLOOR = 1e-8

# Where fp64's own ceiling starts to bind, for a different reason than fp32's -
# see the docstring. Past this the movie stops being about precision.
USEFUL_MAX_ZOOM = 1e6

GUTTER = 6  # pixels between the two panels; labels are drawn inside them

# The footer strip is reserved before the text is drawn, so the caption font has
# to be capped or a large render sizes its text past the strip and clips the last
# line. Both numbers come from the same place, rather than one being a guess at
# the other.
FOOTER_FONT = 12


def hit_class(frame) -> np.ndarray:
    """0 captured, 1 escaped, 2 disk - what the ray hit, ignoring how brightly."""
    out = np.full(frame.radius.shape, 2, dtype=np.int8)
    out[frame.radius == 0.0] = 0
    out[frame.radius < 0.0] = 1
    return out


def broken_mask(frame) -> np.ndarray:
    """Escaped pixels whose ray state went non-finite during integration.

    cuda_ray.h's sky_direction cannot name a direction for a ray whose state has
    stopped being finite, so rather than propagate a NaN into the shading layer
    it returns one fixed direction, (asin(1), 0). ray_step_T also routes its
    explicit isnan/isinf bail-outs through that same call. Every such pixel
    therefore carries a bit-identical sentinel rather than a traced direction,
    and they cluster exactly where this movie spends its time.

    A ray that genuinely escapes along (pi/2, 0) is a measure-zero coincidence
    (it needs dz and dy to vanish exactly), so this identifies the sentinel
    without meaningfully over-counting.
    """
    escaped = frame.radius < 0.0
    fallback_theta = frame.sky_theta.dtype.type(math.asin(1.0))
    fallback_phi = frame.sky_theta.dtype.type(0.0)
    return escaped & (frame.sky_theta == fallback_theta) & (frame.sky_phi == fallback_phi)


def broken_fraction(frame) -> float:
    """Fraction of the whole frame that came back as the non-finite sentinel."""
    return float(broken_mask(frame).mean())


def collapse_fraction(frame) -> tuple:
    """Fraction of adjacent pairs that traced a genuinely bitwise-identical ray.

    Only pairs whose two pixels hit the same kind of target are counted, and the
    comparison uses whichever channels carry continuous information for that
    target: the escape direction for sky, the landing radius and azimuth for the
    disk. Two exclusions, both of which otherwise turn this into a measurement of
    something else entirely:

      - captured rays, because every ray through the horizon reports the same
        zeros, so counting them would report the renderer's sentinel as collapse;
      - the non-finite fallback above, for exactly the same reason. This one is
        not hypothetical: counting it made this function report "collapse" rising
        to 26% at 1e6x, when excluding it gives 0.00% at every zoom in both
        precisions. Those pixels are broken rays, not two pixels sharing one ray,
        and conflating them credits fp64 with fixing a failure that is not the
        one it fixes.

    Returns (fraction, pairs_counted); the fraction is NaN when nothing qualified.
    """
    escaped = frame.radius < 0.0
    disk = frame.radius > 0.0
    usable = escaped & ~broken_mask(frame)

    sky_pair = usable[:, :-1] & usable[:, 1:]
    disk_pair = disk[:, :-1] & disk[:, 1:]

    same = np.zeros(sky_pair.shape, dtype=bool)
    counted = sky_pair | disk_pair
    if sky_pair.any():
        identical = ((np.diff(frame.sky_theta, axis=1) == 0.0)
                     & (np.diff(frame.sky_phi, axis=1) == 0.0))
        same |= sky_pair & identical
    if disk_pair.any():
        identical = ((np.diff(frame.radius, axis=1) == 0.0)
                     & (np.diff(frame.phi, axis=1) == 0.0))
        same |= disk_pair & identical

    total = int(counted.sum())
    if total == 0:
        return float("nan"), 0
    return float(same[counted].mean()), total


def panel(image: Image.Image, label: str) -> Image.Image:
    canvas = image.convert("RGB")
    mc.draw_caption(canvas, [label], origin=(6, 6))
    return canvas


def stitch(left: Image.Image, right: Image.Image, footer: list) -> Image.Image:
    """Two panels side by side, with a shared caption block under them."""
    width = left.width + GUTTER + right.width
    footer_height = (0 if not footer
                     else mc.caption_line_height(FOOTER_FONT) * len(footer) + 8)
    canvas = Image.new("RGB", (width, left.height + footer_height), (0, 0, 0))
    canvas.paste(left, (0, 0))
    canvas.paste(right, (left.width + GUTTER, 0))
    if footer:
        mc.draw_caption(canvas, footer, origin=(6, left.height + 4), width=width,
                        max_size=FOOTER_FONT)
    return canvas


def render_movie(worker, camera: dict, physics, args, target_cy: float) -> int:
    accuracy = {"errormax": args.errormax, "de0": args.de0}
    resolution = (args.width, args.width // 2)

    images = []
    rows = []
    for index in tqdm(range(args.frames), desc="fp32 vs fp64", unit="frame"):
        fraction = index / max(args.frames - 1, 1)
        zoom = args.max_zoom**fraction  # geometric: every frame gains the same factor
        roi = scene.clamp_roi(args.cx, target_cy, 1.0 / zoom)

        started = time.time()
        # Deliberately NOT scaled by scene.zoom_accuracy: both panels must run at
        # one fixed tolerance, or a difference between them could be blamed on
        # the integrator having been given different work rather than on the
        # precision, which is the whole subject.
        f32 = mc.trace(worker, camera, roi, accuracy, resolution, precision="f32")
        f64 = mc.trace(worker, camera, roi, accuracy, resolution, precision="f64")

        agree = float((hit_class(f32) == hit_class(f64)).mean())
        broken32 = broken_fraction(f32)
        broken64 = broken_fraction(f64)
        collapse32, pairs32 = collapse_fraction(f32)
        collapse64, pairs64 = collapse_fraction(f64)

        # Same frame index for both, so the grain finish_frame adds is the same
        # pattern in both panels and the pixel difference below is real.
        left_image = mc.shade(f32, physics, 0.0, index)
        right_image = mc.shade(f64, physics, 0.0, index)
        pixel_gap = float(np.abs(np.asarray(left_image, dtype=np.int16)
                                 - np.asarray(right_image, dtype=np.int16)).mean())

        footer = [
            f"{mc.format_zoom(zoom)}    fp32 agrees with fp64 on {agree * 100:5.1f}% of pixels"
            f"    mean difference {pixel_gap:.1f}/255",
            f"rays that went non-finite:  fp32 {broken32 * 100:5.1f}%   fp64 {broken64 * 100:5.1f}%"
            f"    collapsed neighbours: {collapse32 * 100:.1f}% / {collapse64 * 100:.1f}%",
        ]
        images.append(stitch(panel(left_image, "fp32"), panel(right_image, "fp64"), footer))
        rows.append((zoom, agree, broken32, broken64, collapse32, collapse64,
                     pixel_gap, pairs32, pairs64))

        tqdm.write(f"  {index + 1:3d}/{args.frames} {mc.format_zoom(zoom):>10}  "
                   f"agree {agree:6.3f}  broken {broken32:6.3f}/{broken64:6.3f}  "
                   f"collapse {collapse32:6.3f}/{collapse64:6.3f}  "
                   f"dpix {pixel_gap:5.2f}  [{time.time() - started:5.1f}s]")

    mc.write_gif(images, args.output, args.ms_per_frame, args.pingpong)

    if args.csv:
        args.csv.parent.mkdir(parents=True, exist_ok=True)
        with args.csv.open("w") as handle:
            handle.write("zoom,agree,broken_f32,broken_f64,collapse_f32,collapse_f64,"
                         "mean_pixel_diff,pairs_f32,pairs_f64\n")
            for row in rows:
                handle.write(",".join(f"{value:.6g}" for value in row) + "\n")
        print(f"wrote {args.csv}")

    # The headline, stated from this run's own numbers rather than from the
    # docstring, and separating the failures because they are separate.
    print("\nwhat this run measured:")
    noisy = [row for row in rows if row[1] < 0.90]
    if noisy:
        print(f"  fp32 diverged from fp64 (agreement < 90%) at {mc.format_zoom(noisy[0][0])}, "
              f"with {noisy[0][4] * 100:.1f}% of neighbouring pairs collapsed and "
              f"{noisy[0][2] * 100:.1f}% of rays non-finite -")
        print("    i.e. the rays were still distinct and finite, and individually wrong: chaotic")
        print("    amplification of rounding, not lost resolution.")
    else:
        print(f"  fp32 stayed within 90% agreement all the way to {mc.format_zoom(rows[-1][0])}")

    for label, index in (("fp32", 2), ("fp64", 3)):
        broke = [row for row in rows if row[index] > 0.01]
        if broke:
            print(f"  {label} rays began going non-finite at {mc.format_zoom(broke[0][0])} "
                  f"({broke[0][index] * 100:.1f}%), reaching {rows[-1][index] * 100:.1f}% by "
                  f"{mc.format_zoom(rows[-1][0])}")
        else:
            print(f"  {label} produced no non-finite rays up to {mc.format_zoom(rows[-1][0])}")

    worst32 = max(row[4] for row in rows)
    worst64 = max(row[5] for row in rows)
    if max(worst32, worst64) <= 0.0:
        print(f"  neither precision collapsed a single neighbouring pair anywhere up to "
              f"{mc.format_zoom(rows[-1][0])}:")
        print("    at these resolutions the limit is the integrator losing rays, not the "
              "representation running out of ulps.")
    else:
        print(f"  peak genuine ray collapse: fp32 {worst32 * 100:.2f}%, fp64 {worst64 * 100:.2f}%")
    return 0


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Zoom into the photon ring in fp32 and fp64 side by side.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--output", type=Path, default=Path("web_images/fp64_showdown.gif"))
    parser.add_argument("--frames", type=int, default=36)
    parser.add_argument("--width", type=int, default=320,
                        help="width of EACH panel; the GIF is a little over twice this")
    parser.add_argument("--max-zoom", type=float, default=USEFUL_MAX_ZOOM,
                        help="magnification of the last frame")
    parser.add_argument("--ms-per-frame", type=int, default=140)
    parser.add_argument("--pingpong", action="store_true",
                        help="play the zoom back out again, so the loop is seamless")
    parser.add_argument("--csv", type=Path, default=None,
                        help="also write the per-frame measurements here")

    target = parser.add_argument_group("target")
    target.add_argument("--cx", type=float, default=0.5,
                        help="column of the image plane to zoom down, 0..1")
    target.add_argument("--edge-index", type=int, default=0,
                        help="which shadow-edge crossing down that column to zoom into, top first")
    target.add_argument("--target-cy", type=float, default=None,
                        help="zoom into this row instead, skipping the edge search entirely")
    target.add_argument("--scan-samples", type=int, default=64,
                        help="rows sampled when locating the edges")
    target.add_argument("--bisection-steps", type=int, default=48)

    accuracy = parser.add_argument_group("accuracy")
    accuracy.add_argument("--errormax", type=float, default=TOLERANCE_FLOOR,
                          help="integrator tolerance, held fixed across both panels and every "
                               "frame; also sets the per-ray step budget, int(1/errormax)")
    accuracy.add_argument("--de0", type=float, default=scene.SETTLE_ACCURACY["de0"],
                          help="integrator step ceiling")

    mc.add_scene_args(parser)
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.frames < 2:
        raise SystemExit("--frames must be at least 2")
    if args.max_zoom > scene.MAX_ROI_ZOOM:
        print(f"--max-zoom clamped to the renderer's limit, {scene.MAX_ROI_ZOOM:.0f}x")
        args.max_zoom = scene.MAX_ROI_ZOOM
    if args.max_zoom > USEFUL_MAX_ZOOM:
        print(f"note: past ~{USEFUL_MAX_ZOOM:g}x fp64 hits the integrator's step-tolerance floor")
        print("      rather than any limit of the representation, and the lensed disk leaves the")
        print("      frame. Those frames no longer illustrate precision.")

    camera, physics = mc.resolve_scene(args)
    accuracy = {"errormax": args.errormax, "de0": args.de0}
    worker = render_worker.INTERACTIVE_WORKER

    if args.target_cy is not None:
        target_cy = args.target_cy
        print(f"zooming into the given point cx={args.cx}, cy={target_cy!r} (no edge search)")
    else:
        started = time.time()
        print(f"scanning column cx={args.cx} for the shadow's edge (fp64 probes) ...")
        brackets = mc.scan_column(worker, camera, args.cx, accuracy, args.scan_samples)
        for index, (upper, lower, below) in enumerate(brackets):
            print(f"  edge {index}: between cy {upper:.4f} and {lower:.4f} "
                  f"({'into shadow' if below else 'out of shadow'})")
        if not brackets:
            raise SystemExit(f"no shadow edge on column cx={args.cx}: this column misses the "
                             "shadow. Note cx is a fraction of a frame twice as wide as it is "
                             "tall, so the shadow spans a narrower range of cx than it looks.")
        if not 0 <= args.edge_index < len(brackets):
            raise SystemExit(f"--edge-index {args.edge_index} but this column has "
                             f"{len(brackets)} edge(s), numbered 0..{len(brackets) - 1}")
        target_cy = mc.find_shadow_edge(worker, camera, args.cx, brackets[args.edge_index],
                                        accuracy, args.bisection_steps)
        print(f"  zooming into edge {args.edge_index} at cy = {target_cy!r} "
              f"({time.time() - started:.1f}s)")

    started = time.time()
    result = render_movie(worker, camera, physics, args, target_cy)
    print(f"total {time.time() - started:.1f}s")
    return result


if __name__ == "__main__":
    try:
        sys.exit(main())
    finally:
        render_worker.shutdown_all()
