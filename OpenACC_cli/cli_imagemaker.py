#!/usr/bin/env python3
"""Turn the ray tracer's scalar hit buffer into a cinematic black-hole image.

The C++ renderer remains the source of truth for ray paths.  Its buffer uses
``-1`` for rays that escape, ``0`` for captured rays, and a positive radius for
an accretion-disk hit; alongside that it reports the ray-traced frequency
shift, the disk azimuth, the light travel time, and the direction each escaping
ray was travelling when it left (see cuda_ray.h's RAY_CHANNELS).  Nothing here
re-derives a trajectory: this script only shades what the geodesics already
determined, applying a thermal disk model, the gravitational/Doppler shift, a
lensed sky sampled along the traced escape directions, and a film finish.
"""

from __future__ import annotations

import argparse
import functools
import math
from pathlib import Path
from typing import NamedTuple

import numpy as np
from PIL import Image


DEFAULT_INPUT = Path("./web_images/kep_cli.dat")
DEFAULT_OUTPUT = Path("./web_images/blackhole_cli.png")

RAY_CHANNELS = 6  # keep in step with cuda_ray.h


class RayFrame(NamedTuple):
    """One frame of ray-traced scene data, image-shaped (height, width).

    Straight from the geodesic integrator, in its own units - nothing here is
    a shading decision.  ``radius`` keeps the renderer's hit code (-1 escaped,
    0 captured through the horizon, > 0 the Boyer-Lindquist radius of a disk
    hit); ``redshift`` is the observed/emitted frequency ratio g at a disk hit;
    ``phi`` the disk azimuth there; ``travel_time`` the coordinate time the
    light took to reach the camera; ``sky_theta``/``sky_phi`` the direction an
    escaping ray was travelling on the celestial sphere.
    """

    radius: np.ndarray
    redshift: np.ndarray
    phi: np.ndarray
    travel_time: np.ndarray
    sky_theta: np.ndarray
    sky_phi: np.ndarray


class DiskSamples(NamedTuple):
    """The disk hits of one frame, flattened to just the pixels that hit it.

    Built once per traced frame by `prepare_disk`, because none of it depends
    on animation time: the geodesics fixed which pixel sees which patch of
    disk, and only the patch's brightness changes from tick to tick.
    ``d_phi``/``d_log_radius`` are how far the disk coordinates move between
    neighbouring pixels, which is what the emissivity texture has to be
    band-limited against.
    """

    radius: np.ndarray
    redshift: np.ndarray
    phi: np.ndarray
    travel_time: np.ndarray
    log_radius: np.ndarray
    omega: np.ndarray
    d_phi: np.ndarray
    d_log_radius: np.ndarray
    d_omega: np.ndarray


class ShadingCache(NamedTuple):
    """Everything about a frame that survives from one animation tick to the next.

    The observer does not move while the flow animates, so the lensed sky, the
    silhouette and the disk's traced geometry are all fixed; only the disk's
    own emission is recomputed per tick.
    """

    base: np.ndarray  # the scene minus the disk: lensed background sky
    disk_mask: np.ndarray
    samples: DiskSamples


class DiskPhysics(NamedTuple):
    """Scene parameters the shading model needs that the buffer doesn't carry.

    These are not free artistic knobs: they must match what the renderer was
    actually called with, or the emission model stops describing the disk that
    was traced.  The defaults mirror ffi_bridge.cpp's hardcoded geometry and
    cli_parser.h's Params.
    """

    rs: float = 0.05           # Schwarzschild radius, 2M
    a: float = 0.0             # spin
    Q: float = 0.0             # charge
    inner_radius: float = 0.1  # gyuru_sugar_kicsi: the disk's real inner edge
    outer_radius: float = 0.5  # gyuru_sugar_nagy


# Wall-clock seconds for one full orbit of material at the disk's inner edge.
# This is the only free parameter of the flow animation - the playback speed.
# Everything else follows from the metric: how much slower each larger radius
# turns, and how far the light travel time sets each part of the image back in
# phase, are both fixed by the geometry, not chosen here.
INNER_ORBIT_SECONDS = 14.0

EXPOSURE = 8.0  # scales linear scene radiance into the tonemapper's range


def smoothstep(edge0: float, edge1: float, value: np.ndarray) -> np.ndarray:
    """A clamped smooth interpolation, safe when processing HDR values."""
    t = np.clip((value - edge0) / (edge1 - edge0), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def split_channels(width: int, height: int, values: np.ndarray) -> RayFrame:
    """Split a flat interleaved buffer - the layout szinsaver.h's datasaver_T
    writes and ffi_bridge.cpp's render_frame_f32 fills directly - into
    image-shaped channels."""
    expected = width * height
    if values.size % expected != 0:
        raise ValueError(f"Buffer of {values.size} floats is not a whole number of {expected}-pixel channels.")
    channels = values.size // expected
    if channels not in (1, 3, RAY_CHANNELS):
        raise ValueError(f"Expected 1, 3, or {RAY_CHANNELS} channels per pixel, received {channels}.")

    # C++ indexes [x * height + y], whereas images are row-major [y, x].
    buffer = values.reshape(width, height, channels)
    plane = [buffer[:, :, c].T.copy() for c in range(channels)]
    ones = np.ones((height, width), dtype=np.float32)
    zeros = np.zeros((height, width), dtype=np.float32)

    # Older buffers carry fewer channels.  A radius-only buffer has no
    # frequency information, so it gets an unshifted thermal spectrum; a
    # 3-channel one predates the light-travel-time and escape-direction
    # channels, so it renders without retarded-time lag or a lensed sky.
    if channels == 1:
        return RayFrame(plane[0], ones, zeros, zeros, zeros, zeros)
    if channels == 3:
        return RayFrame(plane[0], plane[1], plane[2], zeros, zeros, zeros)
    return RayFrame(*plane)


def read_hit_buffer(path: Path) -> RayFrame:
    with path.open("rb") as source:
        height = int.from_bytes(source.read(4), byteorder="little", signed=True)
        width = int.from_bytes(source.read(4), byteorder="little", signed=True)
        if width <= 0 or height <= 0:
            raise ValueError(f"Invalid buffer dimensions: {width} x {height}")
        values = np.frombuffer(source.read(), dtype="<f4")
    return split_channels(width, height, values)


def blackbody_lut() -> tuple[np.ndarray, np.ndarray]:
    """Return display-linear sRGB chromaticities from Planck spectra.

    The spectra are numerically integrated against analytic CIE 1931 matching
    curves.  A lookup table avoids allocating a wavelength axis per pixel.
    """
    wavelengths = np.arange(380.0, 781.0, 5.0)
    x = (
        0.362 * np.exp(-0.5 * ((wavelengths - 442.0) / np.where(wavelengths < 442.0, 16.0, 26.7)) ** 2)
        + 1.056 * np.exp(-0.5 * ((wavelengths - 599.8) / np.where(wavelengths < 599.8, 37.9, 31.0)) ** 2)
        - 0.065 * np.exp(-0.5 * ((wavelengths - 501.1) / np.where(wavelengths < 501.1, 20.4, 26.2)) ** 2)
    )
    y = (
        0.821 * np.exp(-0.5 * ((wavelengths - 568.8) / np.where(wavelengths < 568.8, 46.9, 40.5)) ** 2)
        + 0.286 * np.exp(-0.5 * ((wavelengths - 530.9) / np.where(wavelengths < 530.9, 16.3, 31.1)) ** 2)
    )
    z = (
        1.217 * np.exp(-0.5 * ((wavelengths - 437.0) / np.where(wavelengths < 437.0, 11.8, 36.0)) ** 2)
        + 0.681 * np.exp(-0.5 * ((wavelengths - 459.0) / np.where(wavelengths < 459.0, 26.0, 13.8)) ** 2)
    )
    temperatures = np.geomspace(900.0, 100_000.0, 512)
    wavelength_m = wavelengths * 1e-9
    exponent = 1.4387769e-2 / (temperatures[:, None] * wavelength_m[None, :])
    spectral_radiance = 1.0 / (wavelength_m[None, :] ** 5 * np.expm1(exponent))
    xyz = spectral_radiance @ np.stack((x, y, z), axis=1)
    rgb = xyz @ np.array(
        ((3.2406, -1.5372, -0.4986), (-0.9689, 1.8758, 0.0415), (0.0557, -0.2040, 1.0570)),
        dtype=np.float64,
    ).T
    rgb = np.maximum(rgb, 0.0)
    rgb /= np.maximum(rgb.max(axis=1, keepdims=True), 1e-12)
    return temperatures.astype(np.float32), rgb.astype(np.float32)


BLACKBODY_TEMPERATURES, BLACKBODY_RGB = blackbody_lut()


def blackbody_rgb(temperature: np.ndarray) -> np.ndarray:
    """Interpolate the Planck/CIE table, returning (..., 3) linear sRGB."""
    clamped = np.clip(temperature, BLACKBODY_TEMPERATURES[0], BLACKBODY_TEMPERATURES[-1])
    lookup = np.interp(clamped.ravel(), BLACKBODY_TEMPERATURES, np.arange(BLACKBODY_TEMPERATURES.size))
    lower = np.floor(lookup).astype(np.int32)
    upper = np.minimum(lower + 1, BLACKBODY_TEMPERATURES.size - 1)
    mix = (lookup - lower)[:, None]
    blended = BLACKBODY_RGB[lower] * (1.0 - mix) + BLACKBODY_RGB[upper] * mix
    return blended.reshape((*temperature.shape, 3))


# --------------------------------------------------------------------------
# The sky, as a texture on the celestial sphere
# --------------------------------------------------------------------------
#
# The renderer reports the direction each escaping ray was travelling when it
# left the integration volume, so the background can be sampled along that
# direction instead of being pasted onto the frame in screen space.  That is
# what makes the star field actually lensed: the same geodesics that bend the
# disk into a ring drag the background around the shadow, and the strongly
# deflected rays near the photon sphere sweep across the whole sky.

SKY_HEIGHT = 1024  # equirectangular texture: rows span theta, columns span phi
SKY_WIDTH = 2 * SKY_HEIGHT
SKY_LEVELS = 6  # mip levels; the top level is a nearly uniform sky glow


def _sky_base(seed: int) -> np.ndarray:
    """Build the level-0 equirectangular sky plate: dust band plus stars.

    Everything is a function of direction on the sphere, so it stays put in
    space while the camera moves - and, more importantly, it is what the
    lensing map samples, rather than something laid over the result.
    """
    rng = np.random.default_rng(seed)
    theta = (np.arange(SKY_HEIGHT, dtype=np.float32) + 0.5) * (np.pi / SKY_HEIGHT)
    phi = (np.arange(SKY_WIDTH, dtype=np.float32) + 0.5) * (2.0 * np.pi / SKY_WIDTH) - np.pi
    sin_theta = np.sin(theta)[:, None]
    cos_theta = np.cos(theta)[:, None]
    sin_phi = np.sin(phi)[None, :]
    cos_phi = np.cos(phi)[None, :]

    # A galactic plane tilted well off the hole's spin axis, so the dust lane
    # cuts across the frame instead of lying along the disk.
    pole = np.array([0.36, -0.52, 0.77], dtype=np.float32)
    pole /= np.linalg.norm(pole)
    height = pole[0] * sin_theta * cos_phi + pole[1] * sin_theta * sin_phi + pole[2] * cos_theta

    # Position along the band, used to give the dust filamentary structure.
    along = np.arctan2(cos_theta - pole[2] * height, sin_theta * cos_phi - pole[0] * height)
    band = np.exp(-((height / 0.26) ** 2))
    filaments = 0.45 + 0.55 * (
        0.5 + 0.5 * np.sin(3.0 * along + 7.0 * height + 0.9)
    ) * (0.5 + 0.5 * np.sin(-5.0 * along + 11.0 * height + 2.3))
    dust = band * filaments

    sky = np.empty((SKY_HEIGHT, SKY_WIDTH, 3), dtype=np.float32)
    sky[..., 0] = 0.0016 + dust * 0.0125
    sky[..., 1] = 0.0011 + dust * 0.0050
    sky[..., 2] = 0.0038 + dust * 0.0105

    # Stars: uniform on the sphere (uniform in cos theta, not in theta, or they
    # would pile up at the poles), coloured by an actual blackbody temperature
    # rather than a random tint, and concentrated towards the dust band.
    count = 60_000
    star_cos_theta = rng.uniform(-1.0, 1.0, size=count)
    star_phi = rng.uniform(-np.pi, np.pi, size=count)
    star_theta = np.arccos(star_cos_theta)
    row = np.clip((star_theta / np.pi * SKY_HEIGHT).astype(np.int32), 0, SKY_HEIGHT - 1)
    col = np.clip(((star_phi + np.pi) / (2.0 * np.pi) * SKY_WIDTH).astype(np.int32), 0, SKY_WIDTH - 1)

    # A steep luminosity function: overwhelmingly faint stars, a handful bright.
    brightness = 0.02 * (1.0 - rng.uniform(size=count)) ** -0.55
    brightness = np.minimum(brightness, 6.0)
    brightness *= 0.45 + 1.1 * band[row, col]  # the band is crowded with stars too
    temperature = np.clip(2600.0 * (1.0 - rng.uniform(size=count)) ** -0.32, 2600.0, 26_000.0)
    colour = blackbody_rgb(temperature.astype(np.float32))

    for channel in range(3):
        np.add.at(sky[..., channel], (row, col), brightness * colour[:, channel])
    return sky


def _box_downsample(image: np.ndarray) -> np.ndarray:
    """Average 2x2 texel blocks - the mip reduction that keeps mean radiance."""
    height, width = image.shape[:2]
    height -= height % 2
    width -= width % 2
    block = image[:height, :width].reshape(height // 2, 2, width // 2, 2, image.shape[2])
    return block.mean(axis=(1, 3), dtype=np.float32)


@functools.lru_cache(maxsize=4)
def sky_pyramid(seed: int) -> tuple[np.ndarray, ...]:
    """Mip pyramid of the sky plate, cached per seed (it costs ~a second).

    Lensing conserves surface brightness, so where the map squeezes a wide
    patch of sky into one pixel the right answer is the *average* radiance of
    that patch - exactly what a box-reduced mip level holds.  Sampling a lower
    level there is therefore not just an anti-aliasing trick: it is what turns
    the tangle of sky near the photon ring into the smooth bright band it
    physically is, instead of a field of aliased speckle that crawls whenever
    the camera moves.
    """
    levels = [_sky_base(seed)]
    for _ in range(SKY_LEVELS - 1):
        levels.append(_box_downsample(levels[-1]))
    return tuple(levels)


def _sample_level(level: np.ndarray, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Bilinear equirectangular lookup: wrapped in azimuth, clamped at the poles."""
    height, width = level.shape[:2]
    u = (phi / (2.0 * np.pi) + 0.5) * width - 0.5
    v = (theta / np.pi) * height - 0.5
    u0 = np.floor(u)
    v0 = np.floor(v)
    fu = (u - u0)[:, None]
    fv = (v - v0)[:, None]
    u0 = u0.astype(np.int64)
    v0 = v0.astype(np.int64)
    ua, ub = u0 % width, (u0 + 1) % width
    va = np.clip(v0, 0, height - 1)
    vb = np.clip(v0 + 1, 0, height - 1)
    top = level[va, ua] * (1.0 - fu) + level[va, ub] * fu
    bottom = level[vb, ua] * (1.0 - fu) + level[vb, ub] * fu
    return top * (1.0 - fv) + bottom * fv


def _neighbour_steps(values: np.ndarray, mask: np.ndarray, period: float = 0.0) -> np.ndarray:
    """Per-pixel distance to the next pixel's value along each image axis, maxed.

    `values` may be (h, w) or (h, w, c); the last axis of a 3-D input is treated
    as a vector and reduced by its norm.  Pixels outside `mask` come back NaN.

    A non-zero `period` folds each difference into a half-period either side of
    zero.  The disk azimuth needs it: two neighbouring rays straddling the
    camera's meridian went round the hole opposite ways, so the azimuths they
    report differ by nearly a full turn even though they landed on adjacent
    patches of the same disk.  Read literally, that looks like the disk racing
    past in one pixel, and the band limit would erase all texture along a line
    straight down the middle of the frame.
    """
    filled = np.where(mask[..., None] if values.ndim == 3 else mask, values, np.nan)
    span = None
    for axis in (0, 1):
        delta = np.diff(filled, axis=axis)
        if period > 0.0:
            delta = (delta + period / 2.0) % period - period / 2.0
        step = np.linalg.norm(delta, axis=-1) if values.ndim == 3 else np.abs(delta)
        pad = [(0, 0), (0, 0)]
        pad[axis] = (0, 1)
        step = np.pad(step, pad, mode="edge")
        span = step if span is None else np.fmax(span, step)
    return span


def _smooth_footprint(span: np.ndarray, mask: np.ndarray, fallback_scale: float) -> np.ndarray:
    """Average a raw footprint over its 3x3 neighbourhood, ignoring masked-out pixels.

    A one-sided difference between two pixels is a noisy estimator: the adaptive
    integrator's own tolerance jitters the traced coordinates by a fraction of a
    pixel, and taking the worst of two differences turns that jitter into an
    overestimate. Unsmoothed, that overestimate is largest exactly along the
    strongly deflected rays near the shadow, where it read as a blurred seam
    down the middle of the sky. Genuine compression by the lensing is broad
    enough to survive the averaging; single-pixel noise is not.
    """
    valid = mask & np.isfinite(span)
    values = np.where(valid, span, 0.0)
    weights = valid.astype(np.float32)
    total = np.zeros_like(values)
    count = np.zeros_like(weights)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            total += np.roll(np.roll(values, dy, axis=0), dx, axis=1)
            count += np.roll(np.roll(weights, dy, axis=0), dx, axis=1)
    averaged = np.divide(total, count, out=np.zeros_like(total), where=count > 0)
    if not valid.any():
        return np.full_like(span, fallback_scale)
    fallback = float(np.median(averaged[valid]))
    return np.where(valid & (count > 0), averaged, fallback)


def _pixel_solid_angle(theta: np.ndarray, phi: np.ndarray, escaped: np.ndarray) -> np.ndarray:
    """Angle subtended on the sky by one pixel, from the traced escape directions.

    The lensing map's local stretch is exactly the finite difference of
    neighbouring rays' sky directions, so this is measured off the geodesics
    rather than assumed.  Pixels next to the shadow (whose neighbours were
    captured, and so have no direction) fall back to the frame's median.
    """
    direction = np.stack(
        (np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)), axis=-1
    )
    # Chord length stands in for the angle: the two agree to well within a
    # pixel at these separations, and the vector form has no wraparound to
    # trip over where the azimuth crosses +-pi.
    span = _neighbour_steps(direction, escaped)
    return _smooth_footprint(span, escaped, np.pi / max(theta.shape[0], 1))


def lensed_sky(frame: RayFrame, seed: int) -> np.ndarray:
    """Radiance arriving from the background along each escaping ray."""
    escaped = frame.radius < 0.0
    scene = np.zeros((*frame.radius.shape, 3), dtype=np.float32)
    if not escaped.any():
        return scene

    levels = sky_pyramid(seed)
    span = _pixel_solid_angle(frame.sky_theta, frame.sky_phi, escaped)
    texel = 2.0 * np.pi / SKY_WIDTH
    lod = np.clip(np.log2(np.maximum(span / texel, 1e-6)), 0.0, SKY_LEVELS - 1.0)

    theta = frame.sky_theta[escaped]
    phi = frame.sky_phi[escaped]
    lod = lod[escaped]
    sampled = np.zeros((theta.size, 3), dtype=np.float32)
    for index, level in enumerate(levels):
        weight = np.maximum(0.0, 1.0 - np.abs(lod - index))
        hit = weight > 0.0
        if not hit.any():
            continue
        sampled[hit] += _sample_level(level, theta[hit], phi[hit]) * weight[hit][:, None]

    scene[escaped] = sampled
    return scene


# --------------------------------------------------------------------------
# The accretion disk
# --------------------------------------------------------------------------


def orbital_angular_velocity(radius: np.ndarray, physics: DiskPhysics) -> np.ndarray:
    """Coordinate angular velocity d(phi)/dt of the disk material at `radius`.

    This is the *same* prograde circular-orbit expression cuda_ray.h's
    disk_redshift uses to build the emitter's four-velocity.  Reusing it, and
    not a Newtonian stand-in, is what keeps the animation honest: the speed the
    texture is seen to flow at is now the speed the Doppler shift and beaming
    baked into every pixel were computed for.  Because it is a coordinate
    angular velocity it is already the rate a distant static camera measures,
    with the gravitational time dilation of the inner disk included.
    """
    root = np.sqrt(np.maximum(physics.rs * radius / 2.0 - physics.Q**2, 0.0))
    return root / np.maximum(radius**2 + physics.a * root, 1e-12)


def flow_time_scale(physics: DiskPhysics) -> float:
    """Coordinate time per wall-clock second, set so the inner edge orbits in
    INNER_ORBIT_SECONDS.  Only the playback rate is chosen here; every radius's
    speed *relative* to the inner edge comes from the metric."""
    inner_omega = float(orbital_angular_velocity(np.array([physics.inner_radius]), physics)[0])
    if inner_omega <= 0.0:
        return 0.0
    return (2.0 * np.pi / inner_omega) / INNER_ORBIT_SECONDS


def eddy_lifetime(physics: DiskPhysics) -> float:
    """Coordinate time an eddy pattern survives before it is replaced.

    One orbit of the disk's inner edge.  Turbulence in a real disk is not a
    frozen pattern being wound up forever - it is continuously regenerated on
    roughly an orbital timescale, so its structure stays at a fixed physical
    scale instead of being sheared ever finer.  Letting a single realization
    run indefinitely does the opposite: differential rotation drags
    neighbouring radii apart without limit, and everything fine eventually
    winds below what a pixel can resolve, leaving a bland smooth disc after a
    minute or two of watching.  Capping the advection at one lifetime bounds
    how far shear can separate neighbours, so the disk keeps its texture.
    """
    inner_omega = float(orbital_angular_velocity(np.array([physics.inner_radius]), physics)[0])
    return 2.0 * np.pi / inner_omega if inner_omega > 0.0 else 1.0


def _regeneration(coordinate_time: float, physics: DiskPhysics) -> tuple[int, float, float]:
    """Split coordinate time into (generation, time since its epoch, cross-fade).

    The epoch is global, not per-pixel: every radius restarts its eddies at the
    same instant, so the advection two neighbouring pixels have accumulated can
    never drift apart by more than one lifetime's worth of shear.
    """
    lifetime = eddy_lifetime(physics)
    generation = math.floor(coordinate_time / lifetime)
    elapsed = coordinate_time - generation * lifetime
    return generation, elapsed, float(smoothstep(0.0, 1.0, np.asarray(elapsed / lifetime)))


def _sampling_footprint(values: np.ndarray, mask: np.ndarray, period: float = 0.0) -> np.ndarray:
    """How far `values` moves between neighbouring pixels, within `mask`.

    Measured off the traced buffers themselves, so it follows the lensing:
    where the map crams a wide sweep of disk into few pixels - the far side
    seen through the hole, the higher-order images stacked at the photon ring -
    this comes back large, and the texture knows to stop resolving detail it
    cannot sample.
    """
    return _smooth_footprint(_neighbour_steps(values, mask, period), mask, 0.0)


def prepare_disk(frame: RayFrame, physics: DiskPhysics) -> DiskSamples:
    """Flatten the frame's disk hits and derive the time-independent quantities."""
    disk = frame.radius > 0.0
    radius = frame.radius[disk]
    scaled = np.maximum(radius / max(physics.inner_radius, 1e-6), 1.0)
    omega = orbital_angular_velocity(radius, physics)

    image = np.zeros_like(frame.radius)
    image[disk] = np.log(scaled)
    d_log_radius = _sampling_footprint(image, disk)[disk]
    image[disk] = omega
    d_omega = _sampling_footprint(image, disk)[disk]
    image[disk] = omega * frame.travel_time[disk]

    # The azimuth of the material a pixel sees is the traced azimuth minus how
    # far it has turned since the eddies were laid down, so how fast the disk
    # coordinates move across the image has to account for the retarded time
    # too: light travel time varies sharply from pixel to pixel exactly where
    # the lensing is strongest, and that shears the pattern there just as
    # surely as the rotation does.
    d_advect_static = _sampling_footprint(image, disk)[disk]

    return DiskSamples(
        radius=radius,
        redshift=frame.redshift[disk],
        phi=frame.phi[disk],
        travel_time=frame.travel_time[disk],
        log_radius=np.log(scaled),
        omega=omega,
        d_phi=_sampling_footprint(frame.phi, disk, 2.0 * np.pi)[disk] + d_advect_static,
        d_log_radius=d_log_radius,
        d_omega=d_omega,
    )


def disk_turbulence(
    log_radius: np.ndarray,
    phi: np.ndarray,
    orbit_phase: np.ndarray,
    d_phi: np.ndarray,
    d_log_radius: np.ndarray,
    seed: int,
    generation: int = 0,
) -> np.ndarray:
    """A fractal, sheared, evolving turbulence field standing in for MRI eddies.

    Real thin-disk turbulence is wound into trailing filaments by Keplerian
    differential rotation (faster shear at small radius), then cascades into
    progressively finer structure.  This sums several octaves at incommensurate
    (phi, log r) frequencies -- each octave sheared harder than the last, as
    real turbulence would be -- and warps their common domain with a slower
    flow field so the result reads as swirling density structure rather than a
    stack of clean concentric ripples.  It stays a smooth function of the
    ray-traced disk coordinates (r, phi), so it warps correctly under lensing
    instead of being a screen-space overlay.

    `phi` is expected to already be the *material* azimuth: the label of the
    fluid element, with its orbital advection removed by the caller.  What this
    adds on top is the eddies' own life cycle.  A rigidly rotating texture
    reads as a painted disc being spun; real eddies swell, shear out and fade
    on roughly the local orbital timescale, so `orbit_phase` (elapsed local
    orbits, which is short at small radius and long further out, exactly as the
    turnover time is) modulates each octave's amplitude.  Fine octaves churn
    through several cycles per orbit, coarse ones barely change.

    The renderer reports phi as the integrated Boyer-Lindquist azimuth, so a
    strongly deflected ray that wound most of the way around the hole comes
    back with phi well outside (-pi, pi] -- values past +-3*pi show up in a
    normal frame.  Physically that is the same material as the same azimuth one
    turn earlier, so the texture must be exactly 2*pi-periodic in phi, and
    sin(k * phi) only is when k is an integer.  Every coefficient multiplying
    phi here (including inside the warp field) is therefore integral.
    log_radius has no such wraparound and can use any real coefficient; so can
    orbit_phase, which is a phi-independent modulation at fixed radius.

    Each octave is faded out where one pixel spans more than about half a
    period of it (`d_phi`/`d_log_radius` say how far the disk coordinates move
    between neighbouring pixels).  Below that limit a pixel cannot represent
    the octave, only alias it - which is why the tightly wound inner disk and
    the squeezed higher-order images used to break into moire rings, and why
    that moire would crawl and shimmer as soon as the flow started moving.
    Faded octaves settle to their own mean rather than to zero, so detail
    dissolves into smooth surface brightness the way it does when a real
    instrument stops resolving it, instead of leaving dark bands.

    `generation` selects which realization of the eddy field this is; see
    `_regeneration`, which cross-fades successive generations so that no
    pattern is ever advected for longer than eddies actually survive.
    """
    # Masked, not passed straight through: a shutter that opens just before
    # t = 0 puts the generation counter negative, which seeding rejects.
    rng = np.random.default_rng((seed + 97, generation & 0xFFFFFFFF))
    warp_phase = rng.uniform(0.0, 2.0 * np.pi, size=2)
    warp = 0.6 * np.sin(2.0 * phi - 4.3 * log_radius + warp_phase[0]) + 0.4 * np.sin(
        -3.0 * phi + 2.6 * log_radius + warp_phase[1]
    )
    phi_w = phi + 0.22 * warp
    logr_w = log_radius + 0.15 * warp

    # (phi frequency, log-radius frequency, eddy turnovers per local orbit) per
    # octave.  The frequency ratio grows each octave, mimicking how
    # differential rotation shears small eddies into tighter trailing spirals
    # than large ones; the turnover rate grows with it, because smaller eddies
    # live and die faster.
    octaves = (
        (5.0, 11.0, 0.25),
        (9.0, -23.0, 0.5),
        (17.0, 44.0, 1.0),
        (33.0, -81.0, 1.9),
        (61.0, 149.0, 3.4),
        (113.0, -277.0, 6.1),
    )
    phases = rng.uniform(0.0, 2.0 * np.pi, size=len(octaves))
    churn = rng.uniform(0.0, 2.0 * np.pi, size=len(octaves))

    def resolved(k_phi: float, k_r: float) -> np.ndarray:
        """How much of an octave at these frequencies a pixel can still carry.

        The domain warp above stretches the local frequencies somewhat, so the
        step is scaled up a little before being compared against Nyquist.
        """
        step = 1.35 * (abs(k_phi) * d_phi + abs(k_r) * d_log_radius)
        return 1.0 - smoothstep(0.5 * np.pi, np.pi, step)

    fractal = np.zeros_like(phi_w)
    amplitude, total = 1.0, 0.0
    for (k_phi, k_r, turnover), phase, churn_phase in zip(octaves, phases, churn):
        life = 0.62 + 0.38 * np.sin(2.0 * np.pi * turnover * orbit_phase + 0.7 * logr_w + churn_phase)
        detail = resolved(k_phi, k_r) * np.sin(k_phi * phi_w + k_r * logr_w + phase)
        fractal += amplitude * life * (0.5 + 0.5 * detail)
        total += amplitude
        amplitude *= 0.68
    fractal /= total

    # A slower beat between two mid-frequency octaves picks out patchy hot
    # spots/eddies rather than uniformly striping the whole disk.  They drift
    # in and out over about an orbit, the timescale such a clump would survive.
    clump_life = 0.5 + 0.5 * np.sin(2.0 * np.pi * 0.7 * orbit_phase + phases[1])
    clump_detail = resolved(8.0, 13.0)
    clumps = (0.5 + 0.5 * clump_detail * np.sin(8.0 * phi_w + 3.0 * logr_w + phases[0])) * (
        0.5 + 0.5 * clump_detail * np.sin(-13.0 * logr_w + 6.0 * phi_w + phases[2])
    )
    return 0.34 + 1.32 * fractal + 0.8 * (clumps * clump_life) ** 3


def disk_radiance(
    samples: DiskSamples,
    peak_temperature: float,
    seed: int,
    time: float = 0.0,
    physics: DiskPhysics = DiskPhysics(),
) -> np.ndarray:
    """Thermal thin-disk emission transported by the ray-traced redshift.

    Liouville's theorem gives I_nu / nu^3 as an invariant.  For a blackbody,
    this means the observed colour temperature is g*T and the bolometric
    intensity gains g^4, where g is calculated by the C++ geodesic renderer.

    `time` is wall-clock seconds of animation.  Each pixel is evaluated at its
    own *retarded* time - the light reaching it left the disk `travel_time`
    earlier, and the renderer measured that per ray - so the strongly lensed
    images, whose light took a longer way round, correctly show the disk as it
    was further in the past than the direct image does.  With a camera this
    close, that lag is a large fraction of an inner orbit, and it is why the
    lensed arcs run out of step with the disk that casts them instead of
    sliding along in lockstep the way a single-clock texture would.

    Returns one linear-light RGB triple per disk pixel, in `samples` order.
    """
    if samples.radius.size == 0:
        # No ray in this frame hit the disk (e.g. zoomed to a viewpoint where
        # the disk annulus falls entirely outside the camera's field of view).
        return np.zeros((0, 3), dtype=np.float32)

    # The inner edge is a fixed property of the scene the renderer traced, not
    # something to re-estimate per frame from whichever radii happen to be on
    # screen: doing that made the whole temperature and flux profile drift
    # every time the camera moved or zoomed, so the disk changed colour and
    # brightness for no physical reason.
    scaled_radius = np.exp(samples.log_radius)
    no_torque = np.maximum(1.0 - np.sqrt(1.0 / scaled_radius), 0.015)
    emitted_temperature = peak_temperature * scaled_radius**-0.75 * no_torque**0.25
    observed_temperature = emitted_temperature * samples.redshift

    # Retarded time, in coordinate units: the light a pixel shows left the disk
    # `travel_time` before now.  Two generations of eddies are evaluated at it
    # and cross-faded, so the pattern being advected is never older than eddies
    # live (see `_regeneration`).
    lifetime = eddy_lifetime(physics)
    generation, elapsed, blend = _regeneration(time * flow_time_scale(physics), physics)
    texture = None
    for step, weight in ((0, 1.0 - blend), (1, blend)):
        if weight <= 0.0:
            continue
        # Subtracting the turned-through angle from the traced azimuth labels
        # the fluid element, so the texture is carried around by the flow.
        emission_time = elapsed - step * lifetime - samples.travel_time
        turned = samples.omega * emission_time
        # How fast that angle changes from pixel to pixel sets what the frame
        # can still resolve, and it grows with how long the pattern has been
        # advected - which is exactly why that advection is capped.
        d_material_phi = samples.d_phi + abs(elapsed - step * lifetime) * samples.d_omega
        octave_field = disk_turbulence(
            samples.log_radius,
            samples.phi - turned,
            turned / (2.0 * np.pi),
            d_material_phi,
            samples.d_log_radius,
            seed,
            generation + step,
        )
        texture = octave_field * weight if texture is None else texture + octave_field * weight

    # The traced annulus stops dead at gyuru_sugar_nagy, which reads as a
    # knife-edged disc. Fading the emissivity out over the last of it is a
    # statement about where the gas thins out, not about where rays were
    # allowed to hit: the geometry the renderer integrated is untouched. It
    # stays a narrow band on purpose - the renderer still stops every ray that
    # reaches the annulus, so anywhere the emissivity is faded but the geometry
    # is not reads as an opaque dark rim rather than as thinning gas.
    outer_edge = smoothstep(physics.outer_radius, physics.outer_radius * 0.96, samples.radius)

    colour = blackbody_rgb(observed_temperature)
    flux = scaled_radius**-3.0 * no_torque
    intensity = flux * texture * outer_edge * np.clip(samples.redshift, 0.05, 5.0) ** 4
    return colour * (intensity * EXPOSURE)[:, None]


def base_scene(frame: RayFrame, seed: int, physics: DiskPhysics = DiskPhysics()) -> ShadingCache:
    """Everything about the frame that the flow animation must not recompute.

    The observer never moves during flow animation - only the disk's own
    emission should change tick to tick - so the lensed sky and the disk's
    traced geometry are resolved once here and reused, rather than re-sampling
    the whole background on every tick.
    """
    scene = lensed_sky(frame, seed)
    # Captured rays are left black: they carry no light, and the sky already
    # supplies the glow that rings the shadow, from real strongly-deflected
    # geodesics rather than a blur painted over the silhouette.
    return ShadingCache(base=scene, disk_mask=frame.radius > 0.0, samples=prepare_disk(frame, physics))


# --------------------------------------------------------------------------
# Film finish
# --------------------------------------------------------------------------


def _resize(image: np.ndarray, height: int, width: int, resample: int) -> np.ndarray:
    channels = [
        np.asarray(
            Image.fromarray(image[..., c]).resize((width, height), resample=resample),
            dtype=np.float32,
        )
        for c in range(image.shape[2])
    ]
    return np.stack(channels, axis=-1)


def _bloom(scene: np.ndarray, octaves: int = 5) -> np.ndarray:
    """Multi-scale bloom: a sum of successively wider, weaker halos.

    A lens scatters a little light at every angular scale, so a single Gaussian
    reads as a flat sticker around bright areas while a sum over octaves gives
    the tight core and long faint skirt real glare has.  Computed on the linear
    HDR scene, before tonemapping, so the brightest parts bloom in proportion
    to how bright they actually are instead of after being clipped to white.
    """
    height, width = scene.shape[:2]
    current = scene
    total = np.zeros_like(scene)
    weight_sum = 0.0
    weight = 1.0
    for _ in range(octaves):
        if min(current.shape[0], current.shape[1]) <= 2:
            break
        current = _resize(current, max(current.shape[0] // 2, 1), max(current.shape[1] // 2, 1), Image.BOX)
        total += _resize(current, height, width, Image.BILINEAR) * weight
        weight_sum += weight
        weight *= 0.55
    return total / max(weight_sum, 1e-6)


@functools.lru_cache(maxsize=8)
def _vignette(height: int, width: int) -> np.ndarray:
    yy, xx = np.mgrid[0:height, 0:width]
    nx = (xx / max(width - 1, 1) - 0.5) * 2.0
    ny = (yy / max(height - 1, 1) - 0.5) * 2.0
    falloff = 1.0 - 0.30 * np.clip(nx * nx + ny * ny, 0.0, 1.0)
    return (falloff**1.4).astype(np.float32)[..., None]


def finish_frame(scene: np.ndarray, seed: int, frame_index: int = 0) -> Image.Image:
    """Bloom, highlight rolloff, ACES tonemap, vignette and grain: the shared
    finishing pass for a composited linear-light scene. Bloom has to re-run every
    frame (the disk's changing brightness spreads into neighbouring pixels), but
    it runs on a downsampled pyramid - the expensive part (sampling the lensed
    sky) already happened once in base_scene.

    `frame_index` advances the grain so it resamples each frame the way film
    does. Held fixed it stops reading as grain at all and becomes a dirty
    sensor: a static speckle pattern welded to the image.
    """
    height, width = scene.shape[:2]
    scene = scene + _bloom(scene) * 0.20

    # Bright, saturated highlights desaturate towards white as a real sensor's
    # channels saturate one after another; without it the hot inner disk clips
    # to a flat blue-white edge instead of blooming out to white.
    luminance = scene @ np.array([0.2126, 0.7152, 0.0722], dtype=np.float32)
    wash = smoothstep(0.85, 4.5, luminance)[..., None]
    scene = scene * (1.0 - wash) + luminance[..., None] * wash

    # ACES fitted tone map plus a light vignette focuses the frame.
    scene = (scene * (2.51 * scene + 0.03)) / (scene * (2.43 * scene + 0.59) + 0.14)
    scene = scene * _vignette(height, width)
    scene = np.clip(scene, 0.0, 1.0) ** (1.0 / 2.2)
    grain = np.random.default_rng((seed + 1, frame_index)).normal(0.0, 0.006, size=(height, width, 1))
    scene = np.clip(scene + grain, 0.0, 1.0)
    return Image.fromarray((scene * 255.0 + 0.5).astype(np.uint8))


def composite(
    cache: ShadingCache,
    peak_temperature: float,
    seed: int,
    time: float = 0.0,
    physics: DiskPhysics = DiskPhysics(),
    shutter: float = 0.0,
    shutter_samples: int = 1,
) -> np.ndarray:
    """Composite the disk onto the cached base scene at animation time `time`.

    `shutter` is how long the camera's shutter stays open, in the same
    wall-clock seconds as `time`; the disk is averaged over that many instants
    spread across it.  A real camera integrates over its exposure instead of
    sampling one instant, and at the inner edge - which sweeps several pixels
    per frame - that is the difference between material that flows and material
    that strobes from one position to the next.
    """
    scene = cache.base.copy()
    if cache.samples.radius.size == 0:
        return scene
    count = max(1, shutter_samples if shutter > 0.0 else 1)
    accumulated = None
    for index in range(count):
        offset = 0.0 if count == 1 else shutter * ((index + 0.5) / count - 0.5)
        radiance = disk_radiance(cache.samples, peak_temperature, seed, time + offset, physics)
        accumulated = radiance if accumulated is None else accumulated + radiance
    scene[cache.disk_mask] = accumulated / count
    return scene


def cinematic_image(
    frame: RayFrame,
    seed: int,
    peak_temperature: float,
    time: float = 0.0,
    physics: DiskPhysics = DiskPhysics(),
    frame_index: int = 0,
) -> Image.Image:
    cache = base_scene(frame, seed, physics)
    scene = composite(cache, peak_temperature, seed, time, physics)
    return finish_frame(scene, seed, frame_index)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a cinematic image from a ray-tracing hit buffer.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="binary float hit buffer")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="PNG output path")
    parser.add_argument("--seed", type=int, default=23071969, help="fixed star-field seed")
    parser.add_argument("--peak-temperature", type=float, default=12_000.0, help="peak local disk temperature in kelvin")
    parser.add_argument("--time", type=float, default=0.0, help="animation time in seconds")
    parser.add_argument("--spin", type=float, default=0.0, help="black hole spin a the buffer was traced with")
    parser.add_argument("--charge", type=float, default=0.0, help="black hole charge Q the buffer was traced with")
    args = parser.parse_args()

    frame = read_hit_buffer(args.input)
    physics = DiskPhysics(a=args.spin, Q=args.charge)
    image = cinematic_image(frame, args.seed, args.peak_temperature, args.time, physics)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output)
    print(f"Saved cinematic render: {args.output} ({image.width} x {image.height})")


if __name__ == "__main__":
    main()
