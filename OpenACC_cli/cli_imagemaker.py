#!/usr/bin/env python3
"""Turn the ray tracer's scalar hit buffer into a cinematic black-hole image.

The C++ renderer remains the source of truth for ray paths.  Its buffer uses
``-1`` for rays that escape, ``0`` for captured rays, and a positive radius for
an accretion-disk hit.  New buffers additionally contain the ray-traced
frequency shift at every disk hit.  This script applies a thermal disk model,
the gravitational/Doppler shift, and a small deterministic film finish.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageFilter


DEFAULT_INPUT = Path("./web_images/kep_cli.dat")
DEFAULT_OUTPUT = Path("./web_images/blackhole_cli.png")


def smoothstep(edge0: float, edge1: float, value: np.ndarray) -> np.ndarray:
    """A clamped smooth interpolation, safe when processing HDR values."""
    t = np.clip((value - edge0) / (edge1 - edge0), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def read_hit_buffer(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with path.open("rb") as source:
        height = int.from_bytes(source.read(4), byteorder="little", signed=True)
        width = int.from_bytes(source.read(4), byteorder="little", signed=True)
        if width <= 0 or height <= 0:
            raise ValueError(f"Invalid buffer dimensions: {width} x {height}")
        values = np.frombuffer(source.read(), dtype="<f4")

    expected = width * height
    if values.size not in (expected, expected * 2, expected * 3):
        raise ValueError(
            f"Expected {expected}, {expected * 2}, or {expected * 3} float values in {path}, received {values.size}."
        )
    # C++ indexes [x * height + y], whereas images are row-major [y, x].
    if values.size == expected:
        # Compatibility with legacy radius-only buffers.  They have no
        # frequency information, so an unshifted thermal spectrum is used.
        radius = values.reshape(width, height).T.copy()
        return radius, np.ones((height, width), dtype=np.float32), np.zeros((height, width), dtype=np.float32)
    channels = values.size // expected
    buffer = values.reshape(width, height, channels)
    phi = buffer[:, :, 2].T.copy() if channels == 3 else np.zeros((height, width), dtype=np.float32)
    return buffer[:, :, 0].T.copy(), buffer[:, :, 1].T.copy(), phi


def stars_and_nebula(height: int, width: int, seed: int) -> np.ndarray:
    """A restrained, reproducible deep-space plate behind escaped rays."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:height, 0:width]
    x = (xx / max(width - 1, 1) - 0.5) * 2.0
    y = (yy / max(height - 1, 1) - 0.5) * 2.0

    # A low, warm galactic dust lane with a cool violet falloff.  It stays
    # below disk exposure but prevents the background reading as empty black.
    dust = np.exp(-((y + 0.13 * np.sin(x * 2.7)) / 0.40) ** 2)
    dust *= 0.30 + 0.70 * np.sin(x * 5.4 + y * 8.3 + 0.7) ** 2
    haze = np.exp(-((x + 0.25) ** 2 + (y + 0.18) ** 2) / 0.78)
    scene = np.empty((height, width, 3), dtype=np.float32)
    scene[..., 0] = 0.004 + dust * 0.023 + haze * 0.006
    scene[..., 1] = 0.002 + dust * 0.007 + haze * 0.002
    scene[..., 2] = 0.008 + dust * 0.015 + haze * 0.012

    # Draw stars as additive point sprites, with a few coloured bright stars.
    count = max(70, width * height // 500)
    sx = rng.integers(1, max(width - 1, 2), size=count)
    sy = rng.integers(1, max(height - 1, 2), size=count)
    brightness = rng.power(6.5, size=count) * 0.28 + 0.012
    tint = rng.uniform(0.76, 1.22, size=(count, 3))
    tint[:, 0] *= rng.uniform(0.85, 1.22, size=count)
    tint[:, 2] *= rng.uniform(0.88, 1.32, size=count)
    for px, py, intensity, color in zip(sx, sy, brightness, tint):
        scene[py, px] += intensity * color
        if intensity > 0.22:
            scene[py - 1 : py + 2, px] += intensity * color * 0.13
            scene[py, px - 1 : px + 2] += intensity * color * 0.13
    return scene


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


def disk_turbulence(log_radius: np.ndarray, phi: np.ndarray, seed: int) -> np.ndarray:
    """A fractal, sheared turbulence field standing in for MRI eddies.

    Real thin-disk turbulence is wound into trailing filaments by Keplerian
    differential rotation (faster shear at small radius), then cascades into
    progressively finer structure.  This sums several octaves at
    incommensurate (phi, log r) frequencies -- each octave sheared harder
    than the last, as real turbulence would be -- and warps their common
    domain with a slower flow field so the result reads as swirling density
    structure rather than a stack of clean concentric ripples.  It stays a
    smooth function of the ray-traced disk coordinates (r, phi), so it warps
    correctly under lensing instead of being a screen-space overlay.
    """
    # The renderer reports phi wrapped to (-pi, pi], so two neighbouring rays
    # straddling the camera's meridian plane land on opposite sides of that
    # branch cut (e.g. +3.1399 next to -3.1399 -- the same disk azimuth).
    # sin(k * phi) is only exactly 2*pi-periodic, and therefore blind to that
    # cut, when k is an integer, so every coefficient multiplying phi here
    # (including inside the warp field) must stay integral. log_radius has no
    # such wraparound and can use any real coefficient.
    rng = np.random.default_rng(seed + 97)
    warp_phase = rng.uniform(0.0, 2.0 * np.pi, size=2)
    warp = 0.6 * np.sin(2.0 * phi - 4.3 * log_radius + warp_phase[0]) + 0.4 * np.sin(
        -3.0 * phi + 2.6 * log_radius + warp_phase[1]
    )
    phi_w = phi + 0.22 * warp
    logr_w = log_radius + 0.15 * warp

    # (phi frequency, log-radius frequency) per octave.  The ratio grows each
    # octave, mimicking how differential rotation shears small eddies into
    # tighter trailing spirals than large ones.
    octaves = ((5.0, 11.0), (9.0, -23.0), (17.0, 44.0), (33.0, -81.0), (61.0, 149.0))
    phases = rng.uniform(0.0, 2.0 * np.pi, size=len(octaves))
    fractal = np.zeros_like(phi_w)
    amplitude, total = 1.0, 0.0
    for (k_phi, k_r), phase in zip(octaves, phases):
        fractal += amplitude * (0.5 + 0.5 * np.sin(k_phi * phi_w + k_r * logr_w + phase))
        total += amplitude
        amplitude *= 0.56
    fractal /= total

    # A slower beat between two mid-frequency octaves picks out patchy hot
    # spots/eddies rather than uniformly striping the whole disk.
    clumps = (0.5 + 0.5 * np.sin(8.0 * phi_w + 3.0 * logr_w + phases[0])) * (
        0.5 + 0.5 * np.sin(-13.0 * logr_w + 6.0 * phi_w + phases[2])
    )
    return 0.55 + 0.9 * fractal + 0.55 * clumps**3


def disk_radiance(
    radius: np.ndarray, redshift: np.ndarray, phi: np.ndarray, peak_temperature: float, seed: int
) -> np.ndarray:
    """Thermal thin-disk emission transported by the ray-traced redshift.

    Liouville's theorem gives I_nu / nu^3 as an invariant.  For a blackbody,
    this means the observed colour temperature is g*T and the bolometric
    intensity gains g^4, where g is calculated by the C++ geodesic renderer.
    """
    positive = radius[radius > 0.0]
    if positive.size == 0:
        # No ray in this frame hit the disk (e.g. zoomed to a viewpoint where
        # the disk annulus falls entirely outside the camera's field of
        # view). The caller only ever reads this back through the `disk`
        # mask, which is empty in that case too, so any correctly-shaped
        # array is a safe no-op here.
        return np.zeros((*radius.shape, 3), dtype=np.float32)
    inner = float(np.percentile(positive, 1))
    scaled_radius = np.maximum(radius / max(inner, 1e-6), 1.0)
    no_torque = np.maximum(1.0 - np.sqrt(1.0 / scaled_radius), 0.015)
    emitted_temperature = peak_temperature * scaled_radius ** -0.75 * no_torque**0.25
    observed_temperature = np.clip(emitted_temperature * redshift, 900.0, 100_000.0)

    # The emissivity texture lives in (r, phi) on the disk, so it is warped
    # by the same geodesics as the thermal emission.  It represents restrained
    # spiral density/turbulence structure, not a screen-space overlay.
    log_radius = np.log(scaled_radius)
    texture = disk_turbulence(log_radius, phi, seed)

    # Interpolate the Planck/CIE table for each pixel and apply the disk's
    # radial flux profile plus the invariant g^4 brightness transport.
    lookup = np.interp(observed_temperature.ravel(), BLACKBODY_TEMPERATURES, np.arange(BLACKBODY_TEMPERATURES.size))
    lower = np.floor(lookup).astype(np.int32)
    upper = np.minimum(lower + 1, BLACKBODY_TEMPERATURES.size - 1)
    mix = (lookup - lower)[:, None]
    colour = (BLACKBODY_RGB[lower] * (1.0 - mix) + BLACKBODY_RGB[upper] * mix).reshape((*radius.shape, 3))
    flux = scaled_radius ** -3.0 * no_torque
    intensity = flux * texture * np.clip(redshift, 0.05, 5.0) ** 4
    return colour * (intensity * 8.0)[..., None]


def cinematic_image(hits: np.ndarray, redshift: np.ndarray, phi: np.ndarray, seed: int, peak_temperature: float) -> Image.Image:
    height, width = hits.shape
    escaped = hits < 0.0
    captured = hits == 0.0
    disk = hits > 0.0

    scene = np.zeros((height, width, 3), dtype=np.float32)
    scene[escaped] = stars_and_nebula(height, width, seed)[escaped]
    scene[disk] = disk_radiance(hits, redshift, phi, peak_temperature, seed)[disk]

    # Light leaking around the event horizon makes the silhouette read clearly
    # without painting over the physically captured (zero-valued) rays.
    disk_mask = Image.fromarray((disk * 255).astype(np.uint8), mode="L")
    halo = np.asarray(disk_mask.filter(ImageFilter.GaussianBlur(radius=max(1.2, width / 270)))) / 255.0
    ring = halo * captured
    scene[..., 0] += ring * 0.18
    scene[..., 1] += ring * 0.052
    scene[..., 2] += ring * 0.012

    # Bloom only the luminous disk and ring, then use an ACES-like film curve.
    emission = np.clip(scene * 1.8, 0.0, 1.0)
    bloom = np.asarray(
        Image.fromarray((emission * 255).astype(np.uint8), mode="RGB").filter(
            ImageFilter.GaussianBlur(radius=max(1.5, width / 170))
        ),
        dtype=np.float32,
    ) / 255.0
    scene += bloom * 0.34

    # ACES fitted tone map plus a very light vignette focuses the frame.
    scene = (scene * (2.51 * scene + 0.03)) / (scene * (2.43 * scene + 0.59) + 0.14)
    yy, xx = np.mgrid[0:height, 0:width]
    nx = (xx / max(width - 1, 1) - 0.5) * 2.0
    ny = (yy / max(height - 1, 1) - 0.5) * 2.0
    vignette = 1.0 - 0.36 * np.clip(nx * nx + ny * ny, 0.0, 1.0)
    scene *= vignette[..., None]
    scene = np.clip(scene, 0.0, 1.0) ** (1.0 / 2.2)
    grain = np.random.default_rng(seed + 1).normal(0.0, 0.006, size=(height, width, 1))
    scene = np.clip(scene + grain, 0.0, 1.0)
    return Image.fromarray((scene * 255.0 + 0.5).astype(np.uint8), mode="RGB")


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a cinematic image from a ray-tracing hit buffer.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="binary float hit buffer")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="PNG output path")
    parser.add_argument("--seed", type=int, default=23071969, help="fixed star-field seed")
    parser.add_argument("--peak-temperature", type=float, default=12_000.0, help="peak local disk temperature in kelvin")
    args = parser.parse_args()

    hits, redshift, phi = read_hit_buffer(args.input)
    image = cinematic_image(hits, redshift, phi, args.seed, args.peak_temperature)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output)
    print(f"Saved cinematic render: {args.output} ({image.width} x {image.height})")


if __name__ == "__main__":
    main()
