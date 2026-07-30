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


def read_hit_buffer(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open("rb") as source:
        height = int.from_bytes(source.read(4), byteorder="little", signed=True)
        width = int.from_bytes(source.read(4), byteorder="little", signed=True)
        if width <= 0 or height <= 0:
            raise ValueError(f"Invalid buffer dimensions: {width} x {height}")
        values = np.frombuffer(source.read(), dtype="<f4")

    expected = width * height
    if values.size not in (expected, expected * 2):
        raise ValueError(
            f"Expected {expected} or {expected * 2} float values in {path}, received {values.size}."
        )
    # C++ indexes [x * height + y], whereas images are row-major [y, x].
    if values.size == expected:
        # Compatibility with legacy radius-only buffers.  They have no
        # frequency information, so an unshifted thermal spectrum is used.
        return values.reshape(width, height).T.copy(), np.ones((height, width), dtype=np.float32)
    return (
        values.reshape(width, height, 2)[:, :, 0].T.copy(),
        values.reshape(width, height, 2)[:, :, 1].T.copy(),
    )


def stars_and_nebula(height: int, width: int, seed: int) -> np.ndarray:
    """A restrained, reproducible deep-space plate behind escaped rays."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:height, 0:width]
    x = (xx / max(width - 1, 1) - 0.5) * 2.0
    y = (yy / max(height - 1, 1) - 0.5) * 2.0

    # Soft blue-violet dust, deliberately kept below the disk's exposure.
    dust = np.exp(-((y + 0.18 * np.sin(x * 3.0)) / 0.52) ** 2)
    dust *= 0.5 + 0.5 * np.sin(x * 7.0 + y * 11.0)
    scene = np.empty((height, width, 3), dtype=np.float32)
    scene[..., 0] = 0.002 + dust * 0.006
    scene[..., 1] = 0.004 + dust * 0.010
    scene[..., 2] = 0.010 + dust * 0.028

    # Draw stars as additive point sprites, with a few coloured bright stars.
    count = max(40, width * height // 950)
    sx = rng.integers(1, max(width - 1, 2), size=count)
    sy = rng.integers(1, max(height - 1, 2), size=count)
    brightness = rng.power(5.5, size=count) * 0.38 + 0.025
    tint = rng.uniform(0.78, 1.25, size=(count, 3))
    tint[:, 2] *= 1.18
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


def disk_radiance(radius: np.ndarray, redshift: np.ndarray, peak_temperature: float) -> np.ndarray:
    """Thermal thin-disk emission transported by the ray-traced redshift.

    Liouville's theorem gives I_nu / nu^3 as an invariant.  For a blackbody,
    this means the observed colour temperature is g*T and the bolometric
    intensity gains g^4, where g is calculated by the C++ geodesic renderer.
    """
    positive = radius[radius > 0.0]
    inner = float(np.percentile(positive, 1))
    scaled_radius = np.maximum(radius / max(inner, 1e-6), 1.0)
    no_torque = np.maximum(1.0 - np.sqrt(1.0 / scaled_radius), 0.015)
    emitted_temperature = peak_temperature * scaled_radius ** -0.75 * no_torque**0.25
    observed_temperature = np.clip(emitted_temperature * redshift, 900.0, 100_000.0)

    # Interpolate the Planck/CIE table for each pixel and apply the disk's
    # radial flux profile plus the invariant g^4 brightness transport.
    lookup = np.interp(observed_temperature.ravel(), BLACKBODY_TEMPERATURES, np.arange(BLACKBODY_TEMPERATURES.size))
    lower = np.floor(lookup).astype(np.int32)
    upper = np.minimum(lower + 1, BLACKBODY_TEMPERATURES.size - 1)
    mix = (lookup - lower)[:, None]
    colour = (BLACKBODY_RGB[lower] * (1.0 - mix) + BLACKBODY_RGB[upper] * mix).reshape((*radius.shape, 3))
    flux = scaled_radius ** -3.0 * no_torque
    intensity = flux * np.clip(redshift, 0.05, 5.0) ** 4
    return colour * (intensity * 8.0)[..., None]


def cinematic_image(hits: np.ndarray, redshift: np.ndarray, seed: int, peak_temperature: float) -> Image.Image:
    height, width = hits.shape
    escaped = hits < 0.0
    captured = hits == 0.0
    disk = hits > 0.0

    scene = np.zeros((height, width, 3), dtype=np.float32)
    scene[escaped] = stars_and_nebula(height, width, seed)[escaped]
    scene[disk] = disk_radiance(hits, redshift, peak_temperature)[disk]

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
    return Image.fromarray((scene * 255.0 + 0.5).astype(np.uint8), mode="RGB")


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a cinematic image from a ray-tracing hit buffer.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="binary float hit buffer")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="PNG output path")
    parser.add_argument("--seed", type=int, default=23071969, help="fixed star-field seed")
    parser.add_argument("--peak-temperature", type=float, default=12_000.0, help="peak local disk temperature in kelvin")
    args = parser.parse_args()

    hits, redshift = read_hit_buffer(args.input)
    image = cinematic_image(hits, redshift, args.seed, args.peak_temperature)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    image.save(args.output)
    print(f"Saved cinematic render: {args.output} ({image.width} x {image.height})")


if __name__ == "__main__":
    main()
