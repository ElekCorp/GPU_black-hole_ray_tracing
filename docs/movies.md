# Movie scripts

Four scripts under `OpenACC_cli/` render GIFs by driving the geodesic renderer
through `render_worker`, then shading with `cli_imagemaker`. All of them take
their scene constants, accuracy tiers and zoom/precision policy from
`web_server.py`, so a movie always shows the same hole as the interactive view.

Run them from inside `OpenACC_cli/` (the library path and the default output
paths are relative), with an environment that has numpy, Pillow, fastapi and
tqdm — on macOS that is `OpenACC_cli/.venv`:

```bash
cd OpenACC_cli && .venv/bin/python orbit_flyaround.py --help
```

| Script | What it renders |
| --- | --- |
| `orbit_flyaround.py` | Camera moves: inclination sweep, azimuth orbit, gaze cone, dolly. Optional simultaneous true zoom. |
| `fp64_showdown.py` | fp32 and fp64 side by side, zooming until fp32 fails, with the failure measured per frame. |
| `photon_ring_zoom.py` | Zoom into the photon ring one subring at a time. |
| `photon_ring_zoom_tqdm.py` | The same, with a progress bar. |

`movie_common.py` holds the shared plumbing (tracing, shading, shadow-edge
bisection, captions, GIF encoding) for the first two. The `photon_ring_zoom*`
pair predates it and still carries its own copies.

## Which move is worth rendering

`orbit_flyaround.py --mode inclination` is the one to reach for. The camera
swings from looking down on the disk to seeing it edge-on; the far side of the
disk is lensed up and over the shadow and the whole thing closes into the
familiar band with a bright ring around the silhouette.

`--mode azimuth` is a pretty loop but not a new view. Kerr–Newman is
axisymmetric and the disk is an axisymmetric annulus, so orbiting the spin axis
leaves the image geometry **invariant** — only the disk's turbulence texture
rotates. That is the same symmetry that makes `phi0` free to precompute; see
[precompute-investigation.md](precompute-investigation.md) §2.

`--mode look` swings the camera's gaze in a cone, composing a pan and a tilt at
once. That is precisely the case the Rodrigues rotation in `cuda_ray.h` used to
get wrong — a pan alone or a tilt alone both landed in the degenerate set where
the old expression was accidentally correct. Covered now by
`tests/test_camera_rotation.cpp`.

## Two constraints worth knowing before you plan a shot

**The camera cannot back away.** `ffi_bridge.cpp` fixes the outer integration
sphere at `sugar_ki = 1.01`, and `ray_step_T` ends any ray already outside it.
A camera beyond that radius has every ray terminate on its first step and the
frame comes back as a blank starfield — measured, `r0 = 1.01` renders normally
and `r0 = 1.02` is empty. So `--mode dolly` only flies *in*. (`web_server.py`'s
`MAX_R0 = 6.0` does not know this and can dolly the interactive view to a blank
screen.) The scripts refuse such an `r0` rather than rendering nothing.

**A zoom does not track.** `--cx`/`--cy` are fixed for the whole movie. Fine on
their own or with `look`; a trap with `inclination` and `dolly`, which carry the
subject out of the window. Past roughly 30× the end of the movie is a flat
single-colour frame. The script counts featureless frames and says so.

## What the fp64 demo actually establishes

`fp64_showdown.py` zooms both precisions onto the shadow edge and scores fp32
against fp64 frame by frame. It separates two failures that are usually
conflated, because they arrive two decades apart and look nothing alike:

1. **Chaotic amplification**, from around a thousand ×, and this is the main
   event. A ray skimming the photon sphere is exponentially sensitive to where
   it started — the same instability that produces the subring cascade — so it
   amplifies fp32 rounding as efficiently as a real difference in aim. The rays
   stay distinct and finite and are individually wrong; the frame goes to
   salt-and-pepper noise. Measured at 240px panels: agreement with fp64 falls to
   84% at 1,350×, 64% at 2,462× and 40% at 8,185×.
2. **Rays going non-finite**, from around 1e5 ×. Integrations start failing
   outright; `ray_step_T` takes its NaN/inf bail-out and `sky_direction`
   substitutes one fixed sentinel direction, which reads on screen as flat
   patches. fp32 reaches 13.8% of the frame at 1e6× against fp64's 7.1%.

**What is *not* happening**, despite being the obvious guess: adjacent pixels'
launch directions rounding to the same float, so one geodesic gets traced where
the frame asked for two. Measured directly — excluding the sentinel above, over
24 frames from 1× to 1e6×, **fp32 collapsed exactly zero neighbouring pairs and
fp64 collapsed three** (of 12,045) in one frame. `ijk_to_vec_mink_zoom` computes
the launch direction in double even for the f32 trace, which is a large part of
why.

This distinction is easy to get wrong, and worth guarding: counting the sentinel
pixels as collapse makes fp32 look like it is losing resolution — the number
climbs smoothly to 26% at 1e6× — when in fact those are broken rays. The script
reports the two separately for that reason.

Past ~1e6 × fp64 has its own ceiling for a third reason: `dopri54_step` clamps
its step tolerance to `fmax(errormax, 1e-8)`. So the honest claim is bounded —
fp64 buys roughly two decades of usable zoom, after which the integrator rather
than the number format is what needs fixing.

`FLOAT_ZOOM_LIMIT = 512` in `web_server.py` is therefore a reasonable switch
point, but it is guarding against amplified rounding error, not against lost
resolution.
