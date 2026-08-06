// Thin extern "C" entry point for calling the ray tracer as a library instead
// of spawning ./main as a subprocess per frame. Reuses ray_step_T unchanged -
// this file duplicates none of the ray-marching logic, only the small
// argument-assembly main.cpp already does around it (see makeframe_T).
//
// Built as a shared library (see the Makefile's libblackhole.so target) so a
// long-lived process can dlopen it once and keep the CUDA context warm across
// many renders, instead of paying CUDA context init on every single frame.

#include <cstdint>
#include <math.h>

#include "black_hole.h"
#include "cuda_ray.h"

extern "C" {

// Floats per pixel render_frame_f32 writes, so a caller sizes its buffer from
// the library it actually loaded instead of hardcoding cuda_ray.h's constant.
int render_frame_channels(void)
{
    return RAY_CHANNELS;
}

// The virtual screen grid the zoom window is expressed in. ijk_to_vec_zoom
// wants the window as integer indices (ikezd/jkezd/iveg) into a grid of
// SZELESregi x MAGASregi, but the ray directions it builds depend only on the
// *ratio* of those indices to the grid - so the grid size sets nothing but how
// finely the window can be positioned, and a huge one costs nothing.
//
// 2^48 balances the two ways the grid can run out. Too small and the window's
// edges quantize visibly (placement granularity is Z*szeles/SZELESregi pixels
// at magnification Z). Too large and the double arithmetic in
// ijk_to_vec_mink_zoom loses its own resolution, since the ULP there scales
// with SZELESregi: at 2^48 it is 2^-5 of a grid cell, still thousands of
// sub-steps per rendered pixel at a billion times magnification. Both margins
// sit far beyond what the geodesic integrator can resolve.
static uint64_t const ROI_SZELESregi = 1ull << 48;
static uint64_t const ROI_MAGASregi = ROI_SZELESregi / 2;

// Fraction of the unzoomed frame -> index into a grid of `extent` cells.
static uint64_t roi_index(double fraction, uint64_t extent)
{
    if (!(fraction > 0.0)) return 0;  // the negation also catches NaN
    if (fraction >= 1.0) return extent;
    return (uint64_t)(fraction * (double)extent + 0.5);
}

// The zoom window, resolved to the integer indices ray_step_T wants. Shared by
// both precisions so they can never disagree about which region they traced.
struct RoiWindow { uint64_t ikezd, jkezd, iveg; };

static int roi_window(double roi_cx, double roi_cy, double roi_w,
                      uint64_t szeles, uint64_t magas, RoiWindow& window)
{
    if (ROI_SZELESregi * magas != ROI_MAGASregi * szeles) return -1;
    if (!(roi_w > 0.0) || roi_w > 1.0) return -2;

    // Keep the window inside the unzoomed frame; ikezd/jkezd are unsigned, and
    // a window hanging off the edge would only trace rays the camera's screen
    // never covered anyway.
    double const half = roi_w / 2.0;
    double const cx = fmin(fmax(roi_cx, half), 1.0 - half);
    double const cy = fmin(fmax(roi_cy, half), 1.0 - half);

    window.ikezd = roi_index(cx - half, ROI_SZELESregi);
    window.iveg = roi_index(cx + half, ROI_SZELESregi);
    window.jkezd = roi_index(cy - half, ROI_MAGASregi);
    // A window narrower than one grid cell has collapsed to a single ray
    // direction - every pixel would trace the same geodesic.
    if (window.iveg <= window.ikezd) return -3;
    return 0;
}

// Renders one frame into `out`, a caller-allocated buffer of
// RAY_CHANNELS*szeles*magas floats. Layout matches szinsaver.h's datasaver_T
// (the same one cli_imagemaker.read_hit_buffer already parses):
// out[RAY_CHANNELS*(i*magas+j)+c] for pixel (i,j); see cuda_ray.h's
// RAY_CHANNELS for what each channel c holds.
//
// omega_x/y/z is the axis-angle camera-orientation vector ijk_to_vec_zoom
// expects (direction = rotation axis, magnitude = rotation angle in
// radians), applied to the local screen-space ray directions before the
// position-based tetrad transform. main.cpp always hardcodes this to
// (0, pi, 0) - the fixed 180-degree flip that makes the default view face
// the hole; callers that want that exact default should still pass that
// vector explicitly.
//
// roi_cx/roi_cy/roi_w select the sub-rectangle of the camera's full field of
// view to trace, in fractions of the unzoomed frame: (roi_cx, roi_cy) is the
// window's centre and roi_w its width. Its height fraction is roi_w too, since
// the window and the full frame keep the same 2:1 aspect. roi_w == 1 with the
// centre at (0.5, 0.5) is the unzoomed view; halving roi_w doubles the
// magnification. This is a true optical zoom, not a crop-and-upscale: the same
// szeles*magas rays are retraced across the narrower solid angle, so the output
// carries detail the unzoomed frame never sampled. The camera does not move -
// r0/theta0/phi0 are untouched.
//
// Returns 0 on success, -1 if szeles/magas don't keep the required 2:1 aspect
// ratio (mirrors main.cpp's kepernyoSZELES/kepernyoMAGAS check), -2 if roi_w is
// not in (0, 1], -3 if the window has collapsed below one grid cell.
//
// This is the SINGLE-PRECISION trace. Measured against a smooth patch of lensed
// sky, it resolves neighbouring rays cleanly to roughly 500x magnification and
// has visibly degenerated by 4000x - adjacent pixels start landing on identical
// geodesics, and what is left of the pixel-to-pixel variation is integration
// error rather than geometry. Tightening errormax does not move that ceiling,
// because what runs out is the float representation of the ray state itself.
// Past it, use render_frame_roi_f64.
int render_frame_roi_f32(
    double r0, double theta0, double phi0,
    double a, double Q, double rs,
    double errormax, double de0,
    double omega_x, double omega_y, double omega_z,
    double roi_cx, double roi_cy, double roi_w,
    uint64_t szeles, uint64_t magas,
    float* out)
{
    RoiWindow window;
    int const status = roi_window(roi_cx, roi_cy, roi_w, szeles, magas, window);
    if (status != 0) return status;

    float const x[D] = { 0.0f, float(r0), float(theta0), float(phi0) };
    float const Omega[D - 1] = { float(omega_x), float(omega_y), float(omega_z) };

    ray_step_T<float>(out, szeles, magas, x, Omega, float(a), float(Q), float(rs),
                       float(errormax), float(de0),
                       /*kepernyo_high*/ 0.5f, /*kepernyo_tav*/ 0.75f,
                       /*sugar_ki*/ 1.01f, /*gyuru_sugar_kicsi*/ 0.1f, /*gyuru_sugar_nagy*/ 0.5f,
                       ROI_SZELESregi, ROI_MAGASregi, window.ikezd, window.jkezd, window.iveg);
    return 0;
}

// The DOUBLE-PRECISION trace, for magnifications the float path can no longer
// resolve. Identical geometry and identical arguments - only the integrator's
// working precision and the output buffer's element type differ.
//
// `out` is doubles, not floats, and that is not incidental: at deep zoom
// neighbouring pixels' escape directions differ by less than a float's ULP, so
// narrowing the result on the way out would throw away exactly the detail this
// path exists to compute.
//
// Expect this to be markedly slower than the f32 path on hardware with a
// reduced fp64 rate (most consumer GPUs run fp64 at 1/32 or 1/64 of fp32), so
// callers should stay on f32 until the zoom actually needs this.
int render_frame_roi_f64(
    double r0, double theta0, double phi0,
    double a, double Q, double rs,
    double errormax, double de0,
    double omega_x, double omega_y, double omega_z,
    double roi_cx, double roi_cy, double roi_w,
    uint64_t szeles, uint64_t magas,
    double* out)
{
    RoiWindow window;
    int const status = roi_window(roi_cx, roi_cy, roi_w, szeles, magas, window);
    if (status != 0) return status;

    double const x[D] = { 0.0, r0, theta0, phi0 };
    double const Omega[D - 1] = { omega_x, omega_y, omega_z };

    ray_step_T<double>(out, szeles, magas, x, Omega, a, Q, rs, errormax, de0,
                       /*kepernyo_high*/ 0.5, /*kepernyo_tav*/ 0.75,
                       /*sugar_ki*/ 1.01, /*gyuru_sugar_kicsi*/ 0.1, /*gyuru_sugar_nagy*/ 0.5,
                       ROI_SZELESregi, ROI_MAGASregi, window.ikezd, window.jkezd, window.iveg);
    return 0;
}

// The unzoomed whole-frame render, i.e. what this entry point did before the
// zoom window existed.
int render_frame_f32(
    double r0, double theta0, double phi0,
    double a, double Q, double rs,
    double errormax, double de0,
    double omega_x, double omega_y, double omega_z,
    uint64_t szeles, uint64_t magas,
    float* out)
{
    return render_frame_roi_f32(r0, theta0, phi0, a, Q, rs, errormax, de0,
                                omega_x, omega_y, omega_z,
                                /*roi_cx*/ 0.5, /*roi_cy*/ 0.5, /*roi_w*/ 1.0,
                                szeles, magas, out);
}

} // extern "C"
