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

// Renders one frame into `out`, a caller-allocated buffer of
// RAY_CHANNELS*szeles*magas floats. Layout matches szinsaver.h's datasaver_T
// (the same one cli_imagemaker.read_hit_buffer already parses):
// out[RAY_CHANNELS*(i*magas+j)+c] for pixel (i,j); see cuda_ray.h's
// RAY_CHANNELS for what each channel c holds.
//
// Returns 0 on success, -1 if szeles/magas don't keep the required 2:1
// aspect ratio (mirrors main.cpp's kepernyoSZELES/kepernyoMAGAS check).
// omega_x/y/z is the axis-angle camera-orientation vector ijk_to_vec_zoom
// expects (direction = rotation axis, magnitude = rotation angle in
// radians), applied to the local screen-space ray directions before the
// position-based tetrad transform. main.cpp always hardcodes this to
// (0, pi, 0) - the fixed 180-degree flip that makes the default view face
// the hole; callers that want that exact default should still pass that
// vector explicitly.
int render_frame_f32(
    double r0, double theta0, double phi0,
    double a, double Q, double rs,
    double errormax, double de0,
    double omega_x, double omega_y, double omega_z,
    uint64_t szeles, uint64_t magas,
    float* out)
{
    // Matches cli_parser.h's Params defaults for the fields this UI never
    // overrides (kepernyoSZELES/kepernyoMAGAS and the screen/collision geometry).
    uint64_t const SZELESregi = 10240, MAGASregi = 5120;
    if (SZELESregi * magas != MAGASregi * szeles) return -1;

    float const x[D] = { 0.0f, float(r0), float(theta0), float(phi0) };
    float const Omega[D - 1] = { float(omega_x), float(omega_y), float(omega_z) };

    ray_step_T<float>(out, szeles, magas, x, Omega, float(a), float(Q), float(rs),
                       float(errormax), float(de0),
                       /*kepernyo_high*/ 0.5f, /*kepernyo_tav*/ 0.75f,
                       /*sugar_ki*/ 1.01f, /*gyuru_sugar_kicsi*/ 0.1f, /*gyuru_sugar_nagy*/ 0.5f,
                       SZELESregi, MAGASregi, /*ikezd*/ 0, /*jkezd*/ 0, /*iveg*/ SZELESregi);
    return 0;
}

} // extern "C"
