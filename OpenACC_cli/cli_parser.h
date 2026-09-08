#pragma once

#include <cstdint>
#include <charconv>
#include <cstring>
#include <string>
#include <unordered_map>

enum class Precession
{
    Float,
    Double
};

struct Params {
    Precession prec = Precession::Float;
    uint64_t kepernyoSZELES = 10240;
    uint64_t kepernyoMAGAS  = 5120;
    uint64_t SZELES = 640;
    uint64_t MAGAS  = 320;

    uint64_t ikezd = 0;
    uint64_t jkezd = 0;
    uint64_t iveg = 0;

    // Per-step local error tolerance for the DOPRI5 pair.
    double errormax = 0.0001f;
    // Largest affine step the controller may take.
    double de0 = 0.01f;
    // Step budget per ray; a ray that exceeds it is reported as a failure.
    // This used to be derived from errormax as int(1/errormax), which conflated
    // two unrelated knobs (tightening the tolerance shrank the budget) and
    // overflowed int for errormax below about 1e-9.  The default matches the
    // budget the old expression produced for the most demanding preset the
    // cache_warmer scripts use (errormax = 1e-6).  Lowering it bounds worst-case
    // render time; raising it helps rays that orbit near the photon sphere.
    uint64_t max_steps = 1000000;
    double rs = 0.05f;
    double delta_a = 0.0001f;
    double a = 0.0f;
    double Q = 0.0f;

    double t_0 = 0.0f;
    double r_0 = 1.0f;
    double theta_0 = 1.57f + 0.06f;
    double phi_0 = 0.0f;

    double kepernyo_high = 0.5f;
    double kepernyo_tav  = 0.75f;

    double sugar_ki = 1.01f;
    double gyuru_sugar_kicsi = 0.1f;
    double gyuru_sugar_nagy  = 0.5f;
};

void parse_args(int argc, char* argv[], Params& p);
