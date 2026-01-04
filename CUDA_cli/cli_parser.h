#pragma once

#include <string>
#include <unordered_map>

enum class Precession
{
    Float,
    Double
};

struct Params {
    Precession prec = Precession::Float;
    int kepernyoSZELES = 10240;
    int kepernyoMAGAS  = 5120;
    int SZELES = 640;
    int MAGAS  = 320;

    double errormax = 0.001f;
    double de0 = 0.01f;
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
