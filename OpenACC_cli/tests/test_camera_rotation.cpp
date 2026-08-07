// Checks that the camera-orientation rotation inside ijk_to_vec_zoom is an
// actual rotation, by calling the shipped function rather than restating it.
//
// ijk_to_vec_zoom does two things in sequence: it rotates the local screen-space
// ray direction by the axis-angle vector Omega, then maps that through the
// position-dependent tetrad into Boyer-Lindquist components.  Only the first is
// under test, so the tetrad is inverted back out: at a = 0 it is diagonal with
// closed-form entries (see the transform at the end of ijk_to_vec_zoom, with
// rho = r_0 and delta = r_0^2 - rs*r_0), which is what makes the rotated vector
// recoverable from the function's actual output.
//
// Two independent properties are asserted:
//   1. the recovered vector matches textbook Rodrigues applied to the same
//      input, and
//   2. the transform preserves the input's norm - the direct statement that it
//      rotates rather than shears, and the one the old code violated.
//
// Build and run:  make -C .. test

#include <cmath>
#include <cstdint>
#include <cstdio>

#include "../black_hole.h"
#include "../cuda_ray.h"

namespace {

struct Vec3 { double v[3]; };

int failures = 0;

void check(bool ok, char const* what, double measured, double tolerance)
{
    if (ok) return;
    ++failures;
    printf("  FAIL  %-46s  %.3e (tolerance %.1e)\n", what, measured, tolerance);
}

// The pixel -> unit direction map ijk_to_vec_mink_zoom builds, before any
// rotation.  Trivial enough to restate; it is the test's input, not its
// subject.
Vec3 screen_direction(uint64_t i, uint64_t j, uint64_t szeles, uint64_t magas,
                      double kepernyo_high, double kepernyo_tav)
{
    double const y = (kepernyo_high / double(magas)) * (double(magas) / 2 - (double(j) + 0.5));
    double const z = (kepernyo_high / double(magas)) * ((double(i) + 0.5) - double(szeles) / 2);
    double const norm = std::sqrt(kepernyo_tav * kepernyo_tav + y * y + z * z);
    return Vec3{{kepernyo_tav / norm, y / norm, z / norm}};
}

// Textbook Rodrigues, as the reference the shipped code has to reproduce.
Vec3 rodrigues(Vec3 const& u, double phi, Vec3 const& x)
{
    double const c = std::cos(phi), s = std::sin(phi), t = 1.0 - c;
    double const R[3][3] = {
        { c + u.v[0] * u.v[0] * t,        u.v[0] * u.v[1] * t - u.v[2] * s, u.v[0] * u.v[2] * t + u.v[1] * s },
        { u.v[1] * u.v[0] * t + u.v[2] * s, c + u.v[1] * u.v[1] * t,        u.v[1] * u.v[2] * t - u.v[0] * s },
        { u.v[2] * u.v[0] * t - u.v[1] * s, u.v[2] * u.v[1] * t + u.v[0] * s, c + u.v[2] * u.v[2] * t },
    };
    Vec3 out{{0, 0, 0}};
    for (int r = 0; r < 3; ++r)
        for (int k = 0; k < 3; ++k) out.v[r] += R[r][k] * x.v[k];
    return out;
}

double norm3(Vec3 const& x)
{
    return std::sqrt(x.v[0] * x.v[0] + x.v[1] * x.v[1] + x.v[2] * x.v[2]);
}

// One (axis, angle, pixel) sample: run it through the shipped ijk_to_vec_zoom,
// undo the tetrad, and compare against the reference.
void probe(Vec3 const& axis, double phi, uint64_t i, uint64_t j,
           double& direction_error, double& norm_error)
{
    uint64_t const szeles = 64, magas = 32;
    double const rs = 0.05, a = 0.0, Q = 0.0;
    double const r_0 = 1.0, theta_0 = 1.63;
    double const kepernyo_high = 0.5, kepernyo_tav = 0.75;

    double const x[4] = { 0.0, r_0, theta_0, 0.0 };
    double const Omega[3] = { axis.v[0] * phi, axis.v[1] * phi, axis.v[2] * phi };
    kerr_black_hole<double> hole(szeles, magas, x, Omega, a, Q, rs,
                                 /*errormax*/ 1e-4, /*de0*/ 1e-2,
                                 kepernyo_high, kepernyo_tav,
                                 /*sugar_ki*/ 1.01, /*kicsi*/ 0.1, /*nagy*/ 0.5);

    Vec3 traced{{0, 0, 0}};
    for (int k = 0; k < 3; ++k)
        traced.v[k] = ijk_to_vec_zoom<double>(i, j, k + 1, hole, szeles, magas, 0, 0, szeles);

    // Invert the a = 0 tetrad: rho = r_0, delta = r_0^2 - rs*r_0, and the
    // spatial block is diagonal with these three positive entries.
    double const rho = r_0;
    double const delta = r_0 * r_0 - rs * r_0;
    double const scale[3] = {
        std::sqrt(delta) / rho,
        1.0 / rho,
        rho / (std::sin(theta_0) * r_0 * r_0),
    };
    Vec3 recovered{{0, 0, 0}};
    for (int k = 0; k < 3; ++k) recovered.v[k] = traced.v[k] / scale[k];

    Vec3 const input = screen_direction(i, j, szeles, magas, kepernyo_high, kepernyo_tav);
    Vec3 const expected = rodrigues(axis, phi, input);

    direction_error = 0.0;
    for (int k = 0; k < 3; ++k)
        direction_error = std::fmax(direction_error, std::fabs(recovered.v[k] - expected.v[k]));
    norm_error = std::fabs(norm3(recovered) - norm3(input));
}

}  // namespace

int main()
{
    Vec3 axes[] = {
        {{0.0, 1.0, 0.0}},                       // the default view's axis
        {{0.0, 0.0, 1.0}},
        {{0.3, 0.5, 0.81}},
        {{-0.6, 0.2, -0.77}},
    };
    for (Vec3& u : axes) {
        double const n = norm3(u);
        for (int k = 0; k < 3; ++k) u.v[k] /= n;
    }

    double const angles[] = { 0.01, 0.4, 1.2, 2.0, M_PI / 2, M_PI, 2.5 * M_PI };
    uint64_t const pixels[][2] = { {0, 0}, {17, 5}, {31, 16}, {63, 31} };

    double const tolerance = 1e-14;
    printf("ijk_to_vec_zoom camera rotation\n");

    for (Vec3 const& u : axes) {
        for (double phi : angles) {
            double worst_direction = 0.0, worst_norm = 0.0;
            for (auto const& p : pixels) {
                double direction_error, norm_error;
                probe(u, phi, p[0], p[1], direction_error, norm_error);
                worst_direction = std::fmax(worst_direction, direction_error);
                worst_norm = std::fmax(worst_norm, norm_error);
            }
            char label[96];
            snprintf(label, sizeof label, "axis (%.2f %.2f %.2f) angle %.4f",
                     u.v[0], u.v[1], u.v[2], phi);
            printf("  %-46s  direction %.2e  norm %.2e\n", label, worst_direction, worst_norm);
            check(worst_direction < tolerance, "  ^ matches Rodrigues", worst_direction, tolerance);
            check(worst_norm < tolerance, "  ^ preserves norm", worst_norm, tolerance);
        }
    }

    // The default view is a 180-degree flip, where sin(phi) == 0 and the old
    // misplaced parentheses made no difference.  Confirm that directly, so the
    // claim that this fix leaves the default render untouched is checked rather
    // than asserted: reproduce the OLD row-2 expression and show it agrees with
    // the corrected one at pi, and disagrees away from it.
    printf("\ndefault 180-degree view is unaffected by the fix\n");
    for (double phi : { M_PI, 1.2 }) {
        Vec3 u = axes[2];
        Vec3 const in = screen_direction(17, 5, 64, 32, 0.5, 0.75);
        double const c = std::cos(phi), s = std::sin(phi);
        // Verbatim from the pre-fix cuda_ray.h line 768.
        double const old_row2 =
            (u.v[0] * u.v[1] * (1 - c + u.v[2] * s)) * in.v[0]
            + (c + u.v[1] * u.v[1] * (1 - c)) * in.v[1]
            + (u.v[1] * u.v[2] * (1 - c + u.v[0] * s)) * in.v[2];
        double const fixed_row2 = rodrigues(u, phi, in).v[1];
        double const gap = std::fabs(old_row2 - fixed_row2);
        printf("  angle %.4f  |old - fixed| = %.2e\n", phi, gap);
        if (phi == M_PI) check(gap < tolerance, "  ^ identical at pi", gap, tolerance);
        else check(gap > 1e-3, "  ^ differs away from pi (bug was real)", gap, 1e-3);
    }

    printf("\n%s\n", failures == 0 ? "PASS" : "FAILED");
    return failures == 0 ? 0 : 1;
}
