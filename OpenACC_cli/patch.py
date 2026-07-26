import sys

with open("cuda_ray.h", "r") as f:
    content = f.read()

# Chunk 1
c1_target = """template <class FP>
inline void step_size(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP& de);

template <class FP>
inline FP ijk_to_vec_mink_zoom(uint64_t const i, uint64_t const j, uint64_t const k, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg, kerr_black_hole<FP> const& hole);"""

c1_repl = """// step_size removed

template <class FP>
inline FP ijk_to_vec_mink_zoom(uint64_t const i, uint64_t const j, uint64_t const k, uint64_t const SZELESregi, uint64_t const MAGASregi, uint64_t const ikezd, uint64_t const jkezd, uint64_t const iveg, kerr_black_hole<FP> const& hole);"""

content = content.replace(c1_target, c1_repl)

# Chunk 2
c2_target = """template <class FP>
inline void RK6(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de);"""

c2_repl = """template <class FP>
inline FP RKDP45_core(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ x_new, FP* const __restrict__ v_new, FP const de);"""

content = content.replace(c2_target, c2_repl)

# Chunk 3
c3_target = """template <class FP>
inline void step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)//adaptiv step size
{
    //RK38(hole, x, v, de);//RK38 vagy 6
    RK6(hole, x, v, de);
    step_size(hole, x, v, de);
}

template <class FP>
inline void step_size(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)
{
    FP ch[D];
    FP de0 = hole.de0;

    christoffel(hole, x, v, ch);

    FP err = hole.errormax;
    FP sum = 0.0;
    for (size_t i = 0; i < D; ++i)
    {
        sum += fabs(ch[i]);
    }

    de = sqrt(err / sum);

    if (de > de0)
    {
        de = de0;
    }
    else if (isnan(de))
    {
        de = de0 / 10;
    }
    else
    {

    }

}"""

c3_repl = """template <class FP>
inline void step(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)
{
    bool accepted = false;
    FP x_new[D];
    FP v_new[D];
    FP de_new = de;
    int max_retries = 10;
    int retries = 0;

    while (!accepted && retries < max_retries)
    {
        FP err = RKDP45_core(hole, x, v, x_new, v_new, de_new);
        
        if (err < 1e-20) err = 1e-20;
        
        FP factor = 0.9 * pow(hole.errormax / err, 0.2);
        
        if (factor > 5.0) factor = 5.0;
        if (factor < 0.1) factor = 0.1;

        if (err <= hole.errormax) {
            accepted = true;
        }

        FP de_next = de_new * factor;
        if (de_next > hole.de0) de_next = hole.de0;

        if (accepted) {
            for (int i = 0; i < D; ++i) {
                x[i] = x_new[i];
                v[i] = v_new[i];
            }
            de = de_next; 
        } else {
            de_new = de_next; 
            retries++;
        }
    }
    
    if (!accepted) {
        for (int i = 0; i < D; ++i) {
            x[i] = x_new[i];
            v[i] = v_new[i];
        }
        de = de_new;
    }
}"""

content = content.replace(c3_target, c3_repl)

# Chunk 4
c4_target = """template <class FP>
inline void RK6(kerr_black_hole<FP> const& hole, FP* const __restrict__ x, FP* const __restrict__ v, FP& de)
{"""

c4_repl = """template <class FP>
inline FP RKDP45_core(kerr_black_hole<FP> const& hole, FP const* const __restrict__ x, FP const* const __restrict__ v, FP* const __restrict__ x_new, FP* const __restrict__ v_new, FP const de)
{"""

content = content.replace(c4_target, c4_repl)


# Chunk 5
c5_target = """    for (int i = 0; i < D; ++i)
    {
        x[i] = x[i] + (kx1[i] * (5179.0 / 57600.0) + kx3[i] * (7571.0 / 16695.0) + kx4[i] * (393.0 / 640.0) - kx5[i] * (92097.0 / 339200.0) + kx6[i] * (187.0 / 2100.0) + kx7[i] * (1.0 / 40.0)) * de;
        v[i] = v[i] + (kv1[i] * (5179.0 / 57600.0) + kv3[i] * (7571.0 / 16695.0) + kv4[i] * (393.0 / 640.0) - kv5[i] * (92097.0 / 339200.0) + kv6[i] * (187.0 / 2100.0) + kv7[i] * (1.0 / 40.0)) * de;
    }
}"""

c5_repl = """    FP err = 0.0;
    for (int i = 0; i < D; ++i)
    {
        x_new[i] = x[i] + (kx1[i] * (35.0 / 384.0) + kx3[i] * (500.0 / 1113.0) + kx4[i] * (125.0 / 192.0) - kx5[i] * (2187.0 / 6784.0) + kx6[i] * (11.0 / 84.0)) * de;
        v_new[i] = v[i] + (kv1[i] * (35.0 / 384.0) + kv3[i] * (500.0 / 1113.0) + kv4[i] * (125.0 / 192.0) - kv5[i] * (2187.0 / 6784.0) + kv6[i] * (11.0 / 84.0)) * de;
        
        FP x_4th = x[i] + (kx1[i] * (5179.0 / 57600.0) + kx3[i] * (7571.0 / 16695.0) + kx4[i] * (393.0 / 640.0) - kx5[i] * (92097.0 / 339200.0) + kx6[i] * (187.0 / 2100.0) + kx7[i] * (1.0 / 40.0)) * de;
        FP v_4th = v[i] + (kv1[i] * (5179.0 / 57600.0) + kv3[i] * (7571.0 / 16695.0) + kv4[i] * (393.0 / 640.0) - kv5[i] * (92097.0 / 339200.0) + kv6[i] * (187.0 / 2100.0) + kv7[i] * (1.0 / 40.0)) * de;

        FP diff_x = x_new[i] - x_4th;
        FP diff_v = v_new[i] - v_4th;
        err += diff_x * diff_x + diff_v * diff_v;
    }
    return sqrt(err);
}"""

content = content.replace(c5_target, c5_repl)

with open("cuda_ray.h", "w") as f:
    f.write(content)
print("done")
