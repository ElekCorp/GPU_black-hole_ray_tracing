#!/usr/bin/env bash
# Benchmark harness: christoffel version x integrator x precision x region.
#
# Rebuilds ./main four times (christoffel current/old87 x integrator
# dopri54/RK6 via the BENCH_OLD_CHRISTOFFEL/BENCH_RK6 macros in cuda_ray.h),
# and for each build times SZELES x MAGAS renders across {float,double} x
# {wide,zoomed}, REPS times each, into results.csv.
#
# Must run *inside* an environment with nvc++ and the CUDA runtime/GPU, e.g.:
#   distrobox enter nvcc -- bash benchmarks/run_bench.sh
set -euo pipefail

cd "$(dirname "$0")/.."   # repo root: OpenACC_cli/

OUT_DIR="benchmarks/results"
mkdir -p "$OUT_DIR"
OUT="$OUT_DIR/results.csv"
echo "variant,precision,region,rep,seconds" > "$OUT"

W=512
H=256
A=0.02
Q=0
REPS=7

# Fixed huge virtual canvas (same convention ffi_bridge.cpp's roi_window
# uses), so "region" is just a (cx,cy,roi_w) window into it. EDGE_CY below
# was bisected with find_shadow_edge.py against this exact scene
# (a=0.02, Q=0, rs=0.05 default, r0=1.0, theta0=1.63) to sit exactly on the
# shadow's edge, i.e. the photon-sphere-skimming critical-impact-parameter
# ring, rather than blindly zooming into the (empty) shadow centre. Re-run
# find_shadow_edge.py if the scene parameters below ever change.
KSZ=281474976710656   # 2^48
KMG=140737488355328   # 2^47
EDGE_CY=0.33222071879490656

roi_args() {
  local w="$1"
  python3 - "$w" "$EDGE_CY" "$KSZ" "$KMG" <<'EOF'
import sys
w, cy, sz, mg = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
half = w/2
cx = 0.5
cy = min(max(cy, half), 1-half) if w < 1.0 else 0.5
ik = round((cx-half)*sz); iv = round((cx+half)*sz); jk = round((cy-half)*mg)
print(f"{ik} {jk} {iv}")
EOF
}

read -r WIDE_IK WIDE_JK WIDE_IV <<<"$(roi_args 1.0)"
read -r ZOOM_IK ZOOM_JK ZOOM_IV <<<"$(roi_args 1e-6)"

declare -A VARIANTS=(
  [current]=""
  [old87]="-DBENCH_OLD_CHRISTOFFEL"
  [rk6]="-DBENCH_RK6"
  [old87_rk6]="-DBENCH_OLD_CHRISTOFFEL -DBENCH_RK6"
)

for variant in current old87 rk6 old87_rk6; do
  defs="${VARIANTS[$variant]}"
  echo "=== building variant=$variant  defs=[$defs] ===" >&2
  make clean >/dev/null 2>&1
  make -j"$(nproc)" COMMON_FLAGS="-std=c++17 -Wall $defs" >"$OUT_DIR/build_${variant}.log" 2>&1
  if [ ! -x ./main ]; then
    echo "BUILD FAILED for $variant, see $OUT_DIR/build_${variant}.log" >&2
    tail -40 "$OUT_DIR/build_${variant}.log" >&2
    exit 1
  fi

  for prec in float double; do
    precflag="--$prec"
    for region in wide zoomed; do
      if [ "$region" = wide ]; then
        regionargs="--ikezd $WIDE_IK --jkezd $WIDE_JK --iveg $WIDE_IV"
      else
        regionargs="--ikezd $ZOOM_IK --jkezd $ZOOM_JK --iveg $ZOOM_IV"
      fi
      for rep in $(seq 1 $REPS); do
        t=$(./main --SZELES $W --MAGAS $H --kepernyoSZELES $KSZ --kepernyoMAGAS $KMG \
              $regionargs --a $A --Q $Q $precflag \
              2>/dev/null | grep "Time difference" | sed -E 's/.*= ([0-9.]+)\[s\]/\1/')
        echo "$variant,$prec,$region,$rep,$t" | tee -a "$OUT" >&2
      done
    done
  done
done

echo "=== done, results in $OUT ===" >&2
echo "=== rebuilding the plain (unbenchmarked) binary before you leave ===" >&2
make clean >/dev/null 2>&1
make -j"$(nproc)" >/dev/null 2>&1
