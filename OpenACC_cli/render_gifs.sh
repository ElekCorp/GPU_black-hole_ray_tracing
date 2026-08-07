#!/usr/bin/env bash
#
# Build the renderer with make, then render every GIF at a tight integrator
# tolerance. Run from anywhere; it works in its own directory.
#
# WHAT "TIGHTER" MEANS HERE
#
# errormax is the DOPRI5 controller's local error tolerance, and the movie
# scripts default it to web_server's interactive settle tier (1e-3) for the
# flyaround - fine for dragging a camera around, loose for something you are
# going to keep. This script drops it to 1e-8 everywhere.
#
# 1e-8 rather than something smaller because dopri54_step clamps the tolerance
# it actually uses to fmax(errormax, 1e-8): below that the trajectories stop
# changing. What keeps shrinking is the per-ray step BUDGET, which cuda_ray.h
# sets to int(1/errormax) and does not clamp - so a smaller value only buys
# rays more room to keep integrating, which matters near the photon ring where
# they wind many times. Pass --errormax to go lower if you want that budget;
# expect it to cost time and change nothing else.
#
# de0 is deliberately NOT tightened. It is the step CEILING, not the tolerance,
# so shrinking it forces small steps everywhere rather than only where the
# geometry needs them - all cost, little accuracy. The adaptive controller
# already shortens steps where it must.
#
# RESUMABLE
#
# Each GIF is skipped if its output already exists, so an interrupted run picks
# up where it stopped rather than starting over. Use --force to re-render.
# This is also what makes vast_render.sh's preemption handling work.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

ERRORMAX=1e-8
WIDTH=640
FRAMES=48
OUTDIR=web_images
FORCE=0
MOVIES=all
SPIN=0
CHARGE=0

usage() {
    cat <<'USAGE'
usage: render_gifs.sh [options]

  --errormax VALUE   integrator tolerance (default 1e-8, the tightest
                     dopri54_step responds to; smaller only buys step budget)
  --width N          frame width; height is half (default 640). fp64_showdown
                     uses half this per panel, so its output stays this wide
  --frames N         frames per movie (default 48)
  --out DIR          output directory (default web_images)
  --movies LIST      comma-separated: inclination,dolly,ring,fp64 (default all)
  --a VALUE          black hole spin (default 0). The extremal bound is rs/2 =
                     0.025, and the renderer clamps to 0.98 of it, so 0.0245 is
                     as fast as it spins. Non-zero a/Q add a suffix to every
                     output name, so spun renders never overwrite or get
                     mistaken for the Schwarzschild ones
  --Q VALUE          black hole charge (default 0), same units and bound
  --force            re-render even if the output already exists
  --skip-build       do not run make first
  -h, --help         this

Examples:
  ./render_gifs.sh                             # everything, 640 wide, errormax 1e-8
  ./render_gifs.sh --movies fp64 --width 960
  ./render_gifs.sh --a 0.0245                  # near-extremal Kerr, all movies
  ./render_gifs.sh --frames 96 --errormax 1e-9 # more step budget per ray
USAGE
}

SKIP_BUILD=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --errormax) ERRORMAX="$2"; shift 2 ;;
        --width)    WIDTH="$2"; shift 2 ;;
        --frames)   FRAMES="$2"; shift 2 ;;
        --out)      OUTDIR="$2"; shift 2 ;;
        --movies)   MOVIES="$2"; shift 2 ;;
        --a)        SPIN="$2"; shift 2 ;;
        --Q)        CHARGE="$2"; shift 2 ;;
        --force)    FORCE=1; shift ;;
        --skip-build) SKIP_BUILD=1; shift ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

# --------------------------------------------------------------------------
# Interpreter
# --------------------------------------------------------------------------
# .venv first because that is what exists on the developer's machine; then a
# pixi environment, which is how the linux-64 box gets its dependencies; then
# whatever python3 is on PATH. Checked for the imports rather than assumed,
# since a python3 that cannot import PIL fails several minutes into the first
# render instead of here.
pick_python() {
    local candidates=(.venv/bin/python)
    if command -v pixi >/dev/null 2>&1; then
        candidates+=("pixi-run")
    fi
    candidates+=(python3 python)

    for candidate in "${candidates[@]}"; do
        if [[ "$candidate" == "pixi-run" ]]; then
            if pixi run python -c 'import numpy, PIL, tqdm' >/dev/null 2>&1; then
                echo "pixi run python"; return 0
            fi
            continue
        fi
        if [[ -x "$candidate" ]] || command -v "$candidate" >/dev/null 2>&1; then
            if "$candidate" -c 'import numpy, PIL, tqdm' >/dev/null 2>&1; then
                echo "$candidate"; return 0
            fi
        fi
    done
    return 1
}

if ! PYTHON=$(pick_python); then
    echo "error: no python with numpy, Pillow and tqdm available." >&2
    echo "  linux:  pixi install   (pixi.toml declares them)" >&2
    echo "  or:     python3 -m venv .venv && .venv/bin/pip install numpy pillow tqdm fastapi" >&2
    exit 1
fi
echo "python:   $PYTHON"

# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------
# The movie scripts dlopen libblackhole.so through render_worker, so that is
# the target that matters; ./main is the standalone CLI and is not needed here.
if [[ $SKIP_BUILD -eq 0 ]]; then
    echo "building libblackhole.so ..."
    make libblackhole.so
else
    echo "skipping build (--skip-build)"
fi
if [[ ! -f libblackhole.so ]]; then
    echo "error: libblackhole.so missing after build" >&2
    exit 1
fi

mkdir -p "$OUTDIR"
PANEL=$(( WIDTH / 2 ))   # fp64_showdown draws two panels side by side

# Scene parameters go to every movie, and a non-default one renames the outputs.
# Without the rename, rendering a spun set into a directory that already holds
# the Schwarzschild set would skip every movie as "already exists" and quietly
# hand back the wrong hole - the skip-if-exists check is on the filename alone.
# awk does the comparison because these are decimals and bash arithmetic is not.
nonzero() { awk -v v="$1" 'BEGIN { exit !(v + 0 != 0) }'; }
SCENE_ARGS=()
SUFFIX=""
if nonzero "$SPIN";   then SCENE_ARGS+=(--a "$SPIN");   SUFFIX="${SUFFIX}_a${SPIN}"; fi
if nonzero "$CHARGE"; then SCENE_ARGS+=(--Q "$CHARGE"); SUFFIX="${SUFFIX}_Q${CHARGE}"; fi

echo "settings: errormax=$ERRORMAX width=$WIDTH frames=$FRAMES out=$OUTDIR"
if [[ -n "$SUFFIX" ]]; then
    echo "scene:    a=$SPIN Q=$CHARGE -> outputs suffixed '$SUFFIX'"
fi
echo

wants() {
    [[ "$MOVIES" == "all" ]] && return 0
    [[ ",$MOVIES," == *",$1,"* ]]
}

# render <name> <output-file> <script> [args...]
render() {
    local name="$1" output="$2"; shift 2
    if ! wants "$name"; then
        echo "-- $name: not selected, skipping"
        return 0
    fi
    if [[ -f "$output" && $FORCE -eq 0 ]]; then
        echo "-- $name: $output already exists, skipping (--force to redo)"
        return 0
    fi
    echo "-- $name -> $output"
    local started=$SECONDS
    # Rendered to a temporary name and moved into place only on success, so an
    # interrupted run never leaves a half-written GIF that the skip-if-exists
    # check above would then treat as done. The suffix goes BEFORE the
    # extension: Pillow picks its output format from the extension, so a
    # trailing ".partial" makes it refuse to write the file at all.
    local tmp="${output%.gif}.partial.gif"
    $PYTHON "$@" ${SCENE_ARGS[@]+"${SCENE_ARGS[@]}"} -o "$tmp"
    mv "$tmp" "$output"
    echo "   done in $((SECONDS - started))s"
    echo
}

render inclination "$OUTDIR/flyaround_inclination${SUFFIX}.gif" orbit_flyaround.py \
    --mode inclination --frames "$FRAMES" --width "$WIDTH" \
    --errormax "$ERRORMAX" --ms-per-frame 70

render dolly "$OUTDIR/flyaround_dolly${SUFFIX}.gif" orbit_flyaround.py \
    --mode dolly --frames "$FRAMES" --width "$WIDTH" \
    --errormax "$ERRORMAX" --ms-per-frame 70

render ring "$OUTDIR/photon_ring_zoom${SUFFIX}.gif" photon_ring_zoom.py \
    --frames "$FRAMES" --width "$WIDTH" --errormax "$ERRORMAX"

render fp64 "$OUTDIR/fp64_showdown${SUFFIX}.gif" fp64_showdown.py \
    --frames "$FRAMES" --width "$PANEL" \
    --errormax "$ERRORMAX" --csv "$OUTDIR/fp64_showdown${SUFFIX}.csv"

echo "all requested movies present in $OUTDIR:"
ls -lh "$OUTDIR"/*.gif 2>/dev/null || echo "  (none)"
