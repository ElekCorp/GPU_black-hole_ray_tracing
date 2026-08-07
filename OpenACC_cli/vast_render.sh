#!/usr/bin/env bash
#
# Render the GIFs on a vast.ai GPU instance. Wraps render_gifs.sh with the
# things a rented, possibly interruptible box needs: a preflight that fails
# loudly instead of silently producing a CPU build, dependency setup, syncing
# results off the instance as they finish, and optional self-shutdown.
#
# RECOMMENDED IMAGE
#
#   nvcr.io/nvidia/nvhpc:24.5-devel-cuda_multi-ubuntu22.04
#
# The Makefile compiles the ray_step kernels with OpenACC and only offloads to
# the GPU under nvc++ (-acc=gpu). That compiler ships in the NVIDIA HPC SDK, and
# picking an image that already has it is far quicker than installing several
# gigabytes of SDK on every instance. With any other image the Makefile silently
# falls back to g++/clang++ and you get a CPU build - which still renders, just
# orders of magnitude slower, on hardware you are paying GPU rates for. The
# preflight below refuses that by default; --allow-cpu overrides it.
#
# INTERRUPTIBLE INSTANCES
#
# A preempted instance can lose everything on local disk. Two things address it:
#
#   - render_gifs.sh skips any GIF whose output already exists, so re-running
#     this script after a restart resumes rather than starting over.
#   - --sync-cmd runs after EACH finished GIF, not once at the end, so at most
#     one movie's work is ever at risk. Nothing is synced by default, because
#     where your results should go is not something this script can guess.
#
# Usage (on the instance):
#   ./vast_render.sh --check
#   ./vast_render.sh --width 1280 --frames 96
#   ./vast_render.sh --sync-cmd 'rclone copy {} remote:blackhole/' --shutdown
#
# NOTE: written against vast.ai's documented environment but not exercised on a
# live instance from here - there is no GPU on the machine it was authored on.
# The preflight is deliberately thorough for that reason: run --check first and
# it will tell you what is missing before anything long-running starts.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

OUTDIR=""
SYNC_CMD=""
ERRORMAX=1e-8
WIDTH=1280
FRAMES=96
MOVIES=all
CHECK_ONLY=0
ALLOW_CPU=0
INSTALL_NVHPC=0
SHUTDOWN=0
REBUILD=0
DEPS=pip

usage() {
    cat <<'USAGE'
usage: vast_render.sh [options]

  --outdir DIR       where GIFs go (default: /workspace/blackhole_out if
                     /workspace exists, else ./web_images)
  --sync-cmd 'CMD'   run after each finished GIF; {} is replaced by the output
                     directory. e.g. 'rclone copy {} remote:bh/'
                     Strongly recommended on interruptible instances.
  --errormax VALUE   integrator tolerance (default 1e-8)
  --width N          frame width (default 1280)
  --frames N         frames per movie (default 96)
  --movies LIST      inclination,dolly,ring,fp64 (default all)
  --deps pip|pixi|none   how to get numpy/Pillow/tqdm/fastapi (default pip)
  --install-nvhpc    apt-install the NVIDIA HPC SDK if nvc++ is missing
                     (multi-GB; prefer an image that already has it)
  --allow-cpu        proceed without nvc++, accepting a CPU build
  --rebuild          make clean first (only needed when switching compilers)
  --shutdown         power the instance off when finished, to stop billing
  --check            run the preflight and exit
  -h, --help         this
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)   OUTDIR="$2"; shift 2 ;;
        --sync-cmd) SYNC_CMD="$2"; shift 2 ;;
        --errormax) ERRORMAX="$2"; shift 2 ;;
        --width)    WIDTH="$2"; shift 2 ;;
        --frames)   FRAMES="$2"; shift 2 ;;
        --movies)   MOVIES="$2"; shift 2 ;;
        --deps)     DEPS="$2"; shift 2 ;;
        --install-nvhpc) INSTALL_NVHPC=1; shift ;;
        --allow-cpu)     ALLOW_CPU=1; shift ;;
        --shutdown)      SHUTDOWN=1; shift ;;
        --rebuild)       REBUILD=1; shift ;;
        --check)         CHECK_ONLY=1; shift ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "$OUTDIR" ]]; then
    if [[ -d /workspace ]]; then OUTDIR=/workspace/blackhole_out; else OUTDIR=web_images; fi
fi

# vast.ai containers usually run as root with no sudo installed at all, so
# neither "use sudo" nor "never use sudo" is safe to assume.
SUDO=""
if [[ "$(id -u)" -ne 0 ]] && command -v sudo >/dev/null 2>&1; then SUDO="sudo"; fi

say() { printf '%s\n' "$*"; }

# The current phase, so a failure names where it happened. Worth the two extra
# lines: when this runs unattended on a rented box, "it failed" and "it failed
# installing python" are the difference between one diagnostic run and three.
CURRENT_STAGE="startup"
rule() { CURRENT_STAGE="$*"; printf -- '---- %s\n' "$*"; }
trap 'printf "\n!!!! FAILED during: %s\n" "$CURRENT_STAGE" >&2' ERR

# --------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------
problems=0
note_problem() { say "  PROBLEM: $*"; problems=$((problems + 1)); }

rule "preflight"

if command -v nvidia-smi >/dev/null 2>&1; then
    gpu=$(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null | head -1 || true)
    say "  gpu:      ${gpu:-present but nvidia-smi query failed}"
else
    note_problem "nvidia-smi not found - this does not look like a GPU instance."
fi

if command -v nvc++ >/dev/null 2>&1; then
    say "  nvc++:    $(nvc++ --version 2>/dev/null | head -2 | tail -1 | tr -s ' ')"
elif [[ $INSTALL_NVHPC -eq 1 ]]; then
    say "  nvc++:    missing, will install the NVIDIA HPC SDK"
elif [[ $ALLOW_CPU -eq 1 ]]; then
    say "  nvc++:    missing - continuing with a CPU build (--allow-cpu)"
else
    note_problem "nvc++ not found. Without it the Makefile falls back to g++/clang++
           and builds a CPU-only library, which renders far slower on a GPU you
           are paying for. Either use an image that ships the HPC SDK
           (nvcr.io/nvidia/nvhpc:24.5-devel-cuda_multi-ubuntu22.04), pass
           --install-nvhpc, or pass --allow-cpu if you meant to use the CPU."
fi

for tool in make git; do
    command -v "$tool" >/dev/null 2>&1 && say "  $tool:     $(command -v $tool)" \
        || note_problem "$tool not found (apt-get install -y build-essential git)"
done

say "  outdir:   $OUTDIR"
say "  sync:     ${SYNC_CMD:-<none> - results will be LOST if this instance is preempted}"
say "  settings: errormax=$ERRORMAX width=$WIDTH frames=$FRAMES movies=$MOVIES"

if [[ $problems -gt 0 ]]; then
    say ""
    say "preflight found $problems problem(s); not starting."
    exit 1
fi
say "  preflight OK"
say ""

if [[ $CHECK_ONLY -eq 1 ]]; then
    say "--check given, stopping here."
    exit 0
fi

# --------------------------------------------------------------------------
# NVIDIA HPC SDK, if asked for
# --------------------------------------------------------------------------
if [[ $INSTALL_NVHPC -eq 1 ]] && ! command -v nvc++ >/dev/null 2>&1; then
    rule "installing NVIDIA HPC SDK (this is several GB and takes a while)"
    $SUDO apt-get update
    $SUDO apt-get install -y curl gnupg ca-certificates
    curl -fsSL https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK \
        | $SUDO gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
    echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' \
        | $SUDO tee /etc/apt/sources.list.d/nvhpc.list >/dev/null
    $SUDO apt-get update
    $SUDO apt-get install -y nvhpc-24-5
    # The SDK installs outside PATH; find the newest compiler bin it left behind.
    hpc_bin=$(ls -d /opt/nvidia/hpc_sdk/Linux_x86_64/*/compilers/bin 2>/dev/null | sort -V | tail -1 || true)
    if [[ -n "$hpc_bin" ]]; then
        export PATH="$hpc_bin:$PATH"
        say "  added to PATH: $hpc_bin"
    fi
    command -v nvc++ >/dev/null 2>&1 || { say "  nvc++ still not on PATH after install"; exit 1; }
fi

# --------------------------------------------------------------------------
# Python dependencies
# --------------------------------------------------------------------------
# Only what the movie scripts actually import. fastapi and pydantic are in the
# list because the scripts import web_server for its scene constants and
# zoom/precision policy, which pulls FastAPI in even though no server is run.
PY_PACKAGES=(numpy pillow tqdm fastapi pydantic)

case "$DEPS" in
    none) rule "skipping dependency setup (--deps none)" ;;
    pixi)
        rule "installing dependencies with pixi"
        command -v pixi >/dev/null 2>&1 || curl -fsSL https://pixi.sh/install.sh | bash
        export PATH="$HOME/.pixi/bin:$PATH"
        pixi install
        ;;
    pip)
        rule "installing python dependencies"
        # Deliberately several fallbacks. A CUDA/HPC base image is not a Python
        # image: python3 is usually there, but python3-venv and pip frequently
        # are not, and `python3 -m venv` then fails with "ensurepip is not
        # available" - which is a confusing way to lose a rented hour.
        if ! command -v python3 >/dev/null 2>&1; then
            say "  python3 missing, installing"
            $SUDO apt-get update -qq && $SUDO apt-get install -y -qq python3
        fi
        say "  python3:  $(python3 --version 2>&1)"

        deps_ok() { "$1" -c 'import numpy, PIL, tqdm, fastapi' >/dev/null 2>&1; }

        PY=""
        if [[ -x .venv/bin/python ]] && deps_ok .venv/bin/python; then
            PY=.venv/bin/python
            say "  reusing existing .venv"
        fi

        if [[ -z "$PY" ]] && python3 -m venv .venv >/dev/null 2>&1; then
            PY=.venv/bin/python
        fi
        if [[ -z "$PY" ]]; then
            say "  python3 -m venv unavailable, installing python3-venv/pip"
            $SUDO apt-get update -qq || true
            $SUDO apt-get install -y -qq python3-venv python3-pip || true
            if python3 -m venv .venv >/dev/null 2>&1; then PY=.venv/bin/python; fi
        fi

        if [[ -n "$PY" ]]; then
            "$PY" -m pip install --quiet --upgrade pip || true
            "$PY" -m pip install --quiet "${PY_PACKAGES[@]}"
        else
            # No venv even after apt. Fall back to the system interpreter.
            # --break-system-packages is needed on PEP 668 distributions and is
            # simply unknown on older pip, hence the retry.
            say "  no venv available, installing into the system python"
            PY=python3
            python3 -m pip install --quiet --break-system-packages "${PY_PACKAGES[@]}" \
                || python3 -m pip install --quiet "${PY_PACKAGES[@]}"
        fi

        deps_ok "$PY" || { say "  dependency check still failing after install"; "$PY" -c 'import numpy, PIL, tqdm, fastapi'; exit 1; }
        say "  ok: $("$PY" --version 2>&1) at $PY"
        ;;
    *) say "unknown --deps '$DEPS' (use pip, pixi or none)"; exit 2 ;;
esac
say ""

# --------------------------------------------------------------------------
# Build
# --------------------------------------------------------------------------
rule "building"
# No unconditional `make clean`: the library depends on $(HEADERS), so make
# already rebuilds it when a source changes, and cleaning every time would make
# each resume after a preemption pay for a full rebuild it does not need. Use
# --rebuild when switching compilers, where a stale object really can linger.
if [[ $REBUILD -eq 1 ]]; then
    make clean >/dev/null 2>&1 || true
fi
make libblackhole.so
say ""

# --------------------------------------------------------------------------
# Render, syncing after each movie
# --------------------------------------------------------------------------
mkdir -p "$OUTDIR"

run_sync() {
    [[ -z "$SYNC_CMD" ]] && return 0
    local cmd="${SYNC_CMD//\{\}/$OUTDIR}"
    rule "sync: $cmd"
    # A failing sync must not abandon a run that is still producing results;
    # report it and carry on, so the remaining movies still get rendered.
    bash -c "$cmd" || say "  WARNING: sync command failed (exit $?), continuing"
}

# One render_gifs.sh invocation per movie rather than one for all of them, so
# the sync can run between them. Each invocation skips whatever is already
# present, which is also what makes a restart after preemption cheap.
if [[ "$MOVIES" == "all" ]]; then
    selected=(inclination dolly ring fp64)
else
    IFS=',' read -r -a selected <<< "$MOVIES"
fi

started=$SECONDS
for movie in "${selected[@]}"; do
    rule "rendering $movie"
    ./render_gifs.sh --movies "$movie" --errormax "$ERRORMAX" --width "$WIDTH" \
        --frames "$FRAMES" --out "$OUTDIR" --skip-build
    run_sync
done

say ""
rule "done in $((SECONDS - started))s"
ls -lh "$OUTDIR"/*.gif 2>/dev/null || say "  no GIFs produced"

if [[ $SHUTDOWN -eq 1 ]]; then
    say ""
    say "shutting down in 30s to stop billing (ctrl-c to cancel) ..."
    sleep 30
    # vast.ai's own CLI if the instance was given a key, otherwise ask the
    # kernel; on a container the latter ends the instance too.
    if command -v vastai >/dev/null 2>&1 && [[ -n "${VAST_CONTAINERLABEL:-}" ]]; then
        vastai stop instance "${VAST_CONTAINERLABEL#C.}" || $SUDO poweroff
    else
        $SUDO poweroff
    fi
fi
