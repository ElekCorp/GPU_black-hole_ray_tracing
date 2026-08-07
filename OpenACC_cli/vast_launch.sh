#!/usr/bin/env bash
#
# Rent a vast.ai GPU, render the GIFs on it from a fresh GitHub clone, copy the
# results back, and destroy the instance. Runs on YOUR machine and drives the
# vastai CLI; vast_render.sh is what actually runs on the instance.
#
# THIS SPENDS MONEY. It rents a GPU by the hour. Nothing is created without a
# confirmation prompt showing the offer and its price (--yes skips it), and
# --dry-run walks the whole flow without renting anything.
#
# THE INSTANCE IS ALWAYS DESTROYED
#
# A forgotten instance bills until you notice, which is the worst failure this
# script can have - worse than losing a render. So the destroy runs from an EXIT
# trap and fires on success, on error, and on ctrl-c alike. The instance id is
# also written to .vast_instance_id the moment it exists, and printed on the way
# out, so if this script is killed hard enough to skip its own trap you still
# have the one command needed to clean up.
#
#   --keep           never destroy (you clean up by hand)
#   --keep-on-error  destroy on success, keep it for debugging on failure
#
# ON-DEMAND vs INTERRUPTIBLE
#
# On-demand by default, deliberately. This job is short - tens of minutes - and
# a large part of that is pulling a multi-GB image and compiling with nvc++,
# which an interruption makes you pay for twice. At these durations the spot
# discount is worth cents, while the handling it needs is real. Use --bid PRICE
# if you want interruptible anyway; the work resumes (vast_render.sh skips GIFs
# that already exist), but you have to restart the stopped instance yourself.
# Interruptible earns its keep on the long precompute job in
# docs/precompute-investigation.md, not on this one.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

GPU_NAME="A100"
MAX_PRICE="1.50"
DISK=64
IMAGE="nvcr.io/nvidia/nvhpc:24.5-devel-cuda12.4-ubuntu22.04"
REPO_URL="https://github.com/ElekCorp/GPU_black-hole_ray_tracing.git"
GIT_REF=""
BID=""
OUTDIR="vast_results"
RENDER_ARGS="--width 2048 --frames 128"
SSH_KEY="$HOME/.ssh/id_ed25519"
KEEP=0
KEEP_ON_ERROR=0
ASSUME_YES=0
DRY_RUN=0
SMOKE=0
ALLOW_UNPUSHED=0
POLL_SECONDS=30
MAX_WAIT_MINUTES=180
SELF_TEST=0
CLEANUP=0
ID_FILE=".vast_instance_id"

usage() {
    cat <<'USAGE'
usage: vast_launch.sh [options]

  --gpu NAME          GPU model to search for (default A100)
  --max-price USD     highest $/hr to accept (default 1.50)
  --disk GB           instance disk (default 64; the nvhpc image is large)
  --image REF         docker image (default nvcr.io/nvidia/nvhpc:24.5-devel-cuda12.4-ubuntu22.04)
  --repo URL          git URL to clone (default: this repo's origin, as https)
  --ref NAME          branch or tag to clone (default: current branch)
  --bid USD           rent INTERRUPTIBLE at this bid instead of on-demand
  --out DIR           where to copy results (default vast_results)
  --render-args 'A'   passed to vast_render.sh (default '--width 1280 --frames 96')
  --ssh-key PATH      private key to use (default ~/.ssh/id_ed25519)
  --keep              do not destroy the instance at the end
  --keep-on-error     destroy on success only, keep it if something failed
  --allow-unpushed    rent even though local commits are not on GitHub
                      (the instance would run without them)
  --smoke             one tiny movie instead of the full set, to prove the
                      instance can build and render before paying for the rest
  --poll SECONDS      how often to check progress (default 30)
  --max-wait MIN      give up and fetch logs after this long (default 180)
  -y, --yes           do not prompt before renting
  --dry-run           do everything except rent, run and destroy
  --cleanup           list every instance on the account and destroy any this
                      script left running, then exit
  --self-test         check the status-parsing logic and exit (no network)
  -h, --help          this
USAGE
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --gpu)          GPU_NAME="$2"; shift 2 ;;
        --max-price)    MAX_PRICE="$2"; shift 2 ;;
        --disk)         DISK="$2"; shift 2 ;;
        --image)        IMAGE="$2"; shift 2 ;;
        --repo)         REPO_URL="$2"; shift 2 ;;
        --ref)          GIT_REF="$2"; shift 2 ;;
        --bid)          BID="$2"; shift 2 ;;
        --out)          OUTDIR="$2"; shift 2 ;;
        --render-args)  RENDER_ARGS="$2"; shift 2 ;;
        --ssh-key)      SSH_KEY="$2"; shift 2 ;;
        --keep)          KEEP=1; shift ;;
        --keep-on-error) KEEP_ON_ERROR=1; shift ;;
        --poll)         POLL_SECONDS="$2"; shift 2 ;;
        --max-wait)     MAX_WAIT_MINUTES="$2"; shift 2 ;;
        -y|--yes)       ASSUME_YES=1; shift ;;
        --dry-run)      DRY_RUN=1; shift ;;
        --smoke)        SMOKE=1; shift ;;
        --allow-unpushed) ALLOW_UNPUSHED=1; shift ;;
        --self-test)    SELF_TEST=1; shift ;;
        --cleanup)      CLEANUP=1; shift ;;
        -h|--help)      usage; exit 0 ;;
        *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

# A deliberately tiny job: it exercises the whole chain - clone, deps, nvc++
# build, one render, copy back - in a couple of minutes, which is the cheap way
# to find out whether an image works before committing to the full set.
if [[ $SMOKE -eq 1 ]]; then
    RENDER_ARGS="--movies inclination --frames 4 --width 256"
fi

say()  { printf '%s\n' "$*"; }
rule() { printf -- '\n---- %s\n' "$*"; }
die()  { printf 'error: %s\n' "$*" >&2; exit 1; }

# Does this STATUS value mean the job has stopped?
#
# /workspace/STATUS carries BOTH progress markers (RUNNING-CLONE, RUNNING-RENDER)
# and terminal ones (DONE, FAILED-*), so "non-empty" is not the same as
# "finished". Treating any non-empty value as terminal is precisely how a run
# was torn down seconds after it started: the poll saw RUNNING-CLONE on its
# first look, copied an output directory that did not exist yet, and destroyed
# the instance while nvc++ was still compiling. Exercised by --self-test.
is_terminal_status() {
    case "$1" in
        DONE|FAILED-*|TIMEOUT) return 0 ;;
        *) return 1 ;;
    esac
}

# Cheap guard for the above, since getting it wrong costs a rented instance
# rather than a failed assertion.
if [[ $SELF_TEST -eq 1 ]]; then
    fails=0
    for value in DONE FAILED-CLONE FAILED-RENDER TIMEOUT; do
        is_terminal_status "$value" || { echo "FAIL: '$value' should be terminal"; fails=1; }
    done
    for value in "" RUNNING-CLONE RUNNING-RENDER; do
        ! is_terminal_status "$value" || { echo "FAIL: '$value' should NOT be terminal"; fails=1; }
    done
    [[ $fails -eq 0 ]] && echo "self-test: PASS" || echo "self-test: FAILED"
    exit $fails
fi

# One EXIT trap for the whole script, armed before anything it has to clean up
# exists. Separate traps per resource do not compose - each `trap ... EXIT`
# REPLACES the previous one, so the earlier handler silently stops running.
INSTANCE_ID=""
JOB_OK=0
OFFERS_FILE=""
ONSTART=""

cleanup() {
    local status=$?
    # Explicit ifs, not `[[ ... ]] && rm`: under set -e a false test at the top
    # of a function is a failing statement, and the shell would leave cleanup
    # right there - skipping the destroy below, which is the one thing this
    # function exists to guarantee.
    if [[ -n "$OFFERS_FILE" ]]; then rm -f "$OFFERS_FILE"; fi
    if [[ -n "$ONSTART" ]]; then rm -f "$ONSTART"; fi
    if [[ -z "$INSTANCE_ID" ]]; then exit $status; fi
    if [[ $KEEP -eq 1 ]]; then
        say ""
        say "--keep: instance $INSTANCE_ID left running. It is BILLING until you run:"
        say "    $VASTAI destroy instance $INSTANCE_ID"
        exit $status
    fi
    if [[ $KEEP_ON_ERROR -eq 1 && $JOB_OK -eq 0 ]]; then
        say ""
        say "job did not finish and --keep-on-error is set: instance $INSTANCE_ID left running."
        say "It is BILLING until you run:"
        say "    $VASTAI destroy instance $INSTANCE_ID"
        say "    $VASTAI ssh-url $INSTANCE_ID     # to look at /workspace/job.log"
        exit $status
    fi
    rule "destroying instance $INSTANCE_ID"
    $VASTAI destroy instance "$INSTANCE_ID" --yes || {
        say "  DESTROY FAILED - the instance is still billing. Run this yourself:"
        say "    $VASTAI destroy instance $INSTANCE_ID"
    }
    rm -f "$ID_FILE"
    exit $status
}
trap cleanup EXIT INT TERM

# --------------------------------------------------------------------------
# Locate the vastai CLI
# --------------------------------------------------------------------------
# Installed under pixi rather than on PATH in the common case, so look for the
# pixi project too instead of demanding the user activate a shell first.
resolve_vastai() {
    if [[ -n "${VASTAI:-}" ]]; then echo "$VASTAI"; return 0; fi
    if command -v vastai >/dev/null 2>&1; then echo "vastai"; return 0; fi
    local candidate
    for candidate in "$HOME/vastai" "$HOME/.vastai" ../vastai; do
        if [[ -f "$candidate/pixi.toml" ]] && command -v pixi >/dev/null 2>&1; then
            echo "pixi run --manifest-path $candidate/pixi.toml vastai"; return 0
        fi
    done
    return 1
}
VASTAI=$(resolve_vastai) || die "vastai CLI not found. Put it on PATH, set VASTAI='...', or keep its pixi project in ~/vastai"

# --------------------------------------------------------------------------
# --cleanup: find and destroy anything this script left behind
# --------------------------------------------------------------------------
# The safety net for the case the EXIT trap cannot cover - a laptop that slept,
# a killed terminal, a destroy call that failed. Lists EVERY instance on the
# account, not just ones from here, because an orphan you have forgotten about
# is exactly the one still costing money; but only offers to destroy the ones
# this script labelled, so it cannot take out unrelated work.
if [[ $CLEANUP -eq 1 ]]; then
    rule "instances on this account"
    instances=$($VASTAI show instances --raw 2>/dev/null || echo '[]')
    LIST_FILE=$(mktemp -t vast_list.XXXXXX)
    printf '%s' "$instances" > "$LIST_FILE"
    mine=$(python3 - "$LIST_FILE" <<'PY'
import json, sys
try:
    with open(sys.argv[1]) as handle:
        inst = json.load(handle) or []
except Exception:
    inst = []
if not inst:
    print("#none", file=sys.stderr)
    raise SystemExit
ours = []
for i in inst:
    label = i.get("label") or "-"
    mark = "  <- from vast_launch.sh" if label == "blackhole-render" else ""
    if label == "blackhole-render":
        ours.append(str(i.get("id")))
    print("  id=%-10s %-9s %-18s %-14s $%.3f/hr  up %.2fh%s" % (
        i.get("id"), i.get("actual_status") or "?", label,
        i.get("gpu_name") or "?", i.get("dph_total") or 0,
        (i.get("duration") or 0) / 3600.0, mark), file=sys.stderr)
print(" ".join(ours))
PY
)
    rm -f "$LIST_FILE"
    if [[ -z "$mine" ]]; then
        say ""
        say "nothing labelled blackhole-render is running."
        say "If something above is yours and unwanted: $VASTAI destroy instance ID"
        exit 0
    fi
    say ""
    say "labelled blackhole-render: $mine"
    if [[ $ASSUME_YES -eq 0 ]]; then
        read -r -p "Destroy these? [y/N] " reply
        [[ "$reply" =~ ^[Yy]$ ]] || die "cancelled - they are still billing"
    fi
    for id in $mine; do
        say "  destroying $id"
        $VASTAI destroy instance "$id" --yes || say "  FAILED to destroy $id - do it in the web console"
    done
    rm -f "$ID_FILE"
    exit 0
fi

# --------------------------------------------------------------------------
# Preflight
# --------------------------------------------------------------------------
rule "preflight"
say "  vastai:   $VASTAI"
$VASTAI --version >/dev/null 2>&1 || die "'$VASTAI --version' failed"

[[ -f "$HOME/.config/vastai/vast_api_key" ]] || die "no API key at ~/.config/vastai/vast_api_key. Run: $VASTAI set api-key YOUR_KEY"
[[ -f "$SSH_KEY" ]] || die "ssh key $SSH_KEY not found (use --ssh-key)"
[[ -f "$SSH_KEY.pub" ]] || die "public key $SSH_KEY.pub not found"
command -v python3 >/dev/null 2>&1 || die "python3 needed to parse the offer list"

# The instance clones from GitHub, so anything not pushed will simply not be
# there. Cheap to check here, confusing to discover an hour later in the output.
if [[ -z "$GIT_REF" ]]; then
    GIT_REF=$(git rev-parse --abbrev-ref HEAD)
fi
if [[ -z "$REPO_URL" ]]; then
    origin=$(git remote get-url origin 2>/dev/null) || die "no git origin; pass --repo"
    # git@github.com:owner/repo.git -> https://github.com/owner/repo.git, since
    # the instance has no deploy key for the ssh form.
    REPO_URL=$(printf '%s' "$origin" | sed -E 's#^git@([^:]+):#https://\1/#')
fi
say "  repo:     $REPO_URL"
say "  ref:      $GIT_REF"

if git rev-parse --verify --quiet "origin/$GIT_REF" >/dev/null 2>&1; then
    ahead=$(git rev-list --count "origin/$GIT_REF..HEAD" 2>/dev/null || echo 0)
    if [[ "$ahead" -gt 0 ]]; then
        # A hard stop rather than a warning. Renting a GPU to run code that is
        # still sitting on your laptop is the most expensive mistake available
        # here, and it is invisible from the output - the render succeeds, it is
        # just the wrong render.
        say ""
        say "  $ahead local commit(s) on $GIT_REF are NOT pushed:"
        git log --oneline "origin/$GIT_REF..HEAD" | sed 's/^/      /'
        say ""
        say "  The instance clones from GitHub, so it would run without them."
        if [[ $ALLOW_UNPUSHED -eq 0 ]]; then
            die "push first (git push origin $GIT_REF), or pass --allow-unpushed to rent anyway"
        fi
        say "  --allow-unpushed: continuing with the pushed version regardless"
    fi
    say "  in sync with origin/$GIT_REF"
else
    say "  WARNING: origin/$GIT_REF not found locally; cannot verify it exists on GitHub."
fi

# --------------------------------------------------------------------------
# Pick an offer
# --------------------------------------------------------------------------
rule "searching offers"
# The GPU model is filtered here rather than in the query. vast.ai's gpu_name
# values carry spaces ("A100 SXM4", "A100 PCIE", "RTX 4090") and the query
# language wants them underscored and exact, so `gpu_name=A100` silently matches
# nothing at all - measured: 0 results, against 18 for a substring match. A
# substring also lets --gpu A100 pick up both the SXM4 and PCIE variants.
#
# --limit is needed as well: the default page caps the result set (64 offers
# here), which quietly hides most of the market from the price comparison.
QUERY="num_gpus=1 rentable=true verified=true disk_space>=$DISK inet_down>=100 dph<=$MAX_PRICE"
# Via a file, not a shell variable: 500 offers is well past the argv size limit
# ("Argument list too long"), and the heredoc below already occupies stdin.
OFFERS_FILE=$(mktemp -t vast_offers.XXXXXX)
if [[ -n "$BID" ]]; then
    say "  interruptible (bid \$$BID/hr), $GPU_NAME, up to \$$MAX_PRICE/hr"
    $VASTAI search offers "$QUERY" --interruptible --limit 500 --raw > "$OFFERS_FILE" 2>/dev/null || echo '[]' > "$OFFERS_FILE"
else
    say "  on-demand, $GPU_NAME, up to \$$MAX_PRICE/hr"
    $VASTAI search offers "$QUERY" --limit 500 --raw > "$OFFERS_FILE" 2>/dev/null || echo '[]' > "$OFFERS_FILE"
fi

read -r OFFER_ID OFFER_PRICE OFFER_DESC <<<"$(
python3 - "$GPU_NAME" "$OFFERS_FILE" <<'PY'
import json, sys
want = sys.argv[1].lower().replace("_", " ")
try:
    with open(sys.argv[2]) as handle:
        offers = json.load(handle) or []
except Exception:
    offers = []
offers = [o for o in offers if want in (o.get("gpu_name") or "").lower()]
if not offers:
    print("NONE 0 no-offers"); raise SystemExit
# Cheapest that still has decent reliability - a cheap machine that drops the
# job halfway is not cheap. Falls back to the whole list rather than failing if
# nothing clears the bar.
ok = [o for o in offers if (o.get("reliability2") or 0) >= 0.95] or offers
best = min(ok, key=lambda o: o.get("dph_total") or o.get("dph_base") or 1e9)
price = best.get("dph_total") or best.get("dph_base") or 0
desc = "%s|%dGB-disk|%.0fMbps|rel=%.3f" % (
    (best.get("gpu_name") or "?").replace(" ", "-"),
    int(best.get("disk_space") or 0),
    best.get("inet_down") or 0,
    best.get("reliability2") or 0)
print(best.get("id"), "%.4f" % price, desc)
PY
)"

[[ "$OFFER_ID" == "NONE" ]] && die "no '$GPU_NAME' offers matched: $QUERY. Try --gpu, --max-price or --disk."

say "  offer:    $OFFER_ID  \$$OFFER_PRICE/hr  $OFFER_DESC"
say "  image:    $IMAGE"
say "  render:   $RENDER_ARGS"
say "  results:  $OUTDIR"
say ""
say "  Estimated: \$$OFFER_PRICE for the first hour."
say "  Pulling the image and compiling takes a good part of that before any"
say "  frame is rendered, so budget for at least one hour even if the render"
say "  itself is quick."

if [[ $ASSUME_YES -eq 0 && $DRY_RUN -eq 0 ]]; then
    say ""
    read -r -p "Rent this instance and start? [y/N] " reply
    [[ "$reply" =~ ^[Yy]$ ]] || die "cancelled - nothing was rented"
fi

# --------------------------------------------------------------------------
# The onstart script the instance runs on its own
# --------------------------------------------------------------------------
# Written to a file and passed with --onstart rather than --onstart-cmd, so
# nothing has to survive a round trip through shell quoting.
ONSTART=$(mktemp -t vast_onstart.XXXXXX)
cat > "$ONSTART" <<EOF
#!/bin/bash
# mkdir BEFORE the redirect, and this order matters. Bash does not abort on a
# failed \`exec >\` redirect - it prints "No such file or directory" and carries
# on with the redirect simply not in effect - so redirecting into /workspace
# before creating it loses the whole log silently while the script still runs to
# completion and still writes STATUS. That is exactly how a failed run came back
# with a status and no diagnostics at all.
mkdir -p /workspace 2>/dev/null || true
LOGFILE=/workspace/job.log
: > "\$LOGFILE" 2>/dev/null || LOGFILE=/root/job.log
# tee, not plain truncation, so the same output also reaches vast.ai's own
# onstart log and stays reachable through \`vastai logs\` even if the filesystem
# the launcher scps from is not the one this wrote to.
exec > >(tee -a "\$LOGFILE") 2>&1
set -x

# Stage written as it goes, so a failure names its phase even if the log is
# unreachable for any reason.
stage() { set +x; echo "\$1" > /workspace/STATUS; set -x; }

rm -f /workspace/STATUS
stage RUNNING-CLONE
cd /workspace
if [ ! -d repo ]; then
    git clone --depth 1 --single-branch --branch "$GIT_REF" "$REPO_URL" repo || { stage FAILED-CLONE; exit 1; }
fi
cd repo/OpenACC_cli

stage RUNNING-RENDER
if ./vast_render.sh --outdir /workspace/out $RENDER_ARGS; then
    stage DONE
else
    stage FAILED-RENDER
fi
# Kept next to the results so one scp of the output directory brings the log too.
mkdir -p /workspace/out && cp -f "\$LOGFILE" /workspace/out/job.log 2>/dev/null || true
EOF

if [[ $DRY_RUN -eq 1 ]]; then
    rule "dry run - would create instance $OFFER_ID with this onstart:"
    sed 's/^/    /' "$ONSTART"
    rule "dry run - would then poll /workspace/STATUS, copy /workspace/out to $OUTDIR, and destroy"
    say ""
    say "nothing was rented."
    exit 0
fi

# --------------------------------------------------------------------------
# Create, with the destroy trap armed immediately
# --------------------------------------------------------------------------
rule "creating instance"
CREATE_ARGS=(create instance "$OFFER_ID" --image "$IMAGE" --disk "$DISK"
             --onstart "$ONSTART" --ssh --direct --label blackhole-render)
[[ -n "$BID" ]] && CREATE_ARGS+=(--bid_price "$BID")

CREATE_OUT=$($VASTAI "${CREATE_ARGS[@]}" --raw)
INSTANCE_ID=$(printf '%s' "$CREATE_OUT" | python3 -c 'import json,sys; d=json.load(sys.stdin); print(d.get("new_contract") or d.get("id") or "")' 2>/dev/null || true)
[[ -n "$INSTANCE_ID" ]] || { say "$CREATE_OUT"; die "could not read the new instance id from the create response"; }

# Recorded before anything else can fail, so a hard kill still leaves a trail.
printf '%s\n' "$INSTANCE_ID" > "$ID_FILE"
say "  instance: $INSTANCE_ID   (also written to $ID_FILE)"

$VASTAI attach ssh "$INSTANCE_ID" "$(cat "$SSH_KEY.pub")" >/dev/null 2>&1 || \
    say "  note: could not attach the ssh key; relying on the account default"

# --------------------------------------------------------------------------
# Wait for it to boot
# --------------------------------------------------------------------------
rule "waiting for the instance to start"
SSH_HOST=""; SSH_PORT=""
for _ in $(seq 1 60); do
    url=$($VASTAI ssh-url "$INSTANCE_ID" 2>/dev/null || true)
    if [[ "$url" =~ ssh://([^@]+)@([^:]+):([0-9]+) ]]; then
        SSH_HOST="${BASH_REMATCH[2]}"; SSH_PORT="${BASH_REMATCH[3]}"
        if ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
               -o ConnectTimeout=10 -o BatchMode=yes -i "$SSH_KEY" \
               -p "$SSH_PORT" "root@$SSH_HOST" true >/dev/null 2>&1; then
            say "  ssh up: root@$SSH_HOST:$SSH_PORT"
            break
        fi
    fi
    sleep 10
done
[[ -n "$SSH_HOST" ]] || die "instance never became reachable over ssh"

remote() {
    ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
        -o ConnectTimeout=15 -i "$SSH_KEY" -p "$SSH_PORT" "root@$SSH_HOST" "$@"
}

# --------------------------------------------------------------------------
# Wait for the job
# --------------------------------------------------------------------------
rule "rendering (polling every ${POLL_SECONDS}s; the image pull and nvc++ build come first)"
STATUS=""
deadline=$(( SECONDS + MAX_WAIT_MINUTES * 60 ))
while :; do
    STATUS=$(remote 'cat /workspace/STATUS 2>/dev/null' 2>/dev/null || true)
    if is_terminal_status "$STATUS"; then break; fi
    if (( SECONDS > deadline )); then
        say "  giving up after $MAX_WAIT_MINUTES minutes (last stage: ${STATUS:-none})"
        STATUS="TIMEOUT"
        break
    fi
    tail=$(remote 'tail -1 /workspace/job.log 2>/dev/null' 2>/dev/null || true)
    printf '  [%s] %-14s %s\n' "$(date +%H:%M:%S)" "${STATUS:-starting}" "${tail:0:110}"
    sleep "$POLL_SECONDS"
done
say "  status: $STATUS"

# --------------------------------------------------------------------------
# Copy results back BEFORE destroying, whatever the status
# --------------------------------------------------------------------------
# Deliberately unconditional: a failed render often still produced some movies,
# and the log is the only way to find out why the rest did not.
rule "copying results to $OUTDIR"
mkdir -p "$OUTDIR"
scp -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
    -i "$SSH_KEY" -P "$SSH_PORT" -r "root@$SSH_HOST:/workspace/out/." "$OUTDIR/" 2>/dev/null || \
    say "  no /workspace/out to copy (the render did not get that far)"

# Logs from every place one could have ended up. Coming back with a failure
# status and nothing to read is the worst outcome this script can produce, and
# it has happened: a redirect into a directory that did not exist yet sent the
# whole log somewhere else while the run carried on regardless.
fetch_log() {
    scp -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
        -i "$SSH_KEY" -P "$SSH_PORT" "root@$SSH_HOST:$1" "$OUTDIR/$2" 2>/dev/null \
        && { say "  fetched $1"; return 0; }
    return 1
}
fetch_log /workspace/job.log job.log || fetch_log /root/job.log job.log || \
    say "  no job.log on the instance"
# vast.ai's own capture of the onstart script's output, which survives even when
# the script's own redirect did not take effect.
$VASTAI logs "$INSTANCE_ID" > "$OUTDIR/vast_onstart.log" 2>/dev/null && \
    say "  fetched vast.ai instance log -> $OUTDIR/vast_onstart.log" || true

ls -lh "$OUTDIR" 2>/dev/null || true

# The tail of whatever was captured, so a failure is visible right here instead
# of only in a file the user then has to go and find.
if [[ "$STATUS" != "DONE" ]]; then
    for candidate in "$OUTDIR/job.log" "$OUTDIR/vast_onstart.log"; do
        if [[ -s "$candidate" ]]; then
            rule "last 40 lines of $(basename "$candidate")"
            tail -40 "$candidate"
            break
        fi
    done
fi

if [[ "$STATUS" == "DONE" ]]; then
    JOB_OK=1
    rule "finished"
else
    say ""
    say "job reported '$STATUS' - see $OUTDIR/job.log"
fi
# cleanup() runs from the EXIT trap and destroys the instance.
