#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 7 ]]; then
  echo "Usage: $0 CONFIG OUTPUT_DIR WORKERS EXECUTABLE SPECTRA_DIR MANIFEST BASE_SEED" >&2
  exit 2
fi

CONFIGURATION=$1
OUTPUT_DIR=$2
WORKERS=$3
EXECUTABLE=$4
SPECTRA_DIR=$5
MANIFEST=$6
BASE_SEED=$7

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

case "$CONFIGURATION" in
  plug_and_tube)
    export GPD3D_ENABLE_W_PLUG=1
    export GPD3D_ENABLE_W_TUBE=1
    ;;
  plug_only)
    export GPD3D_ENABLE_W_PLUG=1
    export GPD3D_ENABLE_W_TUBE=0
    ;;
  no_plug_no_tube)
    export GPD3D_ENABLE_W_PLUG=0
    export GPD3D_ENABLE_W_TUBE=0
    ;;
  *)
    echo "Unknown geometry configuration: $CONFIGURATION" >&2
    exit 2
    ;;
esac

mkdir -p "$OUTPUT_DIR"
LOG_FILE="$OUTPUT_DIR/screen.log"

cd "$PROJECT_ROOT"

set +e
{
  echo "Configuration: $CONFIGURATION"
  echo "Tungsten plug: $GPD3D_ENABLE_W_PLUG"
  echo "Tungsten tube: $GPD3D_ENABLE_W_TUBE"
  echo "Manifest: $MANIFEST"
  echo "Workers: $WORKERS"
  echo "Base seed: $BASE_SEED"
  echo "Started (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)"

  python3 -u "$PROJECT_ROOT/launch_run.py" \
    --outdir "$OUTPUT_DIR" \
    --exe "$EXECUTABLE" \
    --spectra-dir "$SPECTRA_DIR" \
    --manifest "$MANIFEST" \
    --template-mac "$PROJECT_ROOT/macros/template_bg_sphere_spectrum.mac" \
    --max-workers "$WORKERS" \
    --base-seed "$BASE_SEED"

} 2>&1 | tee "$LOG_FILE"
STATUS=${PIPESTATUS[0]}
set -e
if (( STATUS == 0 )); then
  echo "Finished successfully (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)" |
    tee -a "$LOG_FILE"
else
  echo "Failed with exit status $STATUS (UTC): $(date -u +%Y-%m-%dT%H:%M:%SZ)" |
    tee -a "$LOG_FILE"
fi

echo "$STATUS" > "$OUTPUT_DIR/exit_status"
exit "$STATUS"
