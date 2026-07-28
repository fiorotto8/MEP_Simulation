#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
RUN_TAG=$(date +%Y%m%d_%H%M%S)
CURRENT_USER=$(id -un)

CPU_COUNT=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)
DEFAULT_WORKERS=$((CPU_COUNT / 3))
if (( DEFAULT_WORKERS < 1 )); then
  DEFAULT_WORKERS=1
fi

OUTPUT_ARGUMENT=${1:-"run_outputs_geometry_comparison_$RUN_TAG"}
WORKERS_PER_CONFIGURATION=${2:-$DEFAULT_WORKERS}
BASE_SEED=${GPD3D_BASE_SEED:-20260728}
BUILD_JOBS=${GPD3D_BUILD_JOBS:-$CPU_COUNT}
DOCKER_IMAGE=${GPD3D_DOCKER_IMAGE:-u1804-root62004-g4_1051}
DOCKER_BUILD_DIR_NAME=${GPD3D_DOCKER_BUILD_DIR:-build_docker_geometry}

if [[ "$OUTPUT_ARGUMENT" = /* ]]; then
  OUTPUT_ROOT=$OUTPUT_ARGUMENT
else
  OUTPUT_ROOT="$PROJECT_ROOT/$OUTPUT_ARGUMENT"
fi

if ! [[ "$WORKERS_PER_CONFIGURATION" =~ ^[1-9][0-9]*$ ]]; then
  echo "Workers per configuration must be a positive integer." >&2
  exit 2
fi

if ! [[ "$BASE_SEED" =~ ^[1-9][0-9]*$ ]]; then
  echo "GPD3D_BASE_SEED must be a positive integer." >&2
  exit 2
fi

if ! [[ "$BUILD_JOBS" =~ ^[1-9][0-9]*$ ]]; then
  echo "GPD3D_BUILD_JOBS must be a positive integer." >&2
  exit 2
fi

if ! [[ "$DOCKER_BUILD_DIR_NAME" =~ ^[A-Za-z0-9._-]+$ ]]; then
  echo "GPD3D_DOCKER_BUILD_DIR must be a simple directory name." >&2
  exit 2
fi

if ! command -v screen >/dev/null 2>&1; then
  echo "GNU screen is required on the host." >&2
  exit 1
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "Docker is required on the host." >&2
  exit 1
fi

# Prefer ordinary Docker access. A passwordless/cached sudo path is accepted,
# but never prompt from this launcher because this account may not be a sudoer.
DOCKER=(docker)
if ! docker info >/dev/null 2>&1; then
  if command -v sudo >/dev/null 2>&1 &&
     sudo -n docker info >/dev/null 2>&1; then
    DOCKER=(sudo -n docker)
  else
    cat >&2 <<EOF
The current user '$CURRENT_USER' cannot access the Docker daemon.

From an administrator account, run:

  sudo usermod -aG docker $CURRENT_USER

Then log out and back in (recommended), or start a refreshed shell with:

  newgrp docker

After that, rerun:

  ./scripts/launch_geometry_comparison.sh

Rebuilding the image as this user would not fix access to
/var/run/docker.sock. A rootless rebuild is possible but would create a
separate image store and unnecessarily rebuild ROOT and Geant4.
EOF
    exit 1
  fi
fi

if ! "${DOCKER[@]}" image inspect "$DOCKER_IMAGE" >/dev/null 2>&1; then
  echo "Docker image '$DOCKER_IMAGE' was not found." >&2
  echo "Build it first with:" >&2
  echo "  ${DOCKER[*]} build -t $DOCKER_IMAGE $PROJECT_ROOT" >&2
  exit 1
fi

if [[ -d "$OUTPUT_ROOT" ]] &&
   [[ -n "$(find "$OUTPUT_ROOT" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "Refusing to reuse non-empty output directory: $OUTPUT_ROOT" >&2
  exit 1
fi

echo "Building the shared executable inside Docker image $DOCKER_IMAGE..."
"${DOCKER[@]}" run --rm \
  -e GPD3D_BUILD_JOBS="$BUILD_JOBS" \
  -e GPD3D_BUILD_DIRECTORY="/work/$DOCKER_BUILD_DIR_NAME" \
  -v "$PROJECT_ROOT:/work" \
  -w /work \
  "$DOCKER_IMAGE" \
  bash -lc '
    mkdir -p "$GPD3D_BUILD_DIRECTORY"
    cd "$GPD3D_BUILD_DIRECTORY"
    cmake /work
    cmake --build . -- -j"$GPD3D_BUILD_JOBS"
  '

EXECUTABLE_HOST="$PROJECT_ROOT/$DOCKER_BUILD_DIR_NAME/gpd3d"
EXECUTABLE_CONTAINER="/work/$DOCKER_BUILD_DIR_NAME/gpd3d"
if [[ ! -x "$EXECUTABLE_HOST" ]]; then
  echo "Docker build did not produce executable: $EXECUTABLE_HOST" >&2
  exit 1
fi

# Freeze the inputs once so all configurations use byte-identical spectra and
# manifest contents even if the working copy changes while the runs execute.
INPUT_SPECTRA="$OUTPUT_ROOT/input_spectra"
mkdir -p "$INPUT_SPECTRA"
cp -a "$PROJECT_ROOT/spectra/." "$INPUT_SPECTRA/"
MANIFEST="$INPUT_SPECTRA/manifest.csv"
sha256sum "$MANIFEST" > "$OUTPUT_ROOT/manifest.sha256"

OUTPUT_PLUG_TUBE="$OUTPUT_ROOT/plug_and_tube"
OUTPUT_PLUG_ONLY="$OUTPUT_ROOT/plug_only"
OUTPUT_OPEN="$OUTPUT_ROOT/no_plug_no_tube"
mkdir -p "$OUTPUT_PLUG_TUBE" "$OUTPUT_PLUG_ONLY" "$OUTPUT_OPEN"

CONTAINER_PLUG_TUBE="gpd3d-sim-both-$RUN_TAG"
CONTAINER_PLUG_ONLY="gpd3d-sim-plug-$RUN_TAG"
CONTAINER_OPEN="gpd3d-sim-open-$RUN_TAG"

SESSION_PLUG_TUBE="gpd3d_both_$RUN_TAG"
SESSION_PLUG_ONLY="gpd3d_plug_$RUN_TAG"
SESSION_OPEN="gpd3d_open_$RUN_TAG"

HOST_UID=$(id -u)
HOST_GID=$(id -g)

start_container()
{
  local container_name=$1
  local configuration=$2

  "${DOCKER[@]}" run -d --rm \
    --name "$container_name" \
    --user "$HOST_UID:$HOST_GID" \
    -e HOME=/tmp/gpd3d-home \
    -v "$PROJECT_ROOT:/work" \
    -v "$OUTPUT_ROOT:/outputs" \
    -w /work \
    "$DOCKER_IMAGE" \
    bash -c '
      source /etc/profile.d/root.sh
      source /etc/profile.d/geant4.sh
      mkdir -p "$HOME"
      exec "$@"
    ' \
    bash \
    /work/scripts/run_geometry_configuration.sh \
    "$configuration" \
    "/outputs/$configuration" \
    "$WORKERS_PER_CONFIGURATION" \
    "$EXECUTABLE_CONTAINER" \
    /outputs/input_spectra \
    /outputs/input_spectra/manifest.csv \
    "$BASE_SEED" \
    >/dev/null
}

start_container "$CONTAINER_PLUG_TUBE" plug_and_tube
start_container "$CONTAINER_PLUG_ONLY" plug_only
start_container "$CONTAINER_OPEN" no_plug_no_tube

screen -dmS "$SESSION_PLUG_TUBE" \
  "$SCRIPT_DIR/monitor_geometry_configuration.sh" \
  plug_and_tube "$OUTPUT_PLUG_TUBE" "$CONTAINER_PLUG_TUBE"

screen -dmS "$SESSION_PLUG_ONLY" \
  "$SCRIPT_DIR/monitor_geometry_configuration.sh" \
  plug_only "$OUTPUT_PLUG_ONLY" "$CONTAINER_PLUG_ONLY"

screen -dmS "$SESSION_OPEN" \
  "$SCRIPT_DIR/monitor_geometry_configuration.sh" \
  no_plug_no_tube "$OUTPUT_OPEN" "$CONTAINER_OPEN"

{
  echo "$SESSION_PLUG_TUBE $CONTAINER_PLUG_TUBE plug_and_tube $OUTPUT_PLUG_TUBE"
  echo "$SESSION_PLUG_ONLY $CONTAINER_PLUG_ONLY plug_only $OUTPUT_PLUG_ONLY"
  echo "$SESSION_OPEN $CONTAINER_OPEN no_plug_no_tube $OUTPUT_OPEN"
} > "$OUTPUT_ROOT/screen_sessions.txt"

echo
echo "Started three simulations in parallel Docker containers:"
echo "  $CONTAINER_PLUG_TUBE  (plug + tube)"
echo "  $CONTAINER_PLUG_ONLY  (plug only)"
echo "  $CONTAINER_OPEN  (neither)"
echo
echo "Screen monitors:"
echo "  $SESSION_PLUG_TUBE"
echo "  $SESSION_PLUG_ONLY"
echo "  $SESSION_OPEN"
echo
echo "Output root: $OUTPUT_ROOT"
echo "Workers per configuration: $WORKERS_PER_CONFIGURATION"
echo "Frozen manifest: $MANIFEST"
echo
echo "List screens: screen -ls"
echo "Attach:       screen -r $SESSION_PLUG_TUBE"
echo "Follow log:   tail -f $OUTPUT_PLUG_TUBE/screen.log"
