#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 CONFIGURATION OUTPUT_DIR CONTAINER_NAME" >&2
  exit 2
fi

CONFIGURATION=$1
OUTPUT_DIR=$2
CONTAINER_NAME=$3
LOG_FILE="$OUTPUT_DIR/screen.log"
STATUS_FILE="$OUTPUT_DIR/exit_status"

mkdir -p "$OUTPUT_DIR"
touch "$LOG_FILE"

echo "Monitoring $CONFIGURATION"
echo "Docker container: $CONTAINER_NAME"
echo "Log: $LOG_FILE"
echo

tail -n +1 -F "$LOG_FILE" &
TAIL_PID=$!

cleanup()
{
  kill "$TAIL_PID" 2>/dev/null || true
  wait "$TAIL_PID" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

while [[ ! -f "$STATUS_FILE" ]]; do
  sleep 2
done

STATUS=$(<"$STATUS_FILE")
cleanup
trap - EXIT INT TERM

echo
if [[ "$STATUS" == "0" ]]; then
  echo "$CONFIGURATION completed successfully."
  exit 0
fi

echo "$CONFIGURATION failed with exit status $STATUS."
exit "$STATUS"
