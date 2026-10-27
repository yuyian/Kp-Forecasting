#!/bin/sh
# Run the realtime Kp forecast and store the pages where other people can
# open them. The local copy is always publish/. A second copy goes to the
# University of Michigan Dropbox folder when that disk is mounted, so
# Dropbox syncs it. Override that second folder with KP_SHARE_DIR.
# KP_PUBLISH_DIR still receives an extra copy when it is set.
set -eu
ROOT=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
export PATH="/usr/local/bin:/usr/bin:/bin:${PATH:-}"
PUBLISH="$ROOT/publish"
OUT="$ROOT/Rcode/Main/output/Main_VARX_Realtime"
SHARE_DEFAULT="/Users/yuyian/Library/CloudStorage/Dropbox-UniversityofMichigan/Yu Yian/Kp-Forecasting/publish"
SHARE="${KP_SHARE_DIR:-$SHARE_DEFAULT}"
LOG="$PUBLISH/cron.log"
LOCK="$PUBLISH/.run.lock"
mkdir -p "$PUBLISH"
exec >> "$LOG" 2>&1
echo "==== $(date -u '+%Y-%m-%d %H:%M:%S UTC') ===="
if [ -d "$LOCK" ]; then
  now=$(date +%s)
  if [ "$(uname)" = "Darwin" ]; then
    locked=$(stat -f %m "$LOCK")
  else
    locked=$(stat -c %Y "$LOCK")
  fi
  if [ $((now - locked)) -gt 10800 ]; then
    rmdir "$LOCK"
    echo "Removed a lock older than 3 hours."
  fi
fi
if ! mkdir "$LOCK" 2>/dev/null; then
  echo "A previous run is still going. This hour was skipped."
  exit 0
fi
trap 'rmdir "$LOCK"' EXIT
cd "$ROOT"
# GitHub's own timer has not been creating runs. Start that run from here
# when this script is the Mac job. Skip it inside GitHub Actions, or each
# run would start another one.
if [ -z "${GITHUB_ACTIONS:-}" ]; then
  GH=""
  for c in "$HOME/.local/bin/gh" /tmp/gh-cli/gh_2.103.0_macOS_arm64/bin/gh; do
    if [ -x "$c" ]; then
      GH=$c
      break
    fi
  done
  if [ -n "$GH" ] && "$GH" workflow run forecast.yml --repo yuyian/Kp-Forecasting; then
    echo "Started the GitHub forecast."
  else
    echo "Could not start the GitHub forecast."
  fi
fi
Rscript Rcode/Main/Main_VARX_Realtime.R
FILES="
VARX_Kp_realtime.html
VARX_Kp_baseline.html
VARX_Kp_realtime_latest.csv
VARX_Kp_realtime_27d.csv
issued_forecast_3h.csv
issued_forecast_daily.csv
issued_forecast_27d.csv
VARX_Kp_baseline_metrics.csv
VARX_Kp_baseline_endpoint.csv
VARX_Kp_baseline_27d_metrics.csv
VARX_Kp_baseline_27d.csv
"
publish_to() {
  dest="$1"
  mkdir -p "$dest"
  for name in $FILES; do
    cp "$OUT/$name" "$dest/$name"
  done
  echo "Published to $dest"
}
publish_to "$PUBLISH"
if [ "$SHARE" != "$PUBLISH" ]; then
  if [ -d "$(dirname "$SHARE")" ]; then
    publish_to "$SHARE" || echo "Shared folder copy failed: $SHARE"
  else
    echo "Shared folder is not mounted, skipped: $SHARE"
  fi
fi
if [ -n "${KP_PUBLISH_DIR:-}" ] && [ "$KP_PUBLISH_DIR" != "$PUBLISH" ] && [ "$KP_PUBLISH_DIR" != "$SHARE" ]; then
  publish_to "$KP_PUBLISH_DIR"
fi
