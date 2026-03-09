#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <SPS_HOME> <TEST_RUNNER_PATH>" >&2
  exit 2
fi

SPS_HOME="$1"
TEST_RUNNER="$2"

if [[ ! -d "$SPS_HOME" ]]; then
  echo "Error: SPS_HOME is not a directory: $SPS_HOME" >&2
  exit 2
fi

if [[ ! -x "$TEST_RUNNER" ]]; then
  echo "Error: test runner is not executable: $TEST_RUNNER" >&2
  exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CONVERTER="$REPO_ROOT/tools/convert_spectra_to_fsds.py"

if [[ ! -f "$CONVERTER" ]]; then
  echo "Error: converter script not found: $CONVERTER" >&2
  exit 2
fi

TMP_DIR="$(mktemp -d)"
cleanup() {
  rm -rf "$TMP_DIR"
}
trap cleanup EXIT

H5_FILE="$TMP_DIR/fsps_data_v1.h5"

python3 "$CONVERTER" \
  --sps-home "$SPS_HOME" \
  --spec-lib miles \
  --isoc-lib mist \
  --dust-type DL07 \
  --out "$H5_FILE"

export FSPS_DATA_BACKEND=hdf5
export FSPS_HDF5_DATA_PATH="$H5_FILE"

REF_FILE="$SPS_HOME/tests/data/sps_ref_MILES-1_MIST-1.bin"

if [[ ! -f "$REF_FILE" ]]; then
  echo "Error: reference file not found: $REF_FILE" >&2
  exit 2
fi

set +e
"$TEST_RUNNER" "$REF_FILE" --isoc mist --spec miles
RC=$?
set -e

GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

if [[ $RC -eq 0 ]]; then
  echo -e "${GREEN}SUCCESS: HDF5 backend dual-path validation passed.${NC}"
  exit 0
fi

echo -e "${RED}ERROR: HDF5 backend dual-path validation failed (exit code $RC).${NC}" >&2
exit 1
