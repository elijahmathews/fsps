#!/bin/bash
# The demo's parameter vectors through legacy (upstream) FSPS and through this fork, on the same
# data, settings, SSP weights and outputs: speed and numerical differences (see README.md).
#
# Usage: SPS_HOME=/path/to/legacy/fsps demo/sed_sampling/legacy/run_legacy_compare.sh [N]
#   SPS_HOME    a legacy FSPS checkout (its src/ and data, as legacy FSPS itself uses SPS_HOME)
#   N           number of parameter vectors (default 1000)
# Environment:
#   LEGACY_REF  git revision of SPS_HOME to build (default: its HEAD if that still has the MIST
#               isochrone layout of this fork's data, isoc_z*.dat; otherwise 82a8735, the last
#               upstream revision with it)
#   GPU_LIB, CPU_LIB  fork libraries (default build_gpu32cb/ and build_cpu/libfsps.so)
#   OUT         results directory (default demo/sed_sampling/results_a100/legacy)
#   GPU_ENV, CPU_ENV, LEGACY_ENV  shell snippets run before the fork GPU, fork CPU and legacy
#               steps (e.g. GPU_ENV='module load nvhpc/25.1' LEGACY_ENV='module load gcc/13.4.0'
#               where the two compilers cannot be loaded together); default none
# Needs gfortran (legacy builds with its own Makefile flags, -O3 -cpp -fPIC), python3 with numpy
# and matplotlib, and the fork libraries built. The fork reads its own data/ (SPS_HOME is unset
# for the fork steps: it would take precedence over FSPS_DATA_HOME).

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
LEGACY="$(realpath "${SPS_HOME:?set SPS_HOME to a legacy FSPS checkout}")"
N="${1:-1000}"
OUT="$(realpath -m "${OUT:-$HERE/../results_a100/legacy}")"
GPU_LIB="$(realpath -m "${GPU_LIB:-$ROOT/build_gpu32cb/libfsps.so}")"
CPU_LIB="$(realpath -m "${CPU_LIB:-$ROOT/build_cpu/libfsps.so}")"
WORK="$OUT/work"
mkdir -p "$WORK"

# 1. Legacy source at a revision whose isochrone layout matches the fork's data
if git -C "$LEGACY" rev-parse --git-dir > /dev/null 2>&1; then
    if [ -z "${LEGACY_REF:-}" ]; then
        if git -C "$LEGACY" cat-file -e HEAD:ISOCHRONES/MIST/isoc_zp0.00.dat 2> /dev/null; then
            LEGACY_REF=HEAD
        else
            LEGACY_REF=82a8735
        fi
    fi
    REV="$(git -C "$LEGACY" rev-parse --short "$LEGACY_REF")"
    rm -rf "$WORK/src" && mkdir -p "$WORK/src_tmp"
    git -C "$LEGACY" archive "$REV" src | tar -x -C "$WORK/src_tmp"
    mv "$WORK/src_tmp/src" "$WORK/src" && rmdir "$WORK/src_tmp"
else
    REV="(not a git checkout)"
    rm -rf "$WORK/src" && cp -r "$LEGACY/src" "$WORK/src"
fi
echo "legacy FSPS: $LEGACY at $REV"

# 2. Data: a legacy-layout SPS_HOME of links to this fork's files, verified file by file
#    against the legacy revision's own data (data_provenance.txt)
rm -rf "$WORK/sps_home" && mkdir -p "$WORK/sps_home"
ln -s "$ROOT/data/isochrones" "$WORK/sps_home/ISOCHRONES"
ln -s "$ROOT/data/spectra" "$WORK/sps_home/SPECTRA"
ln -s "$ROOT/data/dust" "$WORK/sps_home/dust"
ln -s "$ROOT/data/nebular" "$WORK/sps_home/nebular"
ln -s "$ROOT/data" "$WORK/sps_home/data"
FSPS_DATA_HOME="$ROOT" python3 "$HERE/compare_legacy.py" verify-data --legacy "$LEGACY" --ref "$REV" --data "$ROOT/data" --dir "$OUT"

# 3. Legacy builds: stock (only the SSP-weight input); aligned (also the parity patches of
#    tools/legacy_parity: the fork's physical constants and upstream's AGN fix); aligned_rv (also
#    Calzetti R_V = 4.05 exactly, one of the fork's documented intentional differences)
build_legacy() {
    local name=$1 dir="$WORK/build_$1"
    rm -rf "$dir" && cp -r "$WORK/src" "$dir"
    sed -i -f "$HERE/weights_csp_gen.sed" "$dir/csp_gen.f90"
    sed -i -f "$HERE/weights_sps_vars.sed" "$dir/sps_vars.f90"
    grep -q 'use_given_weights.eq.1' "$dir/csp_gen.f90" || { echo "weights patch did not apply" >&2; return 1; }
    grep -q 'given_weights=0' "$dir/sps_vars.f90" || { echo "weights patch did not apply" >&2; return 1; }
    if [ "$name" != stock ]; then
        sed -i -f "$ROOT/tools/legacy_parity/patches/constants.sed" "$dir/sps_vars.f90"
        sed -i -f "$ROOT/tools/legacy_parity/patches/agn_dust.sed" "$dir/agn_dust.f90"
        grep -q '4.0479831348E10' "$dir/sps_vars.f90" || { echo "constants patch did not apply" >&2; return 1; }
        grep -q 'pset%dust2 \* attn_curve' "$dir/agn_dust.f90" || { echo "AGN patch did not apply" >&2; return 1; }
    fi
    if [ "$name" = aligned_rv ]; then
        sed -i -f "$HERE/calzetti_rv.sed" "$dir/attn_curve.f90"
        [ "$(grep -c ') + 1.782$' "$dir/attn_curve.f90")" = 4 ] || { echo "R_V patch did not apply" >&2; return 1; }
    fi
    cp "$HERE/legacy_batch.f90" "$dir/"
    (
        cd "$dir"
        rm -f ./*.o ./*.mod
        make -s simple > build.log 2>&1
        objs=$(ls ./*.o | grep -vE '/(simple|lesssimple|autosps|spec_bin)\.o$')
        # shellcheck disable=SC2086
        gfortran -O3 -cpp -fPIC legacy_batch.f90 $objs -o legacy_batch >> build.log 2>&1
    ) || { echo "legacy build $name failed (see $dir/build.log)" >&2; return 1; }
}
pids=()
for name in stock aligned aligned_rv; do
    (eval "${LEGACY_ENV:-}"; build_legacy $name) & pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done

# 4. Fork: the vectors are drawn (and the inputs written) by the first call
rm -f "$OUT/inputs.npz"
(eval "${GPU_ENV:-}"; env -u SPS_HOME FSPS_DATA_HOME="$ROOT" \
    python3 "$HERE/compare_legacy.py" fork --lib "$GPU_LIB" --tag gpu --dir "$OUT" --n "$N")
(eval "${CPU_ENV:-}"; env -u SPS_HOME FSPS_DATA_HOME="$ROOT" \
    python3 "$HERE/compare_legacy.py" fork --lib "$CPU_LIB" --tag cpu --dir "$OUT")

# 5. Legacy, one vector at a time (one core), both builds
for name in stock aligned aligned_rv; do
    (eval "${LEGACY_ENV:-}"; cd "$WORK/build_$name" && SPS_HOME="$WORK/sps_home" ./legacy_batch \
        "$OUT/settings.txt" "$OUT/vectors.bin" "$OUT/legacy_$name.bin") | tee "$OUT/legacy_$name.log"
done

# 6. Report
FSPS_DATA_HOME="$ROOT" python3 "$HERE/compare_legacy.py" report --dir "$OUT"
