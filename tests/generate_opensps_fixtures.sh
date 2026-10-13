#!/bin/bash
# Generates the libopensps parity fixture matrix (OPSF containers) into OUTDIR.
# Usage: tests/generate_opensps_fixtures.sh OUTDIR
#
# Each case line is: name|isoc|spec|dust|flags|key=value overrides
# flags: "csp" omits the SSP spectral grid (it duplicates the case's SSP-family base).

set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$( dirname "$SCRIPT_DIR" )"
OUTDIR="$(realpath -m "${1:?usage: $0 OUTDIR}")"
GENERATOR="$PROJECT_ROOT/build/generate_opensps_fixtures"

export SPS_HOME="$PROJECT_ROOT"
export FSPS_DATA_HOME="$PROJECT_ROOT"

meson compile -C "$PROJECT_ROOT/build" generate_opensps_fixtures > /dev/null
mkdir -p "$OUTDIR"

# Common CSP base: Chabrier, solar MIST, tau model observed at 10 Gyr, no dust
CSP="imf_type=1 sfh=1 tau=1.0 tage=10.0"
DUST="$CSP dust1=0.5 dust2=0.3"

cases=(
  # --- SSP family: quantities that change the SSP basis (full ssp/spec stored) ---
  "ssp_salpeter|mist|miles|DL07||imf_type=0"
  "ssp_chabrier|mist|miles|DL07||imf_type=1"
  "ssp_kroupa|mist|miles|DL07||imf_type=2"
  "ssp_vandokkum|mist|miles|DL07||imf_type=3 vdmc=0.3"
  "ssp_dave|mist|miles|DL07||imf_type=4 mdave=0.7"
  "ssp_kroupa_zlow|mist|miles|DL07||imf_type=2 zmet=4"
  "ssp_mods|mist|miles|DL07||imf_type=1 fbhb=0.2 sbss=2.0 delt=0.05 dell=0.1 pagb=0.5 redgb=1.2 agb=1.5"
  "ssp_plain|mist|miles|DL07||imf_type=1 add_agb_dust_model=0 use_wr_spectra=0 add_stellar_remnants=0"
  "ssp_agbdust2|mist|miles|DL07||imf_type=1 agb_dust=2.0"
  # --- CSP family: SFH, dust, nebular, post-processing (SSP basis = ssp_chabrier) ---
  "csp_tau|mist|miles|DL07|csp|$CSP"
  "csp_tau_allages|mist|miles|DL07|csp|imf_type=1 sfh=1 tau=2.0 tage=0"
  "csp_delayed|mist|miles|DL07|csp|imf_type=1 sfh=4 tau=3.0 tage=8.0"
  "csp_delayed_trunc_burst|mist|miles|DL07|csp|imf_type=1 sfh=4 tau=3.0 tage=12.0 sf_trunc=10.0 fburst=0.2 tburst=6.0 const=0.1"
  "csp_sfh5|mist|miles|DL07|csp|imf_type=1 sfh=5 tau=2.0 tage=12.0 sf_trunc=8.0 sf_slope=-0.5"
  "csp_sfstart|mist|miles|DL07|csp|imf_type=1 sfh=1 tau=1.0 sf_start=2.0 tage=10.0"
  "csp_ssp_offgrid|mist|miles|DL07|csp|imf_type=1 sfh=0 tage=3.3"
  "csp_dust0|mist|miles|DL07|csp|$DUST dust_type=0"
  "csp_dust0_full|mist|miles|DL07|csp|$DUST dust_type=0 dust3=0.2 frac_nodust=0.1 frac_obrun=0.1 dust1_index=-0.7 dust_index=-0.5 dust_tesc=7.5"
  "csp_dust1|mist|miles|DL07|csp|$DUST dust_type=1 uvb=0.5 mwr=3.5"
  "csp_dust2|mist|miles|DL07|csp|$DUST dust_type=2"
  "csp_dust3|mist|miles|DL07|csp|$DUST dust_type=3"
  "csp_dust4|mist|miles|DL07|csp|$DUST dust_type=4 dust_index=-0.3"
  "csp_dust5|mist|miles|DL07|csp|$DUST dust_type=5"
  "csp_dust6|mist|miles|DL07|csp|$DUST dust_type=6"
  "csp_duste_off|mist|miles|DL07|csp|$DUST dust_type=2 add_dust_emission=0"
  "csp_duste_qpah|mist|miles|DL07|csp|$DUST dust_type=2 duste_qpah=1.0 duste_umin=5.0 duste_gamma=0.1"
  "csp_nebular|mist|miles|DL07|csp|imf_type=1 sfh=1 tau=0.5 tage=0.5 add_neb_emission=1 gas_logu=-3.0 gas_logz=-0.5 dust1=0.5 dust2=0.2"
  "csp_igm|mist|miles|DL07|csp|$CSP zred=3.0 add_igm_absorption=1"
  "csp_agn|mist|miles|DL07|csp|$CSP add_agn_dust=1 fagn=0.1 agn_tau=20.0"
  "csp_smooth|mist|miles|DL07|csp|$CSP sigma_smooth=150.0 smooth_velocity=1"
  # --- Library configurations ---
  "lib_padova_basel|pdva|basel|DL07||imf_type=1 zmet=20 sfh=1 tau=1.0 tage=10.0 dust2=0.3"
  "lib_mist_c3k|mist|c3k_afe+0.0|DL07||imf_type=1 sfh=1 tau=1.0 tage=10.0 dust2=0.3"
  "lib_themis|mist|miles|themis|csp|$DUST dust_type=2"
  "lib_bpass|bpss|bpass|DL07||imf_type=1 sfh=1 tau=1.0 tage=10.0 dust2=0.3"
)

for line in "${cases[@]}"; do
  IFS='|' read -r name isoc spec dust flags overrides <<< "$line"
  extra=()
  [[ "$flags" == *csp* ]] && extra+=(--csp-only)
  echo "== $name"
  # shellcheck disable=SC2086
  (cd "$OUTDIR" && "$GENERATOR" --isoc "$isoc" --spec "$spec" --dust "$dust" "${extra[@]}" \
      --out "$name.opsf" $overrides > "$name.log" 2>&1) || { echo "FAILED: $name (see $OUTDIR/$name.log)"; continue; }
  printf '%s|%s|%s|%s|%s|%s\n' "$name" "$isoc" "$spec" "$dust" "$flags" "$overrides" >> "$OUTDIR/MANIFEST.tmp"
done

mv "$OUTDIR/MANIFEST.tmp" "$OUTDIR/MANIFEST.psv"
echo "Fixtures written to $OUTDIR"
