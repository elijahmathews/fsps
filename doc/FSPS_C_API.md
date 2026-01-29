# FSPS C Driver API

This document describes the C driver API exposed by FSPS for use in other language bindings.

## Build

From the repo root:

- Build the shared library:
  - `make shared`
- Build the C test driver:
  - `make test_c`

The shared library is created at build/libfsps.so. The C header is include/fsps.h.

## Initialization

- `fsps_initialize(int zin)`
- `fsps_initialize_full(int zin, int compute_vega_mags, int vactoair_flag, const char *isoc, const char *spec, const char *dust)`

Notes:
- `zin = -1` loads the full metallicity grid (required for full-Z interpolation and full-Z SSP dumps).
- `isoc`, `spec`, `dust` are optional (empty strings use defaults).

## Parameter setting

- `fsps_set_int(key, value)`
- `fsps_set_float(key, value)`
- `fsps_set_str(key, value)`
- `fsps_set_mag_compute(n_bands, mask)`
- `fsps_set_ssp_gen_age(n_age, mask)`

## Driver version and locking

- `fsps_get_driver_version(major, minor, patch)`
- `fsps_lock(status)` / `fsps_unlock(status)`

Note: The lock is a cooperative in-process guard only. It is not a real mutex and does not provide thread or interprocess safety.

Parameter names match FSPS variable names. Use `fsps_validate_params()` to catch common invalid settings.

## Computation

- `fsps_compute(spec)` — SSP/CSP depending on `sfh`.
- `fsps_compute_csp(zcontinuous)` — convenience wrapper for z-grid CSP computation (0–3).
- `fsps_compute_ssp(zin)` / `fsps_compute_ssps()`
- `fsps_compute_zdep(ztype)` — lower-level z-interpolation call.

## Output

- `fsps_get_spec(spec)` / `fsps_get_spec_peraa(spec)`
- `fsps_get_mags(zred, mags)` / `fsps_get_mags_mask(zred, mags, mask)`
- `fsps_get_stats(age, mass, lbol, sfr, mdust, mformed, emlines)`
- `fsps_get_indices(spec, indices)`

## Metadata

- `fsps_get_dims(nspec, ntfull)` and direct getters `fsps_get_nspec`, `fsps_get_ntfull`, `fsps_get_nbands`, `fsps_get_nz`, `fsps_get_nemline`
- Size helpers: `fsps_get_nt`, `fsps_get_nm`, `fsps_get_ntabmax`, `fsps_get_nindx`
- `fsps_get_lambda`, `fsps_get_emlambda`, `fsps_get_res`, `fsps_get_filter_data`
- `fsps_get_zlegend`, `fsps_get_zsol`, `fsps_get_timefull`
- `fsps_get_isochrone_dimensions`, `fsps_get_nmass_isochrone`
- `fsps_get_libraries`

## Error handling

The driver records the last error or warning internally:

- `fsps_get_last_error(&status, msg, msg_len)`
- `fsps_clear_error()`
- `fsps_set_debug(1)` to print error/warning messages to stdout.

`status == 0` means no error recorded. Non-zero values correspond to driver-specific warning/error codes.

## Array layout

All arrays are passed as flat C buffers. The driver uses Fortran order internally, so treat 2D buffers as column-major (first index varies fastest):

- `spec` is `[nspec, ntfull]`
- `mags` is `[nbands, ntfull]`
- `emlines` is `[nemline, ntfull]`
- `ssp_spec` is `[nspec, ntfull, nz]`

## Thread safety

FSPS uses global state and is not thread-safe. Do not call these APIs concurrently from multiple threads.

## Pkg-config

A pkg-config file is provided at fsps.pc. You may need to adjust `prefix` and install locations to match your environment.
