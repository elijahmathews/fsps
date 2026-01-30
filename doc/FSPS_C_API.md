# FSPS C Driver API

This document describes the C driver API exposed by FSPS for use in other language bindings.

## Build

From the repo root:

- Build the shared library:
  - `make shared`
- Build the C test driver:
  - `make test_c`

The shared library is created at build/libfsps.so.* (with symlinks build/libfsps.so and build/libfsps.so.<soname>). The C header is include/fsps.h.

## Initialization (context-based)
- `fsps_context_create(int *handle, int *status)`
- `fsps_context_setup(int zin, const char *isoc, const char *spec, const char *dust, int handle, int *status)`
- `fsps_context_destroy(int handle, int *status)`

Notes:
- `zin = -1` loads the full metallicity grid (required for full-Z interpolation and full-Z SSP dumps).
- `isoc`, `spec`, `dust` are optional (empty strings use defaults).

## Parameter setting

- `fsps_context_set_int(handle, key, value, status)`
- `fsps_context_set_float(handle, key, value, status)`
- `fsps_context_set_str(handle, key, value, status)`

## Driver version and locking

- `fsps_get_driver_version(major, minor, patch)`
- `fsps_lock(status)` / `fsps_unlock(status)`

Note: The lock is a cooperative in-process guard only. It is not a real mutex and does not provide thread or interprocess safety.

Parameter names match FSPS variable names.

## Computation

Context-based SSP compute (requires prior `fsps_context_setup`):
- `fsps_context_compute_ssp(int handle, double *spec, double *mass, double *lbol, int *status)`

Note: `spec` is column-major with dimensions `[nspec, ntfull]`.

## Output

Context-based output access is currently limited to SSP arrays returned from `fsps_context_compute_ssp`.

## Metadata

- `fsps_context_get_dims(handle, &nspec, &ntfull, &status)`
- `fsps_context_get_nspec(handle, &nspec, &status)`
- `fsps_context_get_ntfull(handle, &ntfull, &status)`
- `fsps_context_get_nbands(handle, &nbands, &status)`
- `fsps_context_get_nindx(handle, &nindx, &status)`
- `fsps_context_get_nz(handle, &nz, &status)`
- `fsps_context_get_nemline(handle, &nemline, &status)`
- `fsps_context_get_zsol(handle, &zsol, &status)`
- `fsps_context_get_paths(handle, sps_home, sps_len, data_home, data_len, output_home, out_len)`

## Error handling

The driver records the last error or warning internally:

- `fsps_get_last_error(&status, msg, msg_len)`
- `fsps_clear_error()`
- `fsps_set_debug(1)` to print error/warning messages to stdout.

`status == 0` means no error recorded. Non-zero values correspond to driver-specific warning/error codes.

Context-based setters return status directly. A non-zero status indicates an unknown key or invalid handle.

## Array layout

All arrays are passed as flat C buffers. The driver uses Fortran order internally, so treat 2D buffers as column-major (first index varies fastest):

- `spec` is `[nspec, ntfull]`
- `mags` is `[nbands, ntfull]`
- `emlines` is `[nemline, ntfull]`
- `ssp_spec` is `[nspec, ntfull, nz]`

## Thread safety

FSPS is not thread-safe. Do not call these APIs concurrently from multiple threads.

## Minimal example (C)

```c
int h = 0, status = 0;
fsps_context_create(&h, &status);
fsps_context_setup(-1, "mist", "miles", "DL07", h, &status);
fsps_context_set_int(h, "zmet", 1, &status);
fsps_context_compute_ssp(h, spec, mass, lbol, &status);
fsps_context_destroy(h, &status);
```

## Pkg-config

A pkg-config file is provided at fsps.pc. `make install` will install it to your pkg-config directory; otherwise adjust `prefix` and install locations to match your environment.
