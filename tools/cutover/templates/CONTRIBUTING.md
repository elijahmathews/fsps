# Contributing to libopensps

Thank you for helping. libopensps carries the physics of FSPS; changes are held to the standard
of a referee: every physics decision needs a reason, and every result a test. `AGENTS.md` is the
detailed guide (it is written for AI coding agents but applies to everyone).

## Build and test

```sh
python3 tools/fetch_data.py ~/opensps-data && export FSPS_DATA_HOME=~/opensps-data
meson setup build --buildtype=release                              # FP64
meson setup build32 --buildtype=release -Dkernel_precision=single  # FP32 kernels
meson compile -C build && meson compile -C build32
meson test -C build --suite unit --suite regression --print-errorlogs
meson test -C build32 --suite unit --suite regression --print-errorlogs
```

## Before opening a pull request

- The unit and regression suites pass in both the FP64 and the FP32 build. Report failures with
  their output; never claim a test passed without running it.
- Regression references (`tests/data/sps_ref_*.bin`, `tests/data/param_ref/`) change only for an
  intended, explained physics change; say why in the commit message.
- Physics changes are compared with legacy FSPS (`tools/legacy_parity/`); every deliberate
  difference is listed with its cause in `tools/legacy_parity/README.md`.
- The single CSP (`compute_csp_scenario`) and the batch (`compute_csp_batch`) stay physically
  identical: a physics change goes into both (`test_csp_batch` compares them).
- C API additions: declared in `include/opensps.h` with the `opensps_` prefix, tested from C
  (`tests/test_c_*.c`), documented in `doc/C_API.md`, with the API version raised.

## Numerical rules

- No `-ffast-math`, nvfortran `-fast` or any flag that reassociates floating point: it deletes
  the compensated sums.
- FP32 kernels follow `doc/MIXED_PRECISION.md`: compensated or pairwise long sums, no
  differences of nearly equal integrals, grid-only factors formed in FP64, FP32 test tolerances
  in terms of the kernel precision.
- Device code (OpenACC, OpenMP target) follows the rules in `AGENTS.md` and `doc/GPU_PLAN.md`;
  GPU builds must match the CPU build (`tools/compare_builds`).

## Git

- Branch from `main`, one topic per pull request, small verified commits.
- Commit messages: an imperative summary line, then what changed and why, with numbers
  (timings, maximum differences) where they matter.
- Pull requests target `main`; CI runs the suites on several compilers. A maintainer merges.
- Add yourself to `AUTHORS` with your first contribution.

Contributions are released under the MIT license of the project (`LICENSE.md`).
