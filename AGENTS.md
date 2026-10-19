# Agent guide

Instructions for AI coding agents (and a quick reference for people). Read this first; the
linked documents hold the detail. Keep this file short and current: when a rule here changes,
change it here.

## The project

A modernized fork of FSPS (Flexible Stellar Population Synthesis; Conroy et al.), in Fortran
2008+, built with Meson. Goals: physics that is correct and traceable to legacy FSPS or to
independent identities, and fast forward models (FP32 kernels, GPU offload, batched CSPs for
SED fitting).

- `src/core` (contexts, types, precision), `src/physics` (IMF, isochrones, SSP, CSP, dust,
  nebular, AGN, batched CSP), `src/spectra` (libraries, filters, indices, batched photometry),
  `src/math`, `src/io`, `src/bindings` (C API; `include/`).
- Two real kinds (`src/core/fsps_precision.f90`): `WP` is always FP64 (setup, isochrones, IMF
  and SFH weights, the public API); `KP` is the kernel precision, FP64 by default and FP32
  with `-Dkernel_precision=single`.
- The user is an astronomer and an SPS expert: explain software matters plainly; justify every
  physics decision with a reason that would satisfy a referee.

## Build and test

```
export SPS_HOME=$PWD FSPS_DATA_HOME=$PWD
meson setup build --buildtype=release                              # CPU, FP64
meson setup build32 --buildtype=release -Dkernel_precision=single  # FP32 kernels
meson compile -C build
meson test -C build --suite unit --suite regression --print-errorlogs
```

GPU builds (details and per-machine notes: `doc/GPU_PLAN.md`, results: `doc/GPU_RESULTS.md`):

- OpenACC: `FC=nvfortran CC=nvc meson setup b -Daccel=gpu -Dgpu_arch=cc80` (add
  `-Dcublas=true` for cuBLAS batch products); `-Daccel=multicore` runs the same code on CPU cores.
- OpenMP target: `-Daccel=omp_gpu` with nvfortran or LLVM flang (AMD build recipe:
  `doc/AMD_GPU.md` on branch `gpu/780m-omp`). Run with `OMP_TARGET_OFFLOAD=MANDATORY` so a silent host
  fallback fails.

Other tools: `tests/benchmark_phases` (SSP/CSP timings), `tests/benchmark_csp_batch`,
`tools/compare_builds/` (one build against another over the 36-case matrix, e.g. GPU vs CPU),
`tools/legacy_parity/` (this fork against legacy FSPS on identical data).

## Correctness gates

- All of `--suite unit --suite regression` pass in FP64 and FP32 builds before a PR.
- Regression references (`tests/data/sps_ref_*.bin`, `tests/data/param_ref/`) are regenerated
  only for an intended, explained change (`tests/generate_test_data`,
  `test_param_regression ... --generate`); say why in the commit message.
- Physics changes are checked against legacy FSPS (`tools/legacy_parity/run_matrix.sh`); every
  deliberate difference is listed with its cause in `tools/legacy_parity/README.md`.
- GPU builds match the CPU build: FP64 to ~1e-10 per pixel, FP32 to ~5e-7
  (`tools/compare_builds`).
- Never claim a test passed without having run it; report failures with their output.

## Hard rules

- Never `-ffast-math`, nvfortran `-fast`, or any flag that reassociates floating point: it
  deletes the compensated sums.
- FP32 kernels follow `doc/MIXED_PRECISION.md` ("Rules for FP32 kernels"): compensated or
  pairwise long sums; no differences of nearly equal integrals; no FP32 logarithms that are
  exponentiated later; grid-only factors formed in FP64; physical units; floors of
  `tiny(1.0_KP)`. Test tolerances in FP32 builds use `KP` (e.g. 100 eps of `KP`).
- Device code (OpenACC and OpenMP target):
  - `!$acc routine seq` on every routine called from device code;
  - no libm calls (e.g. `log_gamma`) without a device-safe fallback;
  - `contiguous` dummies (no runtime copy-in/out);
  - no array constructors or other temporaries (LLVM flang mallocs them on the device);
  - map a derived type either whole or by component, never both;
  - map everything a kernel touches explicitly: NVIDIA GPUs may tolerate a missing map,
    AMD GPUs fault.
  More in `doc/GPU_PLAN.md` ("Rules learned so far").
- Keep `compute_csp_scenario` and `compute_csp_batch` physically identical; a physics change
  goes into both and the batch test (`test_csp_batch`) checks them against each other.

## Git workflow

- Work on a topic branch from `develop`; commit small, verified steps and push often:
  `git push origin <branch>` (GitHub). Never push to other remotes.
- Commit unsigned (signing needs a hardware key): `git -c commit.gpgsign=false commit`. The
  user is the author; end the message with a `Co-Authored-By:` trailer naming the agent.
- Open pull requests against `develop` (CI runs only there); do not merge them, the user does.
  Stacked branches: say the merge order in the PR.
- Commit messages: an imperative summary line, then what changed and why, with numbers
  (timings, max differences) where they matter.

## Machines

- ICDS cluster: `module load nvhpc/25.1` (CUDA 12; P100 is `cc60`, A100 `cc80`). HPC SDK 26
  dropped Pascal and Volta.
- The user's laptop has a Radeon 780M iGPU (gfx1103) that also drives the display. A GPU fault
  or hang there resets the GPU and kills the desktop session. Do not run GPU offload tests on
  it without the user's go-ahead; debug with host runs (`OMP_TARGET_OFFLOAD=DISABLED`) and by
  reading the device IR first. Unset `HSA_OVERRIDE_GFX_VERSION` for gfx1103 builds.
