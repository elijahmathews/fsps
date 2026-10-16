# GPU results (P100 session)

Running log of the first GPU session (plan: `doc/GPU_PLAN.md`, "On a P100 node"). Newest
status at the top; dated entries below.

## Status

**Done**
- Environment checks: GPU visible (`nvidia-smi`: Tesla P100-PCIE-12GB, driver 580.178.04);
  `git push --dry-run origin HEAD` succeeds (HTTPS through the `gh` credential helper).
- Build fix: Fortran executables now link with the Fortran compiler (nvc 25.1 does not add
  the Fortran runtime). Branch `gpu/p100-baseline`.
- Builds: `build_cpu` (gfortran 13.4, FP64), `build_gpu` (nvfortran 25.1, cc60, FP64),
  `build_gpu32` (FP32 kernels).
- Device data-residency fixes (context lifecycle, AGN kernel, one unit test): all builds pass
  31/31 (`build_cpu`, `build_gpu`, `build_gpu32`).
- Case matrix (tools/compare_builds): FP64 GPU = nvfortran CPU to 1.6e-11 per pixel; FP32 GPU
  vs FP64 CPU 3.3e-7 per pixel, 2e-7 mag (same as the laptop FP32 numbers). PASS.

**In progress**
- Goal 1 benchmarks: `build_*/benchmark_phases 5 200` for cpu, gpu, gpu32, then nsys.

**Open problems**
- gfortran 13.4 and nvfortran 25.1 CPU builds differ by up to 8.5e-10 per pixel (4e-10 mag)
  for the power-law IMFs (Salpeter, Kroupa, Dave); Chabrier agrees to 2e-13. Not a GPU issue
  (the GPU matches nvfortran CPU to 2e-13 on those cases), but an FP64 GPU build compared
  against a *gfortran* CPU build fails `--rtol 1e-10` there; compare FP64 GPU builds against
  an nvfortran CPU build (`build_nvcpu`), or use `--rtol 1e-9`.

## Entries

### 2026-10-15: node and toolchain

- Node `p-gc-3150` (ICDS), 2x Xeon E5-2680 v4 (5 cores allocated), 251 GB RAM.
- GPU: Tesla P100-PCIE-12GB, driver 580.178.04 (CUDA 13.0 driver API).
- Compilers: `module load nvhpc/25.1` (nvfortran/nvc 25.1-0, Nsight Systems 2024.7.1);
  CPU reference built with `module load gcc/13.4.0` (the system gfortran is 8.5).
  meson 1.10.1, numpy 1.26.4.
- Build fix: meson linked the test executables with `nvc` (libfsps has one C source) and
  nvc 25.1 does not link `libnvf`, so every Fortran program failed with
  `undefined reference to pgf90io_*`. Added `link_language : 'fortran'` to every Fortran
  executable (the library already had it). CPU builds are unaffected.
- CPU (`build_cpu`, gfortran 13.4): `meson test --suite unit --suite regression`, 31/31 OK, 83 s.

### 2026-10-15: device data-residency bugs (GPU FP64 suite 25/31 -> 31/31)

First run of `meson test -C build_gpu --suite unit --suite regression` on the P100:
6 failures, all host/device mapping bugs that OpenACC multicore (CI) cannot see:

1. `fsps_destroy` never removed a context from the device. Its freed host arrays kept their
   present-table entries; the next context's arrays, allocated at overlapping host addresses,
   then failed with "variable in data clause is partially present"
   (`ctx%state%gaussnebarr`, `fsps_context_move_to_device`). Broke param_regression_* and
   physics_oracles. Fix: `ctx%on_device` flag; destroy, setup (which repoints and reallocates
   state arrays) and `sps_takedown` call `fsps_context_remove_from_device` first; a second
   `move_to_device` refreshes instead of stacking reference counts.
2. `fsps_context_remove_from_device` on a context never mapped, as the first OpenACC operation
   in the process, segfaults inside the 25.1 runtime (`__pgi_uacc_select_devid`). Hit by the
   photometry unit tests' teardown-before-setup. Fix: it returns when `.not. ctx%on_device`.
3. `apply_agn_dust_emission` had `present(settings)` for the caller's host `pset`. Fix:
   `copyin(settings)`; the kernel also now waits for queue 1 (`wait(1)`), since the preceding
   CSP post-processing kernel writing the spectrum runs `async(1)` and the AGN kernel ran on
   the default queue (a race on a GPU, harmless on multicore).
4. Unit test `test_remnants_ignore_appended_stragglers` passed a host-only isochrone buffer to
   `compute_integrated_properties` (device reductions). Fix: map it, as the neighbouring
   remnant test already does.

Results after the fixes (meson test, unit + regression):

| Build | Result | Wall time (4 jobs) |
|---|---|---|
| build_cpu (gfortran 13.4, FP64) | 31/31 | 83 s (5 jobs) |
| build_gpu (nvfortran 25.1, cc60, FP64) | 31/31 | ~170 s |
| build_gpu32 (FP32 kernels) | 31/31 | |

### 2026-10-15: case matrix, GPU vs CPU (tools/compare_builds)

```
tools/compare_builds/run_cases.sh build_cpu   out/cpu      # gfortran 13.4, FP64
tools/compare_builds/run_cases.sh build_nvcpu out/nvcpu    # nvfortran 25.1, no accel, FP64
tools/compare_builds/run_cases.sh build_gpu   out/gpu
tools/compare_builds/run_cases.sh build_gpu32 out/gpu32
python3 tools/compare_builds/compare_builds.py out/gpu out/nvcpu --rtol 1e-10   # PASS
python3 tools/compare_builds/compare_builds.py out/gpu32 out/cpu                # PASS
```

Worst case over the 36 cases:

| Comparison | dlogM | dlogL | spec (912 A-2.5 um) | spec (all) | mag | indx | lines |
|---|---|---|---|---|---|---|---|
| FP64 GPU vs nvfortran CPU | 5.4e-13 | 7e-12 | 1.6e-11 | 1.2e-10 | 1.3e-11 | 2.3e-11 | 1.2e-12 |
| FP32 GPU vs gfortran CPU (FP64) | 1.2e-12 | 1.8e-11 | 3.3e-7 | 3.4e-7 | 2.0e-7 | 4.1e-6 | 2.3e-7 |

Note the FP64 GPU path computes SSPs per star (the OpenACC kernel) while the CPU builds use the
template sum; they agree to 2e-13 per pixel on the SSP cases.

FP64 GPU vs gfortran CPU fails `--rtol 1e-10` on four SSP cases (ssp_salpeter 5.1e-10,
ssp_kroupa 4.7e-10, ssp_dave 4.7e-10, ssp_kroupa_zlow 8.5e-10). Isolated with two more nvfortran
builds on the same commit (`build_nvcpu`: no accel; `build_mc`: `-Daccel=multicore`): both
show the same 5e-10 against gfortran, and the GPU matches both to 2e-13. It is a
gfortran-vs-nvfortran difference in the power-law IMF weights (Chabrier agrees to 2e-13),
not a device effect.
