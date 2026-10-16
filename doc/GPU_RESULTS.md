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
- CPU suite: 31/31 pass.

**In progress**
- GPU FP64 suite: 25/31 pass; failing: unit_ssp, unit_photometry (SIGSEGV),
  param_regression_{ssp,csp,lib}, physics_oracles. Investigating (device data residency?).

**Open problems**
- (none yet)

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
