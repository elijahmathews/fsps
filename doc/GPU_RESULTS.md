# GPU results (P100 and A100 sessions)

Running log of the GPU sessions (plan: `doc/GPU_PLAN.md`). Newest status at the top; dated
entries below. The A100 session (2026-10-16) is first; the P100 session's status and entries follow.

## Questions for Elijah

(Decisions that are yours; the agent took the conservative option in each case and kept working.)

1. *(none yet)*

## Status (A100 session, 2026-10-16)

If restarted: read this section, then `git log` of the newest `gpu/a100-*` branch.

**Setup**
- Node `p-gc-3009`, NVIDIA A100-PCIE-40GB, driver 580.178.04 (CUDA 13.0 driver API), 2x Xeon
  Gold 6248R (Cascade Lake, 3.0 GHz), 12 cores allocated. Slurm job 56476035, 9 h.
- Compiler: `module load nvhpc/25.1` (nvfortran/nvc 25.1-0, CUDA 12.6). `module avail nvhpc`
  lists 21.9, 23.1 and 25.1; the driver could run HPC SDK 26.x (CUDA 13), but 25.1 is the newest
  installed and supports cc80. CPU references: nvfortran 25.1 without accel (`build_a100_nvcpu`,
  `build_a100_nvcpu32`) and gfortran 13.4 (`module load gcc/13.4.0`; `build_a100_cpu`,
  `build_a100_cpu32`).
- Builds (all `--buildtype=release`, GPU ones `-Dgpu_arch=cc80`): `build_a100_gpu` / `gpu32`
  (OpenACC FP64 / FP32), `build_a100_gpucb` / `gpu32cb` (`-Dcublas=true`), `build_a100_omp` /
  `omp32` (`-Daccel=omp_gpu`, run with `OMP_TARGET_OFFLOAD=MANDATORY`). The P100 session's
  `build_*` directories target cc60 and are not used.

**Done**
- Step 0: the P100 session's uncommitted change (cuBLAS band product only for >= 16 spectra)
  moved onto branch `gpu/a100-cublas-threshold` from `develop` (9dc0c19). (A stale worktree of
  `develop` from an earlier session, whose directory no longer existed, was pruned.)
  `test_csp_batch` now also runs a 24-model batch: the 12-model batch never reached the
  many-spectra band kernels (OpenACC `band_fluxes` switches at 16 spectra too) nor, with the
  threshold, the cuBLAS product.
- Item 1 (all configurations validated on the A100, branch `gpu/a100-cublas-threshold`):
  unit + regression 32/32 in all ten builds (OpenMP ones with `OMP_TARGET_OFFLOAD=MANDATORY`:
  the first device run of the OpenMP target code); `test_csp_batch` 32/32 everywhere. Case
  matrix (35 cases) against `build_a100_nvcpu`, worst case:

  | Build | spec (912 A-2.5 um) | spec (all) | mag | indx | lines | result |
  |---|---|---|---|---|---|---|
  | OpenACC FP64 (`--rtol 1e-10`) | 1.6e-11 | 3.4e-11 | 1.3e-11 | 8.5e-12 | 3.3e-16 | PASS |
  | cuBLAS FP64 (`--rtol 1e-10`) | 1.6e-11 | 3.4e-11 | 1.3e-11 | 8.5e-12 | 3.3e-16 | PASS |
  | OpenMP target FP64 (`--rtol 1e-10`) | 2.2e-16 | 7.1e-15 | 1.4e-14 | 2.4e-16 | 0 | PASS |
  | OpenACC FP32 | 4.2e-7 | 4.2e-7 | 5.1e-7 | 3.0e-6 | 3.5e-7 | PASS |
  | cuBLAS FP32 | 4.2e-7 | 4.2e-7 | 5.1e-7 | 3.0e-6 | 3.5e-7 | PASS |
  | OpenMP target FP32 | 4.2e-7 | 4.4e-7 | 1.7e-7 | 3.4e-6 | 3.5e-7 | PASS |

  (OpenMP FP64 is near bit-identical to the CPU because in that mode only the template product
  and the batch kernels are offloaded; the OpenACC FP64 builds run the per-star SSP kernel and
  device dust and nebular kernels.)

**In progress / queue** (user-added items first, in order)
- Nebular emission (continuum, lines) and AGN dust emission in `compute_csp_batch`, plus
  indices and line luminosities as outputs (branch `gpu/a100-batch-nebular`, written, being
  validated; line attenuation still missing).
- Then: velocity smoothing in the batch and a device kernel for the single call; per-model
  metallicity (log Z interpolation between resident bases); the full real-world benchmark
  (dust1 0.5, dust2 0.3, DL07, nebular, IGM at z = 0.5, fagn 0.1, indices, sigma_smooth 150);
  then the FP32 speed work, TF32 / FP64 tensor-core experiments and the write-up.

**Open problems**
- Nsight Compute: `RmProfilingAdminOnly: 1` on this node (GPU performance counters for admins
  only); to be confirmed with a run.
- OpenACC FP64 GPU builds: `test_csp_batch` spectra differ from `compute_csp_scenario` by up to
  8.6e-11 (tolerance 1e-10) and mass / L_bol / SFR by 7.7e-13, against <= 6e-16 and 0 in the
  CPU and OpenMP builds. The batch computes SFH weights on the device, the single-output
  scenario on the host (since the P100 session), so the difference is device vs host libm /
  FMA in the SFH weights; not yet investigated (small margin to the tolerance).

## P100 session status (final)

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

- Goal 1 baseline benchmarks and Nsight Systems profile (below): GPU SSP grid is 3-5x slower
  than one CPU core (per-star kernel = 95% of GPU time); CSP call 2.8 ms (FP64) vs 4.0 ms CPU.

- Goal 2 (branch `gpu/p100-ssp-templates`, stacked on `gpu/p100-baseline`): the SSP template
  sum runs in OpenACC builds; the main-library product is a device kernel (3.2 ms per grid).
  SSP grid 2.48 s -> 0.38 s (FP64 GPU), 1.33 s -> 0.39 s (FP32 GPU); nvfortran CPU 0.45 s.
  What remains is host work (weight aggregation, dust-shell AGB corners, small libraries,
  IMF weights). Fixed an FP32 subnormal problem in the AGB template sums (nvfortran only).
  All builds 31/31; case matrix PASS.

- Goal 3 (branch `gpu/p100-batched-csp`, stacked on `gpu/p100-ssp-templates`): batched CSP
  entry point `compute_csp_batch` + `tests/test_csp_batch` + `tests/benchmark_csp_batch`.
  GPU FP32 at N = 10000: 14.6 us per model (OpenACC), 8.1 us (cuBLAS, `-Dcublas=true`),
  vs 5.9 ms for one `compute_csp_scenario` call (1.8 ms FP64). Also fixed on this branch:
  stale device copies of changed settings (dust_type, IGM, ... ignored on the GPU) and
  per-thread device heap allocations (crash at N = 10000). All builds 32/32; case matrix PASS.

- Continuation items 2-3 (branches `gpu/p100-csp-allages`, `gpu/p100-device-mags`, stacked):
  all-ages CSP sum collapsed over outputs (8x); magnitudes of all outputs computed on the
  device (`fsps_batch_photometry`, shared with the batch); single-output SFH weights on the
  host. P100 `compute_csp_scenario`: single call 1.93 -> 0.66 ms (FP64), 5.8 -> 0.66 ms (FP32);
  all-ages call (107 outputs) 94 -> 6.0 ms (FP64), 500 -> 5.0 ms (FP32). All builds 32/32.

- Nebular emission in the CSP call (branch `gpu/p100-device-mags`): no longer re-sends the
  SSP grid per young age; all young ages in a few kernels. P100 CSP call with nebular
  emission 33 -> 2.8 ms (FP64), 22 -> 2.4 ms (FP32). All GPU builds 32/32; case matrix PASS.

- Continuation item 1 (branch `gpu/p100-ssp-small-libs`): the other libraries' templates and
  the dust-shell AGB corners are summed on the device with the main library; per-age buffer
  uploads removed. P100 SSP grid 0.38 -> 0.27 s (FP64), 0.33 -> 0.27 s (FP32); what remains is
  host work, half of it IMF integrals (see open questions). Fixed the template-vs-per-star
  oracle, which compared the template path with itself in nvfortran builds.

- Continuation item 4 (branch `gpu/p100-omp-target`): OpenMP target offload versions of the
  template product, the batched age sums and the band fluxes (`-Daccel=omp_gpu`). Compiles;
  passes 32/32 with `OMP_TARGET_OFFLOAD=DISABLED` (host fallback). Not runnable on the P100
  (see open problems).

- Item 3 remainder (branch `gpu/p100-pset-resident`): `pset` mapped once per CSP call
  (0.668 -> 0.639 ms, FP64); spectral indices of all outputs computed on the device
  (all-ages call with indices 26 -> 6.4 ms); fixed a 25% gfortran FP64 SSP slowdown introduced
  with the FP32 AGB rescaling (now FP32-only; gfortran FP64 output bit-identical to `develop`).
- Session ended before item 5 (final write-up) and before restricting cuBLAS to >= 16 spectra
  (single-spectrum magnitudes with cuBLAS lose accuracy, 5.2e-6, within tolerance, no speed gain).

**Verified on the laptop (2026-10-16, branch tip)**
- gfortran FP64 and FP32: unit + regression 32/32; FP64 outputs bit-identical to `develop` on
  all 36 matrix cases; CPU CSP call with nebular emission 17.9 -> 5.8 ms (Ryzen 7840U).
- nvfortran 26.1 OpenACC multicore FP64 32/32; FP32 failed `csp_batch` at 3-4.7e-6 against a
  2e-6 tolerance (multicore FP32 reductions are plain per-core sums); tolerance now 100 eps(KP).
- `-Daccel=gpu` (cc80,cc90, SDK 26.1) compiles in both precisions.

**In progress**
- Next session (A100): see "Next" in doc/GPU_PLAN.md and the open problems below.

**Open problems**
- *OpenMP offload cannot run on this node.* nvfortran's OpenMP GPU offload needs compute
  capability >= 7.0 ("OpenMP GPU Offload is available only on systems with NVIDIA GPUs with
  compute capability >= cc70"); the P100 is cc60. No offload-capable gfortran/LLVM is installed
  (gcc 12.2-13.4 have no OFFLOAD_TARGET_NAMES; Intel LLVM targets Intel GPUs). `-Daccel=omp_gpu`
  is built for cc70,cc80 and validated with host fallback only; the device run with
  `OMP_TARGET_OFFLOAD=MANDATORY` (which aborts here, as it should) is for an A100 node.
- *Open question (physics/parity, needs a decision):* half of the remaining SSP-grid time on the
  GPU is host IMF integrals (Romberg: IMF weights ~127 ms, remnant masses ~44 ms of ~270 ms).
  Branch `perf/imf-analytic` (exact erf integrals for Chabrier) is marked "needs a parity
  decision"; not merged or used here (conservative choice).
- Branches are stacked (each later one contains the earlier ones' commits), because the
  device-residency fixes on `gpu/p100-baseline` are needed for any GPU run; merge them in order
  (baseline, ssp-templates, ...), or rebase after the first is merged.
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

### 2026-10-15: baseline benchmarks (goal 1)

`build_<x>/benchmark_phases 5 200` (MIST + MILES + DL07, Chabrier; seconds; one CPU core of the
Xeon E5-2680 v4 for the CPU builds; GPU builds use the P100 with the host orchestration on the
same core):

| Metric | gfortran CPU FP64 | nvfortran CPU FP64 | GPU FP64 | GPU FP32 kernels |
|---|---|---|---|---|
| setup | 9.47 | 6.41 | 7.17 | 7.55 |
| first spectrum | 9.96 | 6.86 | 9.71 | 8.94 |
| SSP grid (one call) | 0.477 | 0.424 | 2.475 | 1.332 |
| CSP call | 0.00404 | 0.00481 | 0.00283 | 0.00594 |
| SSP grid, nebular | 0.477 | 0.415 | 2.477 | 1.333 |
| CSP call, nebular | 0.0327 | 0.0461 | 0.0380 | 0.0293 |

(Laptop Ryzen 7840U for reference: SSP 0.147 s, CSP 1.7 ms, CSP+nebular 16 ms; this Broadwell
core is ~3x slower.) Setup is file I/O on GPFS and varies run to run.

Nsight Systems, `nsys profile -t cuda,openacc build_gpu32/benchmark_phases 2 50`
(2+1 SSP grids, 101 CSP calls, 50 with nebular emission), GPU time:

| Kernel / operation | Total | Per call | Share |
|---|---|---|---|
| `accumulate_spectrum` line 711: per-star spectra (one thread per star, seq over pixels) | 5.55 s | 10.4 ms x 535 ages | 94.8% of kernel time |
| `accumulate_spectrum` line 753: weighted sum over stars | 54 ms | 0.10 ms | 0.9% |
| `integrate_csp_step` line 367: CSP age sum, nt x nz per output | 85 ms | 0.84 ms per CSP call | 1.5% |
| nebular kernels (`add_lines_to_spectrum` etc.) | ~100 ms | ~45 us x 2350 | 1.7% |
| dust attenuation + emission | ~11 ms | ~0.1 ms per CSP | 0.2% |
| H2D copies | 0.66 s, 11098 copies | median 1.5 us, max 43 ms | |
| D2H copies | 10 ms, 5303 copies | | |

`cuStreamSynchronize` accounts for 6.2 s of host API time (the host waits on every kernel).
Conclusions: (1) SSP time is entirely the per-star kernel: each thread loops over the full
spectrum of one star, with ~1000-2000 threads per age, so the P100 is latency bound; this is
goal 2 (template GEMM). (2) The CSP call is ~1 ms of GPU work per call plus host orchestration
and many small transfers; the age sum (line 367) is the largest kernel.

### 2026-10-16: goal 2, SSP template sum on the GPU (branch `gpu/p100-ssp-templates`)

**Change.** `generate_ssp_grid` now uses the template path in OpenACC builds as well: per age,
`accumulate_by_templates` aggregates the main-library weights on the host in FP64
(`w_main_all`, templates x ages) and sums the other libraries on the host; after the age loop
`add_main_library_templates` runs on the device: one thread per (pixel, age) sums the used
templates of the resident `speclib(:, zmet, :, :)` in the host order, compensated (Kahan) at
FP32 kernel precision, plain at FP64. Results come back in a private array and are added on
the host, so callers that map `spec_grid` themselves (test_runner does) stay consistent; the
routine also refreshes caller-mapped outputs at exit (`update device ... if_present`). The
per-star OpenACC kernel stays as the reference (`ssp_per_star = .true.`, used by
`physics_oracles`).

**Bug found on the way (nvfortran FP32, CPU and GPU).** The first FP32 GPU run failed
`param_regression_ssp` (ssp_plain: 104 far-IR pixels off by up to 1%). The same error, to the
last digit, appeared in an nvfortran FP32 *CPU* build, and not in gfortran FP32: the template
path had never been compiled by nvfortran in a non-OpenACC build (CI's nvfortran job is
OpenACC multicore, which used the per-star kernel). Dumping the per-age remainder showed the
O/C-rich AGB term: the template path folds the IMF weight in before lambda^2, so at 1-2.7 mm
(lambda^2 ~ 7e14 A^2) the sum is ~1e-37, at the bottom of the FP32 range. gfortran keeps
subnormals; nvfortran's runtime and GPUs flush them. (Ruled out first: optimization level,
FMA, vectorization, and the `-Mflushz/-Mnoflushz` main-program flag, which did not change the
result.) Fix: those sums carry an exact factor 2^64, removed after lambda^2 (bit-identical in
FP64). This is rule 5 of `doc/MIXED_PRECISION.md` (keep values in range) for an intermediate.

**Validation.** Unit + regression 31/31 in all four builds (gfortran FP64 and FP32, GPU FP64 and
FP32). Case matrix, worst over 36 cases:

| Comparison | spec (912 A-2.5 um) | spec (all) | mag | indx | lines |
|---|---|---|---|---|---|
| FP64 GPU vs nvfortran CPU (rtol 1e-10) | 1.6e-11 | 1.2e-10 | 1.3e-11 | 2.3e-11 | 1.2e-12 |
| FP32 GPU vs gfortran CPU FP64 (defaults) | 3.1e-7 | 4.1e-7 | 2.0e-7 | 4.1e-6 | 2.3e-7 |

**Speed** (`benchmark_phases 5 200`, seconds per SSP grid):

| Build | before | after |
|---|---|---|
| GPU FP64 | 2.475 | 0.383 |
| GPU FP32 | 1.332 | 0.387 |
| nvfortran CPU FP64 (reference) | | 0.446 |

Nsight: `add_main_library_templates` kernel 3.16 ms per grid (FP32, Kahan, one thread per
pixel x age, 1292 templates x 107 ages x 5994 pixels). The SSP call is now host-bound. gprof
of the FP32 GPU build (`-pg`, 7 grids): per grid ~0.2 s in `accumulate_by_templates`
(aggregation, the per-age small-library sums `add_column`, and above all the dust-shell AGB
corners, where each corner slot does several full-spectrum passes per age) and ~0.09 s in the
IMF weights (`get_imf_value` / Romberg `refine_trapezoid`). Next for SSP: move the per-age
small-library and dust-corner sums to the device as products over all ages (continuation item
1); the IMF weights need the analytic integrals (branch `perf/imf-analytic`) to leave the host.
A batch of grids with different IMF parameters would currently be host-bound at ~0.38 s per
grid, so it was not timed separately.

### 2026-10-16: goal 3, batched CSP (branch `gpu/p100-batched-csp`)

**New entry point** `compute_csp_batch(ctx, psets, mags, status [, spec, mass, lbol, mdust, sfr])`
in `src/physics/fsps_csp_batch.f90` (`compute_csp_scenario` unchanged). N single-age models
(`tage > 0`) against the resident SSP basis (one metallicity), sharing zred / igm_factor;
SFH, dust attenuation (dust1/2/3, indices, birth-cloud age, frac_obrun, frac_nodust, curve
parameters) and dust emission (U_min, q_PAH, gamma) vary per model. Status 2 for nebular
emission, AGN, smoothing, redshifted colors, light ages, tage <= 0 (decided to keep v1 to the
stages that are matrix products or element-wise; these are listed as open items). Stages:

1. SFH weights per model on the device in FP64 (`sfh_weights_core`, one thread per model;
   factored out of `compute_sfh_weights` with caller-provided scratch), split into young/old
   at each model's `dust_tesc`, cast once to KP; mass and L_bol from them (FP64).
2. Age sums: one product basis(pixels x ages) . [W_young W_old](ages x 2N).
3. Attenuation + absorbed light density: element-wise over (pixel, model), as in `fsps_dust`.
4. Dust emission: per-model reductions (gang per model), Draine & Li template interpolated per
   model on the device from the resident grid, added with the energy-balance normalization.
5. Magnitudes: observed-frame interpolation, then per-(band, model) filter integrals over each
   band's support (Kahan at FP32), or with `-Dcublas=true` a split-K product (256-pixel blocks,
   blocks summed in FP64).

**Validation.** `tests/test_csp_batch.f90` (regression suite, all builds): 12 models x 4
configurations (Calzetti z = 0; power law z = 0.5 with IGM; Cardelli z = 0; Kriek & Conroy
z = 1.2; SFH types 0, 1, 4, 5, bursts, constant, sf_start/trunc/slope, dust emission on)
against `compute_csp_scenario` model by model:

| Build | spectrum (rel.) | magnitudes | mass, L_bol, SFR | dust mass |
|---|---|---|---|---|
| gfortran CPU FP64 | 0 (bit-identical) | 1.4e-14 | 0 | 0 |
| GPU FP32, OpenACC | 2.6e-7 | 2.7e-7 | 1.1e-16 | 1.1e-7 |
| GPU FP32, cuBLAS | 6.8e-7 | 4.4e-7 | 1.1e-16 | 1.1e-7 |

All builds pass unit + regression (32 tests, including the new one); case matrix unchanged
(FP64 GPU = nvfortran CPU to 1.6e-11; FP32 GPU vs FP64 CPU 3.1e-7 / 2e-7 mag).

**Bugs found while building it** (all fixed on this branch):
- *Stale settings on the device.* Kernels read scalar settings (`dust_type`,
  `add_igm_absorption`, ...) from the device copy of `ctx`, written only by `move_to_device`.
  A setting changed afterwards had no effect on the GPU: the test's IGM case came out without
  IGM absorption, and all dust types used the first one. This affects `compute_csp_scenario`
  too (any caller that changes a setting after moving the context). Fix: the setters mark the
  context dirty and the compute entry points refresh the device copy
  (`fsps_context_sync_device_settings`).
- *Device heap exhaustion.* The batch crashed at N = 10000 (illegal address; fine at 1000;
  `compute-sanitizer` pointed at local memory in `compute_ssp_weights`; a larger stack did not
  help, a 256 MB device heap did). Device code allocated temporaries on the 8 MB device heap
  in every thread: the context's state arrays are pointers, so passing them to `contiguous`
  dummies (`find_interval`, ...) needed a runtime contiguity copy, and the tabular SFH code
  passed the strided row `sfh_tab(1, :)`. Fix: state pointer arrays are declared `contiguous`
  (they always point to whole allocated arrays) and the table search runs in place. A scan of
  the PTX of every module (`-gpu=keep`, grep `malloc`) now finds no device allocations.
- An out-of-bounds read in the new observed-frame kernel at z = 0 (found with a
  `-fcheck=all` gfortran build).

**Speed** (`build_x/benchmark_csp_batch 10000 200`; microseconds per model; tau, age, dust1/2
and U_min vary; `mags` = magnitudes only, `spec` = also copying the spectra back):

| N | GPU FP32 OpenACC | GPU FP32 cuBLAS | GPU FP64 OpenACC | GPU FP64 cuBLAS |
|---|---|---|---|---|
| one `compute_csp_scenario` call | 5837 | 5754 | 1832 | 1916 |
| 1 | 5859 | 7932 | 6827 | 10267 |
| 100 | 89 | 119 | 110 | 119 |
| 1000 | 51 | 22 | 64 | 22 |
| 10000 (mags) | 41.9 | 10.9 | 57.6 | 12.5 |
| 10000 (spec) | 62.2 | 31.4 | 86.9 | 40.8 |

(The FP32 single call is slower than FP64 because its band integrals run on the host,
compensated and serial per band; `doc/MIXED_PRECISION.md` notes the same on the CPU.)
Against the forecast of ~4 us per model (A100, FP32, large N), the P100 reaches 10.9 us with
cuBLAS. Where the rest goes, N = 10000, OpenACC FP32 before the last kernel changes (timers in
the driver): host stage (indices, SFR) and transfers small; device stages: age sums +
attenuation + dust emission ~128 ms; band fluxes ~237 ms; host magnitude assembly 22 ms;
spectra copy-back ~200 ms (240 MB pageable D2H).

**Kernel work after the first version** (OpenACC FP32 unless noted; us per model, mags only):

| Step | N = 1 | N = 100 | N = 1000 | N = 10000 |
|---|---|---|---|---|
| first version (host SFH weights) | 5753 | 99 | 63 | 57.5 |
| SFH weights on the device, no device heap allocations | 5859 | 89 | 51 | 41.9 |
| band fluxes per (band, model block), transposed observed spectrum | 6190 | 79 | 30 | 21.4 |
| age sums: one pass (young or old per age), 4 models per thread | 6665 | 79 | 25 | 17.4 |
| keep transmission + dust template from the attenuation pass | 6070 | 73 | 23 | 15.4 |
| cache grid-only data per redshift (on the device) | 1688 | 31 | 16 | 14.6 |
| same, cuBLAS FP32 | 1076 | 24 | 9.4 | **8.1** |

Profile of the OpenACC FP32 build at N = 10000 before the last two steps (per call): age sums
81 ms (before the one-pass change), band fluxes 36 ms, attenuation 26 ms, dust reductions
18 ms, dust emission 13 ms, SFH weights 3.3 ms, observed frame 4 ms. With spectra copied back
(`spec`), add ~20 us per model at N = 10000: 240 MB through pageable host memory.

### 2026-10-16: continuation items 2 and 3, compute_csp_scenario on the device

Timing harness: `compute_csp_scenario` with tage = 0 (107 outputs, tau model, dust), and
`benchmark_phases` for single calls; ms per call on the P100 unless noted.

| Change | single call FP64 / FP32 | all-ages FP64 / FP32 |
|---|---|---|
| before (branch `gpu/p100-batched-csp`) | 1.93 / 5.8 | 86 / 500 |
| item 2: multi-output age sums, one thread per (output, pixel) | 1.93 / 5.8 | 89 / 499 (kernel 4.6 -> 0.57 ms) |
| item 3: magnitudes of all outputs on the device | 1.93 / 2.0 | 6.1 / 5.4 |
| single output: SFH weights on the host; band reductions for < 16 spectra | **0.66 / 0.66** | **6.0 / 5.0** |
| nvfortran CPU, one Broadwell core (reference) | 4.8 | 278 |

- Item 2: the multi-output spectra and emission-line sums in `integrate_csp_step` ran one gang
  per output with the ages sequential and a loop-carried accumulator (`-Minfo`: serialized).
  They now use the single-output kernel collapsed over (output, pixel); same summation order,
  same results. (The batched GEMM form of the plan is what `compute_csp_batch` does; here the
  weights differ per output only in the age range, and the collapsed kernel is 0.57 ms.)
- Item 3: the per-output host loop computed magnitudes on the host after a device round trip
  of each spectrum (FP32: compensated serial band sums, 4.6 ms per output). Now, without
  smoothing, redshifted colors or light-weighted ages, `compute_magnitudes_columns` computes all
  outputs' magnitudes from `out_csp_spec` on the device (OpenACC builds only; CPU results
  unchanged), and indices take each spectrum directly. The batch grid cache had to be released
  after the call for contexts mapped by hand (unit tests), which never reach
  `remove_from_device`.
- Single output: `compute_sfh_weights` ran on one device thread (0.65 ms); the host takes
  ~20 us. Band sums for one spectrum ran one thread per band (0.65 ms); a gang reduction per
  (band, spectrum) is used below 16 spectra (pairwise tree, as the dust integrals).
- Profile of the remaining single call (FP32, 0.66 ms): device kernels ~0.3 ms (dust
  attenuation and emission ~0.1 ms, age sum 50 us, post-processing 35 us, ...), the rest
  launches, small transfers and host packing. With nebular emission: 22-33 ms per call, which
  is the next item (47 kernel launches per call and host work).

### 2026-10-16: nebular emission in the CSP call

`benchmark_phases` "CSP call, nebular" (gas_logu varies per call): 33 ms (FP64) / 22 ms (FP32)
on the GPU, against 16 ms on the laptop CPU. Two causes:

1. `process_ionizing_radiation` ran `!$acc update device(spec_in)` for every young age, and
   `spec_in` is the whole SSP grid: ~50 x 5 MB host-to-device per call (`apply_nebular_emission`
   did it once more on entry). Callers either map it (the resident basis) or rely on `pcopyin`,
   so the updates were redundant. Removed: 33 -> 7.5 ms.
2. Per young age, ~6 kernel launches and a synchronous copy-back of Q. `nebular_all_ages` (OpenACC
   only; the CPU loop is unchanged) does EUV escape, Q per age (gang reduction), continuum, line
   luminosities and line profiles each in one kernel over (pixel or line, age), with Q kept on
   the device: 7.5 -> 2.8 ms (FP64), 2.4 ms (FP32).

Case matrix after this: FP64 GPU vs nvfortran CPU, csp_nebular 7.9e-14 per pixel, lines 2e-16;
FP32 GPU vs FP64 CPU, csp_nebular 3.1e-7 per pixel, 3.5e-7 mag, lines 2.3e-7 (all cases PASS;
worst magnitude 4.6e-7, csp_ssp_offgrid, from the single-spectrum band reduction).

Current `benchmark_phases 5 200` on the P100 (s): SSP grid 0.38 / 0.33, CSP call 0.00066 /
0.00066, SSP grid with nebular 0.40 / 0.33, CSP with nebular 0.0028 / 0.0024 (FP64 / FP32).

### 2026-10-16: continuation item 1, the rest of the SSP template sum on the device

Timers in `generate_ssp_grid` (FP32 GPU, per grid, before): load 9 ms, isochrone physics (IMF
weights) 127 ms, buffer upload + integrated properties 56 ms, template aggregation + other
libraries + dust corners 133 ms (of which: aggregation 56, other libraries ~10, dust corners
75), main-library kernel 3-5 ms; total 330 ms.

1. Buffer uploads: the isochrone buffer is only needed on the device by the per-star kernel;
   with the template path the 9 per-age uploads and the device reductions of
   `compute_integrated_properties` are skipped (`if(ssp_per_star)` clauses). The remaining
   44 ms there is `add_remnant_mass` (an IMF integral, host).
2. Other libraries (WMBasic, post-AGB, WN/WC, carbon stars, O/C-rich AGB): in OpenACC builds
   their template weights are deferred as extra rows of the all-ages weight matrix and summed
   in the main-library kernel (O/C-rich apart, times lambda^2, with the exact 2^64 factor).
3. Dust-shell AGB corners: each age's corners are recorded (transfer-function indices and the
   corner's template weights, CSR-style) and the device kernel rebuilds each corner's template
   sum (smoothed beyond 3 um; O/C-rich below it times lambda^2) times its transfer function,
   compensated at FP32 as on the host. 75 ms -> ~0 on the host.

CPU builds keep the per-age host sums (unchanged results). SSP grid, P100: FP64 0.38 -> 0.27 s,
FP32 0.33 -> 0.27 s (nvfortran CPU: 0.45 s; gfortran CPU 0.48 s, one Broadwell core).

**Test gap found:** `physics_oracles` (template sum = per-star sum) set `fsps_ssp::ssp_per_star`
by assignment from the test program. With nvfortran, a program that assigns a module variable
of a shared library changes its own copy (the library keeps its own): both runs used the
template path, the check printed 0.00 and passed vacuously in every nvfortran build
(including CI's multicore job, where it compared per-star with per-star). New setter
`fsps_ssp_set_per_star`; the oracle now compares the two paths on the GPU: FP64 8.9e-12 per
pixel (tolerance 1e-10), FP32 within KP_RTOL above 1e-6 of the peak (the printed maximum,
1.8e5, is from pixels near FP32 underflow, ~1e-30 of the grid maximum, which the GPU flushes
to zero).

**nvfortran quirk:** assigning `dusty%slot_ptr(1:n)` (a component of an optional
`intent(in)` derived-type dummy) to an unallocated local allocatable left it unset (garbage
indices, illegal address in the kernel); explicit `allocate` fixed it.

Current `benchmark_phases 5 200` (s), P100, FP64 / FP32: SSP 0.270 / 0.270, CSP 0.00067 /
0.00067, SSP + nebular 0.269 / 0.278, CSP + nebular 0.0027 / 0.0024. Case matrix: FP64 GPU vs
nvfortran CPU 1.6e-11 per pixel, 1.3e-11 mag; FP32 GPU vs FP64 CPU 3.8e-7, 4.6e-7 mag. PASS.

### 2026-10-16: continuation item 4, OpenMP target offload

New build mode `-Daccel=omp_gpu` (nvfortran `-mp=gpu`, `-DOMP_TARGET`, OpenACC off): OpenMP
`target teams distribute parallel do` versions of the main-library template product
(`main_templates_omp` in `fsps_ssp.f90`), the batched age sums and the band fluxes (same
arithmetic as the OpenACC kernels; arrays mapped per call with `map(to:/from:)`, no device
residency yet). Everything else runs on the host in this mode.

```
FC=nvfortran CC=nvc meson setup build_omp32 --buildtype=release -Daccel=omp_gpu \
    -Dgpu_arch=cc70,cc80 -Dkernel_precision=single
OMP_TARGET_OFFLOAD=MANDATORY build_omp32/test_csp_batch   # aborts on the P100 (cc60)
OMP_TARGET_OFFLOAD=DISABLED meson test -C build_omp32 --suite unit --suite regression  # 32/32
```

With `-Dgpu_arch=cc60` nvfortran 25.1 refuses the target regions: OpenMP GPU offload needs
cc70+. So on this node the OpenMP path is validated numerically through host fallback only
(csp_batch: spectra 5.9e-7, magnitudes 1.5e-7 against `compute_csp_scenario`); the device run
is pending an A100 (V100 or newer).
