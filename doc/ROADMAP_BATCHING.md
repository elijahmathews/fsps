# Roadmap: batched forward modelling, SSP-level parameters and gradients

Planning notes (2026-10-16) for using the fork as a fast forward model in SED fitting:
projected speeds, what the batched path supports for real pipelines, the design of a batched
SSP (IMF and other SSP parameters free per model), and the effort to add gradients. Status
references: `doc/GPU_RESULTS.md` (measured GPU work), `doc/MIXED_PRECISION.md` (FP32
evidence), `doc/AMD_GPU.md` (AMD build).

## 1. Timing projections (full physics)

Configuration: MIST + MILES (5994 pixels), Chabrier IMF (exact integrals); dust attenuation
(birth cloud + diffuse), DL07 dust emission, nebular lines and continuum, IGM at z = 0.5, AGN
torus emission, velocity smoothing at 150 km/s, spectral indices, 159 filter magnitudes;
metallicity varying per model (log Z interpolation between the two bracketing grid
metallicities, as in resolved fitting). Bold = measured; the rest projected. GPU cells are
FP32 / FP64.

| Hardware | Single SSP grid | Batched SSP+CSP per θ, naive | Batched SSP+CSP per θ, folded | Single CSP | Batched CSP per θ |
|---|---|---|---|---|---|
| Ryzen 7 7840U, 1 core | **0.11 s** | ~0.23 s | ~60 ms | ~8 ms (**6.7 ms** at fixed Z) | ~8 ms |
| 7840U, 8 cores (1) | ~15 ms | ~30 ms | ~8 ms | ~1.0 ms | ~1.0 ms |
| AMD Radeon 780M, iGPU of the 7840U (4) | ~0.11 s (2) | ~1.4 / ~30 ms | ~0.6–1 / ~4–6 ms | ~3–5 ms | ~130–260 µs / ~1.5–2.5 ms |
| Intel Arc Pro B70 | ~0.1 s (2) | ~200 / ~3000 µs | ~100 / ~400 µs | ~3–4.5 ms | ~22 / ~180 µs |
| AMD Radeon AI Pro R9700 | ~0.1 s (2) | ~150 / ~2000 µs | ~80 / ~300 µs | ~3–4.5 ms | ~20 / ~120 µs |
| NVIDIA P100 | **0.27 s** (2) | ~400 / ~700 µs | ~100 / ~140 µs | ~3.5 ms | ~18 / ~33 µs |
| NVIDIA V100 | ~0.1 s (2) | ~280 / ~500 µs | ~70 / ~100 µs | ~2.5–3.5 ms | ~13 / ~22 µs |
| NVIDIA A40 | ~0.1 s (2) | ~130 / ~4000 µs | ~60 / ~500 µs | ~2.5–3.5 ms | ~16 / ~130 µs |
| NVIDIA A100 | ~0.1 s (2) | ~160 / ~160 µs (3) | ~40 / ~50 µs | ~2.5 ms | ~8.5 / ~13 µs |

(1) Throughput from one model per core. (2) Limited by scalar work on the host CPU.
(3) FP64 tensor cores run GEMMs at the FP32 rate; TF32 would give ~30–40 µs at ~1e-4 net
error (to be measured). (4) Shares DDR5-5600 memory (~60–70 GB/s) and the 15–28 W budget
with the CPU; practically FP32-only. Measured: 32/32 tests with `OMP_TARGET_OFFLOAD=MANDATORY`;
batched CSP ~450 µs/model (light configuration, 3 kernels offloaded) vs ~1.3 ms on one core.

Columns: *single SSP grid* is one metallicity (a fixed-basis fit needs all ~12 once, ~1.3 s,
plus ~3 s setup). *Batched CSP* is the usual SED fit with fixed SSP parameters. *Batched
SSP+CSP* frees SSP parameters per θ: *naive* builds each θ's SSP grids (two metallicities) and
runs the CSP; *folded* combines the SFH with the template weights first (section 3), so no SSP
grid is built. Measured basis: P100 batched CSP 8.1 µs/model (FP32, cuBLAS, light config);
single CSP on one 7840U core 6.2 ms (full physics, fixed Z) and 6.7 ms with smoothing.
Uncertainty: ±2x NVIDIA, ±3x AMD/Intel (need the full OpenMP port) and the folded column (not
implemented). Other libraries scale with pixel count: C3K ~1.8–2x, BaSeL ~0.4–0.5x; BPASS has
no SSP generation but 2.5x the pixels. FP32 reproduces FP64 to ~3e-7 per pixel and 2e-7 mag.

## 2. Batched CSP against real fitting pipelines

Free parameters in current pipelines: unresolved fits (`logmass`, `logzsol`, `dust2`,
`dust_ratio`, `dust_index`, SFH incl. non-parametric `logsfr_ratio_i`, `fagn`, `agn_tau`,
`gas_logz`, `gas_logu`, `duste_umin/gamma/qpah`, `zred`); resolved fits per pixel (`logmass`,
`logzsol`, `dust2`, `dust_index`, SFH) with some global parameters (`gas_logz`, normalizations).

Per-model support in `compute_csp_batch` (develop, 2026-10-16):

| Parameter | Per model in a batch? | To add |
|---|---|---|
| `logmass` | yes (normalization) | — |
| parametric SFH | yes | — |
| non-parametric SFH | no (tabular SFH is shared context state) | accept an SSP-weight matrix (ages x N) directly: faster and general |
| `dust2`, `dust1`/`dust_ratio`, `dust_index` | yes (`dust_type` shared) | — |
| `duste_umin`, `duste_gamma`, `duste_qpah` | yes | — |
| `logzsol` | no (one metallicity basis) | all metallicities resident; per-model log Z interpolation (+15–30%) |
| `gas_logz`, `gas_logu` | no (nebular rejected) | nebular from the interpolated-Z ionizing spectrum |
| `zred` | no (shared) | per-model resampling to the fixed observed grid, per-model IGM, shared filter matrix (+10–15%) |
| `fagn`, `agn_tau` | no (AGN rejected) | one extra pass |
| `sigma_smooth` | no | one extra pass before the filter integrals; device kernel for single calls |

Also needed: a C/Python entry point taking parameter arrays (struct of arrays, length N) rather
than N parameter structures (building those on the host would rival the GPU time).

Priority: per-model `logzsol`; nebular; SSP-weight-matrix input; per-model `zred`; AGN;
smoothing; C/Python binding. Then both pipelines batch fully.

## 3. Batched SSP-level parameters (IMF and others)

An SSP spectrum is templates . W (W = template weights x ages). Folding the SFH first gives a
CSP = templates . (W . sfh_weights): per θ one or two weight vectors (young/old for birth-cloud
dust) and one batched GEMM; the ionizing photon rate, AGB dust shells and dust emission stay
linear in the templates. Per θ ~3e7 flops instead of ~1e9 for building SSP grids.

IMF parameters touch only the IMF weights: for each age precompute a sparse S (templates x
stars, ~4 entries per star, independent of the IMF; one per metallicity), then per θ:
1. IMF weights for all isochrone points (closed forms; ~2e5 evaluations, ~1e7 flops);
2. W = S . w_IMF(θ) (sparse, ~1e6 flops);
3. fold with the SFH; shared template GEMM; rest of the CSP.
So IMF-free fits cost little more than batched CSPs. IMF *type* best shared per batch.
Isochrone modifications (`sbss`, `fbhb`, `delt`, `dell`, `agb`, `pagb`) change which
templates stars use, so S is rebuilt per θ (the ~40–100 µs folded case). Narrow mass bins:
evaluate as midpoint value x width with a correction (no erf/power-law cancellation) so IMF
weights work in FP32 on consumer GPUs; IMF weights are GPU-compatible since Romberg was
removed (all IMF integrals are closed forms).

## 4. Gradients

Batched CSP = SFH weights -> GEMM -> element-wise attenuation -> energy-balance integrals and
dust emission -> nebular/AGN/IGM -> filter GEMM: linear algebra plus smooth element-wise steps.
Hand-written reverse mode (transposed GEMMs, element-wise adjoints) gives all gradients of a
likelihood for ~2–4x a forward evaluation (A100 ~20–40 µs per θ, full physics), storing a few
intermediate spectra per model; FP32 gradients (~1e-5 relative) suffice for HMC/NUTS.

| Parameter | Difficulty | Notes |
|---|---|---|
| `logmass`, `fagn` | trivial | linear |
| SFH parametric | easy | analytic derivatives of the SFR integrals |
| SFH non-parametric (weight matrix) | trivial | linear in weights; chain to `logsfr_ratio` in the sampler |
| `dust2`, `dust1`/`dust_ratio`, `dust_index` | easy | analytic exp(-tau curve) derivatives; linear energy balance |
| `duste_gamma` | easy | linear mix |
| `duste_umin`, `duste_qpah`, `agn_tau` | easy, kinked | linear grid interpolation: piecewise-constant gradient |
| `gas_logz`, `gas_logu` | moderate, kinked | nebular grids; chain through the ionizing rate (linear in the SSP) |
| `logzsol` | moderate, kinked | kinks at the ~12 grid metallicities |
| `zred` | moderate | smooth distance/IGM; resampling needs the spectral derivative (noisy per pixel, smooth band-integrated) |
| `sigma_smooth` | easy | one convolution with the kernel's sigma derivative |
| IMF parameters (batched SSP) | moderate | closed forms incl. ∫ m^-α ln m; through S . w and the fold |
| isochrone modifications | hard / ill-posed | discontinuous; finite differences or fixed |

Kinks from linear grid interpolation are tolerable for HMC (as in DSPS); an optional smooth
(cubic) interpolation mode removes them at a small, validatable physics change.

Implementation options: (1) hand-written reverse mode in the batched GPU kernels (recommended:
fastest, one codebase, portable OpenMP path; validate each gradient against finite differences
in the FP64 build; build alongside the remaining batch features); (2) AD tools on the Fortran
(Enzyme on flang IR: experimental, immature with GPU offload; Tapenade: no OpenACC/OpenMP, serial
output); (3) a JAX re-implementation of the batched CSP on the fork's precomputed data (free
gradients, but duplicate physics and no portable Fortran path).

Effort: CSP-level parameters moderate (comparable to the remaining batch features); IMF
parameters moderate once the batched SSP exists; isochrone modifications not worth it.
