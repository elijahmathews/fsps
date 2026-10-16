# AMD GPUs with LLVM flang (OpenMP target offload)

Tested on a Radeon 780M (gfx1103, Ryzen 7 7840U APU), Fedora 44: flang 22.1.8, ROCm 7.1.1
(packages flang, clang, clang-tools-extra, lld, llvm, libomp-devel, rocm-device-libs,
rocm-runtime).

```
D=/usr/lib64/rocm/llvm/lib/clang/20/lib/amdgcn/bitcode      # ROCm device libraries
L="['--rocm-device-lib-path=$D','-Xoffload-linker','--library-path=/usr/lib64/amdgcn-amd-amdhsa','-lomp','-lomptarget']"
FC=flang-new CC=clang meson setup build_amd --buildtype=release -Daccel=omp_gpu \
    -Dgpu_arch=gfx1103 -Dkernel_precision=single \
    "-Dfortran_args=--rocm-device-lib-path=$D" "-Dc_args=--rocm-device-lib-path=$D" \
    "-Dfortran_link_args=$L" "-Dc_link_args=$L"
meson compile -C build_amd
env -u HSA_OVERRIDE_GFX_VERSION OMP_TARGET_OFFLOAD=MANDATORY meson test -C build_amd --suite unit --suite regression
```

- Link arguments must use meson's array syntax, and the device-runtime path must be given as
  `--library-path=` (meson reorders `-L`, which separates it from `-Xoffload-linker`).
- `-lomp -lomptarget` are needed explicitly for the shared library.
- `HSA_OVERRIDE_GFX_VERSION=11.0.0` (a common workaround) makes the 780M appear as gfx1100;
  unset it for a gfx1103 build, or build for gfx1100 with it set.

Results (2026-10-16): unit + regression 32/32 with `OMP_TARGET_OFFLOAD=MANDATORY`. Batched CSP
(benchmark_csp_batch, FP32, N >= 100): ~450 us per model on the 780M vs ~1.3 ms on one CPU core.
Only the template product, the batched age sums and the band fluxes have OpenMP target
versions so far; the other batch stages run on the host and data are mapped per call, so this
is a lower bound on what a full OpenMP port gives.
