# Ecosystem Impact Analysis (Build, CI, Packaging)

## 1. Dependency options

## Preferred
- Optional HDF5 dependency (`-Ddata_backend=hdf5`) with default `legacy`.
- Keep zero-new-dependency path for conservative environments.

## Alternatives
- NetCDF4 backend (also HDF5-backed) can be added later.
- Pure-manifest/raw backend for minimal builds, but not recommended as canonical runtime.

## 2. Meson changes (`meson.build`, `meson.options`)

### Add options
- `data_backend` (combo: `legacy`, `hdf5`, `auto`) default `legacy`
- `with_hdf5` (boolean) default `false`
- `strict_schema` (boolean) default `true`
- `hdf5_fortran_module_dir` (string) default `''`
- `hdf5_fortran_link_args` (array) default `[]`
- `hdf5_fortran_include_args` (array) default `[]`

### Add dependency wiring
- Attempt layered discovery strategy:
  1. `dependency('hdf5-fortran', required : false)`
  2. `dependency('hdf5', modules : ['fortran'], required : false)` (where supported)
  3. `dependency('hdf5', required : false)` + manual module include args
  4. fallback to user-provided `hdf5_fortran_include_args` + `hdf5_fortran_link_args`

If `-Dwith_hdf5=true` and all discovery paths fail, emit a clear configure error with exact override instructions. If `-Ddata_backend=auto`, fall back to `legacy` with warning.

- Add new source files for backend abstraction and HDF5 implementation.
- Conditional compile definitions:
  - `-DFSPS_HAS_HDF5=1`
  - `-DFSPS_STRICT_SCHEMA=1` when enabled.

### Recommended Meson fallback behavior
- `data_backend=legacy`: ignore HDF5 probing entirely.
- `data_backend=auto`: use HDF5 if probing succeeds; otherwise continue with legacy backend.
- `data_backend=hdf5`: require successful HDF5 Fortran resolution, else hard fail.

This prevents brittle CI failures across gfortran/nvfortran and distro differences.

### Installation impact
- If distributing standardized data file(s), install under `${datadir}/fsps/standard/`.
- Keep existing `install_subdir('data', ...)` behavior during migration.

## 3. CI changes (`.github/workflows/ci.yml`)

Add matrix variants:
1. **legacy mode** (existing tests unchanged)
2. **hdf5 mode**
   - install HDF5 dev package
   - configure with `-Dwith_hdf5=true -Ddata_backend=hdf5`
  - add one job variant that uses manual Meson overrides (`hdf5_fortran_include_args/link_args`) to verify fallback path

Regression strategy:
- Run same regression suite for both modes.
- Optional parity check job compares selected outputs (`.spec/.mags`) with tolerance.

Container notes:
- NVHPC and Flang jobs may need explicit HDF5 Fortran package availability; if unavailable, keep those jobs on `legacy` mode only.
- For NVHPC specifically, module include path injection may be required even when linker discovery succeeds.

## 4. RPM impact (`packaging/fsps.spec`)

## BuildRequires additions (conditional preferred)
- `hdf5-devel` (or distro equivalent with Fortran modules)

## Spec toggles
- `%bcond_with hdf5` to control optional dependency.
- Pass `%meson -Dwith_hdf5=%{?with_hdf5:true}%{!?with_hdf5:false}`.

## Subpackage considerations
- Optional package: `fsps-data-standard` for standardized HDF5 datasets, separated from `libfsps` runtime.
- Maintain current stance (external data payload optional) for lightweight packaging.

## 5. Debian impact (`debian/control`, `debian/rules`)

### `debian/control`
- Add optional Build-Depends when enabling HDF5 build profile:
  - `libhdf5-dev` or `libhdf5-fortran-dev` (distribution-specific naming)

### `debian/rules`
- Pass Meson option flags for backend selection.
- Keep current removal of `/usr/share/fsps` if binary package remains data-free.

## 6. Developer tooling impact
- Add schema validation command/tool in `tests/` or `tools/`.
- Add conversion pipeline CI check: legacy->FSDS conversion smoke test.

## 6.1 Memory and performance impact of alpha-axis expansion
Adding an explicit `afe` axis can multiply spectral-library memory usage by approximately `N_afe`:

$$
M_{new} \approx M_{old} \times N_{afe}
$$

Implications:
- Larger RAM requirements during initialization and cache residency.
- Higher I/O pressure if full cubes are loaded eagerly.

Therefore, chunking/slicing is a requirement (not optional):
- Enable slice reads for only needed `(z, afe)` slabs when possible.
- Prefer chunk layouts that keep `lambda` contiguous for compute kernels.
- Consider lazy loading or staged cache population for very large model families.
- Include memory regression checks in CI for `hdf5` mode to prevent accidental eager full-cube loads.

## 7. Risk assessment
- **Low risk** if HDF5 is optional and legacy remains default.
- **Medium risk** for compiler matrix where HDF5 Fortran modules differ.
- **Medium-high risk** for memory pressure on large `afe` grids if eager loading is retained.
- **Mitigation**: backend abstraction + compile-time feature gates + dual CI lanes.

## 8. Recommendation
Adopt HDF5 backend incrementally with strict optionality in build system until all major CI/compiler targets confirm stability.
