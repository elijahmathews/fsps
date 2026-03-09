# Refactoring Ingestion Workflow (`fsps_io` + `fsps_initialization`)

## 1. Current pain points observed
- Library and file selection is hard-coded with `select case` + string conventions.
- Logic is distributed across many routines with implicit assumptions about dimensions and units.
- Data interpretation (e.g., transforms, floors, interpolation expectations) is embedded in procedural code, not metadata.

## 2. Target workflow (metadata-driven)

```text
fsps_initialize_data
  -> resolve requested model tuple
  -> create backend (legacy|hdf5)
  -> read and validate manifest/schema
  -> query axes and locate bracketing indices
  -> allocate dimensions from manifest + options
  -> lazily load required stellar hyperslabs
  -> load interstellar/nebular/dust
  -> load photometry/calibration
  -> compute derived quantities
  -> cache publish
```

## 3. New responsibilities

### `fsps_initialization`
- Keep orchestration only.
- Stop constructing filenames.
- Decide slice policy (`eager_full` for small libs, `lazy_slice` for large libs).
- Resolve requested `(z, afe)` coordinates to index bounds before spectral load.
- Request roles from loader API:
  - `load_role("isochrones")`
  - `load_role("spectral_base")`
  - `load_role("dust_emission")`
  - etc.

### `fsps_io`
- Become facade over backend + mapper.
- Provide role-centric public API:
  - `fsps_data_open()`
  - `fsps_data_load_stellar()`
  - `fsps_data_load_interstellar()`
  - `fsps_data_load_photometry()`
  - `fsps_data_close()`

## 4. Refactor map from existing routines
- Keep existing loaders as **legacy adapter** internals:
  - `read_isochrone_database`
  - `read_spectral_binary`
  - `load_nebular_grid`
  - `load_dust_emission_table`
  - etc.
- Add a thin dispatcher deciding backend path.

## 5. Dynamic dimensioning strategy
Today dimensions are hard-coded by library names (`set_simulation_dimensions`).

Target:
1. Read dimensions from schema for selected library tuple.
2. Validate against compile-time maxima where required.
3. Populate `ctx%state%nt/nz/nspec/nzinit` from data (except where physics flags explicitly override).

## 5.1 Hyperslab and lazy-load requirements
For large spectral libraries, ingestion must avoid full 5D reads.

Required control flow:
1. Query axis vectors (`z`, `afe`, `logt`, `logg`).
2. Find bracketing indices for requested `(z*, afe*)`.
3. Read only required slab(s):
  - single-slice path: `flux[:, iz, iafe, :, :]`
  - interpolation neighborhood path: `flux[:, iz:iz+1, iafe:iafe+1, :, :]`
4. Cache slices/neighborhoods for reuse in SSP loops.

## 6. Mapper shape-shifting responsibilities
`fsps_data_mapper.f90` is the authoritative reshaping layer between backend payloads and fixed-rank internal arrays.

Required behaviors:
1. Read `dims`/`dims_csv` (e.g., `"lambda,z,logt,logg"` or `"lambda,z,afe,logt,logg"`).
2. Populate canonical internal spectral storage `[lambda,z,afe,logt,logg]`.
3. Insert degenerate `afe=1` dimension automatically when absent.
4. Apply `missing_value` masking rules and transform policy.
5. Support both dense-ND and flat-model representations (when configured).
6. Materialize bounded dense neighborhoods from sparse `flat_models` subsets when needed.

## 7. Interpolation/transforms policy
Centralize all transform declarations in metadata:
- `none`
- `log10`
- `log10_plus_floor:x`
- `angstrom_to_hz` (derived)

Interpolation remains in Fortran compute layer but is selected by metadata + policy flags.

Interpolation upgrade requirement:
- Extend `src/math/fsps_interpolation.f90` call paths to support multilinear interpolation over `(z, afe, logt, logg)` when `size(speclib_flux,3) > 1`.
- Preserve fast legacy path for degenerate alpha dimension (`afe=1`).
- Sparse strategy: structured interpolation is preferred; if sparse neighborhoods cannot be densified, use explicit nearest-neighbor fallback in parameter space and emit quality flags.

## 8. Backward compatibility contract
- Default mode initially: `legacy`.
- Opt-in mode: `fsds_hdf5`.
- Dual-mode regression tests must produce equivalent outputs within existing tolerance (`FSPS_TEST_RTOL`).

## 9. Suggested API sketch

```fortran
call fsps_data_open(ctx, source_uri, backend_mode)
call fsps_data_select(ctx, isoc=..., spec=..., dust=...)
call fsps_data_load_stellar(ctx, zin)
call fsps_data_load_interstellar(ctx)
call fsps_data_load_photometry(ctx)
call fsps_data_close(ctx)
```

## 10. Rollout phases
- **Phase A**: backend abstraction + legacy adapter.
- **Phase B**: HDF5 read path for spectral core with hyperslab reads only.
- **Phase C**: sparse `flat_models` subset + neighborhood materializer.
- **Phase D**: all roles migrated.
- **Phase E**: deprecate direct filename helpers (retain in compatibility module).

## 11. Definition of done
- No direct file path logic in orchestration code.
- New dataset can be loaded by metadata-only registration.
- Large `afe` libraries initialize without full-cube resident memory requirement.
- Full unit + regression suite green for both backend modes.
