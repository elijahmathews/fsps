# Modern Fortran Architecture for Self-Describing FSPS Data

## 1. Architectural intent
Refactor ingestion into a metadata-driven pipeline with a backend abstraction. `fsps_initialization` should request **roles** (isochrones, spectra, nebular, etc.), not hard-coded files.

## 1.1 Practical Fortran constraint
Fortran requires array rank to be known at compile time. A runtime-polymorphic `real, allocatable :: flux(::)` style container is not suitable for high-performance physics kernels expecting fixed-rank arrays.

## 1.2 Practical memory constraint
Large `afe`-expanded libraries cannot be loaded as full 5D arrays on all targets. Architecture must support lazy loading/hyperslab reads as a first-class feature.

## 2. Proposed module split
- `src/io/fsps_data_schema.f90`
  - Data descriptor types (`dataset_desc_t`, `axis_desc_t`, `library_manifest_t`)
- `src/io/fsps_data_backend.f90`
  - Abstract backend interface (`data_backend_t`)
- `src/io/fsps_data_backend_hdf5.f90`
  - HDF5 implementation
- `src/io/fsps_data_backend_legacy.f90`
  - Adapter to existing text/binary loaders
- `src/io/fsps_data_mapper.f90`
  - Maps descriptors to `fsps_context_state_t` arrays
- `src/io/fsps_data_registry.f90`
  - Library aliases, compatibility checks, selection logic

## 3. Core derived types

```fortran
type :: axis_desc_t
  character(len=:), allocatable :: name
  character(len=:), allocatable :: unit
  character(len=:), allocatable :: path
  integer :: n = 0
  real(WP), allocatable :: values(:)
end type

type :: dataset_desc_t
  character(len=:), allocatable :: role
  character(len=:), allocatable :: path
  character(len=:), allocatable :: dims_csv
  character(len=:), allocatable :: unit
  character(len=:), allocatable :: transform
  character(len=:), allocatable :: dtype
  integer, allocatable :: shape(:)
  logical :: required = .true.
end type

type :: library_manifest_t
  character(len=:), allocatable :: isoc_name
  character(len=:), allocatable :: spec_name
  character(len=:), allocatable :: dust_name
  type(dataset_desc_t), allocatable :: datasets(:)
  type(axis_desc_t), allocatable :: axes(:)
end type

type :: spectral_grid_t
  ! Canonical shape for dense in-memory materialization (optional).
  ! [lambda, z, afe, logt, logg], with degenerate afe (=1) for 4D libraries.
  real(WP), allocatable :: flux(:,:,:,:,:)
  ! Optional validity mask for sentinel-safe missing-node handling.
  logical, allocatable :: valid(:,:,:,:)
  real(WP), allocatable :: axis_lambda(:)
  real(WP), allocatable :: axis_z(:)
  real(WP), allocatable :: axis_afe(:)
  real(WP), allocatable :: axis_logt(:)
  real(WP), allocatable :: axis_logg(:)
end type

type :: spectral_slice_t
  ! Runtime working set for lazy mode: one (z,afe) slab.
  ! Shape: [lambda, logt, logg]
  real(WP), allocatable :: flux(:,:,:)
  logical, allocatable :: valid(:,:)
  integer :: iz = 0
  integer :: iafe = 0
end type

type :: sparse_model_table_t
  ! For representation=flat_models kept sparse in memory.
  ! flux shape: [lambda, n_models]
  real(WP), allocatable :: flux(:,:)
  real(WP), allocatable :: params(:,:)    ! [n_models, n_param], e.g. z,afe,logt,logg
  logical, allocatable :: valid_model(:)
end type
```

## 4. Backend abstraction

```fortran
type, abstract :: data_backend_t
contains
  procedure(open_backend_if), deferred :: open
  procedure(close_backend_if), deferred :: close
  procedure(read_manifest_if), deferred :: read_manifest
  procedure(read_real_1d_if), deferred :: read_real_1d
  procedure(read_real_2d_if), deferred :: read_real_2d
  procedure(read_real_4d_if), deferred :: read_real_4d
  procedure(read_real_5d_if), deferred :: read_real_5d
  procedure(query_axis_if), deferred :: query_axis
  procedure(find_bracketing_if), deferred :: find_bracketing_indices
  procedure(read_spectral_slice_if), deferred :: read_spectral_slice
  procedure(read_spectral_neighborhood_if), deferred :: read_spectral_neighborhood
  procedure(read_sparse_subset_if), deferred :: read_sparse_subset
  procedure(read_int_nd_if),  deferred :: read_int_nd
  procedure(has_path_if),     deferred :: has_path
end type
```

This keeps `fsps_initialization` backend-agnostic while preserving rank-specific Fortran performance.

Interface intent:
- `read_spectral_slice`: load only `flux[:, iz, iafe, :, :]` into `spectral_slice_t`.
- `read_spectral_neighborhood`: load minimal dense block around bracketing `(iz, iafe)` (typically 2x2 in z/afe).
- `read_sparse_subset`: for flat-model data, return only models in bounded parameter neighborhood.

## 5. Mapping strategy to `fsps_context_state_t`
The mapper owns all conversions:
- dtype promotion (`float32 -> WP`)
- transforms (`log10_plus_floor`, unit normalization)
- axis compatibility checks
- shape assertions against dynamic dimensions

### 5.1 Degenerate-dimension strategy for alpha axis
- Internal canonical spectral storage is rank-5: `[lambda, z, afe, logt, logg]`.
- If input data is 4D (`lambda,z,logt,logg`), mapper allocates `afe=1` and writes to `flux(:, :, 1, :, :)`.
- Interpolation code branches on `size(flux, 3)`:
  - `==1`: bypass alpha interpolation
  - `>1`: perform multilinear interpolation over `(z, afe, logt, logg)`

This avoids flattening/index arithmetic in physics kernels.

### 5.2 Sparse-grid strategy (selected)
Chosen strategy: **bounded dense neighborhood materialization**.

- On disk, sparse libraries may remain `flat_models`.
- In memory, mapper/backend assemble only a local dense neighborhood around requested `(z, afe)` into `spectral_slice_t` or a tiny `(z,afe)` block.
- Interpolation in `logt/logg` then uses existing structured-grid kernels on this local dense block.
- If neighborhood cannot be densified (insufficient nodes), fallback is explicit nearest-neighbor in parameter space with warning/quality flag.

Rationale: preserves fast structured interpolation paths while avoiding full dense 5D expansion.

Existing state arrays can remain in place initially; loader behavior changes first.

## 6. Runtime flow (OO)
1. `factory_create_backend(mode)` (`legacy`/`hdf5`)
2. `backend%open(path)`
3. `backend%read_manifest(manifest)`
4. `registry_select_libraries(manifest, requested_names)`
5. `mapper_prepare_axes_and_indexing(ctx, backend, manifest)`
6. `mapper_load_required_slices(ctx, backend, requested_z, requested_afe)`
7. `backend%close()`

## 7. Error model
Introduce structured errors (instead of broad `error stop`):
- `FSPS_IO_ERR_NOT_FOUND`
- `FSPS_IO_ERR_SCHEMA`
- `FSPS_IO_ERR_DIM_MISMATCH`
- `FSPS_IO_ERR_UNIT_MISMATCH`
- `FSPS_IO_ERR_BACKEND`

Return status + message; optionally escalate to `error stop` in strict mode.

## 8. Caching implications
Current setup cache remains valid. Add cache key components:
- backend mode
- data file path
- schema version
- selected library IDs
- slice policy (`eager_full`, `lazy_slice`)
- requested `(z, afe)` indices or neighborhood signature

This prevents stale cache collisions across formats.

## 9. Threading and safety
- Backend handles are context-local.
- Shared cache remains read-mostly.
- If future parallel init is introduced, backend reads should be serialized per handle or use independent handles.

## 10. Minimal integration sequence
1. Add backend abstraction + legacy adapter (no behavior change).
2. Move current `get_*filename` routing into registry/adapter.
3. Add HDF5 backend with `query_axis` + `read_spectral_slice` first.
4. Introduce lazy slice cache in `fsps_context_state_t`.
5. Add sparse `flat_models` subset path and neighborhood materializer.
6. Expand to isochrones, dust, nebular, auxiliary libraries.

## 11. Interpolation-module impact
`src/math/fsps_interpolation.f90` and callers must be extended for optional alpha-axis interpolation. The new interpolation path should remain branch-free for the common degenerate case (`afe=1`) as much as possible.

Additional impacts:
- add neighborhood interpolation helper for `(z,afe)` bilinear weights over slice blocks,
- add explicit nearest-neighbor fallback for sparse subset failure,
- propagate interpolation quality flags upward for diagnostics/tests.
