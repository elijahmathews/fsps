# FSPS Unified Data Standard (FSDS) — Format Specification

## 1. Goal
Define a self-describing, high-performance, versioned data format so FSPS can ingest default and user-provided datasets without hard-coded filename logic in `fsps_io`/`fsps_initialization`.

## 2. Design Decision
### Recommended primary format: **HDF5 container + in-file schema metadata**
- **Why**: native n-D arrays, chunking, compression, attributes, mature tooling, parallel I/O pathways.
- **Fortran fit**: stable HDF5 Fortran API; clean mapping to allocatable/pointer arrays.
- **Performance**: contiguous reads for full-grid initialization; chunked reads for selective metallicity/time slices.
- **Memory safety**: design for hyperslab access so FSPS can read only requested `(z, afe)` slices instead of full 5D cubes.

### Optional sidecar for authoring: `manifest.toml`
- Not required at runtime.
- Used by conversion tooling to generate HDF5 consistently.

## 3. File Layout (single file: `fsps_data_v1.h5`)

```text
/
  attrs:
    fsds_version = "1.0.0"
    producer = "fsps-data-tools"
    produced_at = "2026-02-24T00:00:00Z"
    compatibility_min_fsps = "3.2.0"
    default_library_set = "mist+miles+dl07"

  /registry
    /libraries
      isochrone_default = "mist"
      spectral_default  = "miles"
      dust_default      = "DL07"

  /axes
    /wavelength_angstrom
    /zlegend/{mist,pdva,prsc,bsti,gnva,bpss,miles,basel,c3k_afe+0.0,...}
    /alpha_fe/{miles,basel,c3k_afe+0.0,...}
    /time_logyr/{mist,pdva,prsc,...}
    /logt/{speclib,wmb,wr,pagb,agb_*}
    /logg/{speclib,wmb}
    /nebular/{logz,logu,age_logyr,line_lambda}

  /libraries
    /isochrones/<name>
      /tracks
        mini[z,t,m]
        mact[z,t,m]
        logl[z,t,m]
        logt[z,t,m]
        logg[z,t,m]
        ffco[z,t,m]
        phase[z,t,m]
        lmdot[z,t,m]
        nmass[z,t]
        timestep_logyr[z,t]
      attrs:
        zsol
        source
        native_units

    /spectra/<name>
      /base
        lambda_angstrom[nspec]
        resolution_fwhm_angstrom[nspec]
        spectral_grid_nd[...]
      attrs:
        flux_unit = "Lsun/Hz/Msun"   # canonical FSPS internal convention tag
        storage_dtype = "float32"    # optional if file precision differs
        zsol_spec

    /auxiliary
      /wmbasic
      /wr
      /agb
      /post_agb
      /xrb

    /dust
      /emission/<DL07|THEMIS>
      /attenuation/{wg_h,wg_c,g03_smc}
      /agn/nenkova08
      /dusty_agb/{orich,crich}

    /nebular/<isoc_name>/<WD|ND>
      continuum_logflux[nspec,nebnz,nebna,nebnu]
      lines_logflux[nline,nebnz,nebna,nebnu]
      axes references via attrs

    /photometry
      filters
      vega_sed
      sun_sed
      indices
```

## 4. Dataset metadata contract
Each dataset must define attributes:
- `role` (e.g., `isochrone_track`, `spectral_cube`, `nebular_continuum`)
- `dims` (ordered dim names, e.g., `"lambda,z,afe,logt,logg"` or `"lambda,z,logt,logg"`)
- `unit` (UDUNITS-style text)
- `dtype` (`float32`, `float64`, `int32`)
- `missing_value` (required for non-rectangular physical grids embedded in rectangular arrays; use numeric sentinel, not `NaN`)
- `transform` (e.g., `none`, `log10`, `log10_plus_floor:1e-95`)
- `axis_refs` (HDF5 paths to axes)
- `provenance` (source publication/version/checksum)
- `representation` (`dense_nd` or `flat_models`)

For datasets with missing nodes, one of these must be present:
- `missing_value = -1.0e99` (or similarly unphysical finite value), and/or
- companion validity mask dataset via `valid_mask_ref` (preferred for long-term robustness).

### 4.1 Arbitrary N-dimensional parameter spaces
The schema is axis-driven, not hardcoded by path or fixed rank. Any spectral grid rank is legal if:
1. `dims` and `axis_refs` are complete and ordered,
2. all non-`lambda` axes have corresponding axis datasets,
3. `missing_value` semantics are defined.

### 4.2 Alpha-abundance axis
Add an explicit abundance axis for spectral libraries that support it:
- Axis name: `afe` (alias `alpha_fe`)
- Recommended path: `/axes/alpha_fe/<library>`
- Unit: `dex`

### 4.3 Sparse / non-rectangular grid handling
Because grids such as C3K can be physically sparse across `(z, afe, logt, logg)`, FSDS supports two representations:

1. **Dense ND with sentinel nodes**
  - Dataset shape example: `[lambda, z, afe, logt, logg]`
  - Missing nodes are filled with finite `missing_value` (e.g., `-1.0e99`).
  - Optional companion mask: `valid_mask[z,afe,logt,logg]` (`uint8` 0/1).

2. **Flat model table for highly sparse grids**
  - `flux[lambda, n_models]`
  - `model_parameters[n_models, n_param]` with ordered columns, e.g. `(z, afe, logt, logg)`
  - `model_parameter_names` attribute to declare column mapping.

### 4.4 Required I/O granularity for large spectral grids
For spectral datasets with `representation=dense_nd`, providers must chunk to permit efficient hyperslab reads of:
- `flux[:, z_i, afe_j, :, :]` and
- small neighborhoods around `(z_i, afe_j)`.

For `representation=flat_models`, providers must include index datasets enabling bounded subset selection, e.g.:
- `index_by_z_afe` or equivalent sorted-key lookup metadata.

### 4.5 Compiler-safe missing data policy
FSPS may be compiled with aggressive floating-point optimization (`-ffast-math`, `-Ofast`). Therefore:
- `NaN`-dependent logic is forbidden in FSDS-required behavior.
- Missing-node detection must use explicit finite sentinel comparisons and/or mask arrays.

Backends may store either representation; mappers must expose a common internal form to FSPS.

## 5. Canonical naming and normalization
Replace hard-coded name branching with canonical names in metadata:
- Isochrones: `mist`, `pdva`, `prsc`, `bsti`, `gnva`, `bpss`
- Spectra: `miles`, `basel`, `bpass`, `c3k_*`, `ckc14`
- Dust: `DL07`, `THEMIS`

Aliases (`padova`, `parsec`, `themis`) stored in `/registry/aliases`.

## 6. Precision and storage policy
- Keep physically large cubes as `float32` by default (matches current `.spectra.bin` footprint intent).
- Promote to FSPS working precision (`WP`) during load where needed.
- Keep critical axes and metadata as `float64`.

## 7. Chunking and compression policy
- Spectral cubes: chunk major axis by `lambda` slabs for fast contiguous reads.
- Isochrones: chunk by `(z,t)` tiles.
- Compression: `gzip` level 1–4 (default 2) to balance speed and size.
- Provide a `no-compression` profile for HPC local scratch.

## 8. Versioning and compatibility
- Semantic version in root attr: `fsds_version`.
- Backward-compatible additions allowed as new groups/attrs.
- Breaking changes require major bump and explicit reader adapter.

## 9. Validation rules
- All required groups/datasets exist for selected model tuple.
- Axes monotonic where expected (`lambda`, `time`).
- Shape checks against declared `dims`.
- Unit/transform compatibility checks before ingestion.
- Missing-node policy checks (`missing_value` finite and/or mask present and shape-compatible).
- Checksum verification (optional strict mode).

## 10. Why not NetCDF-only or JSON+raw-only?
- **NetCDF4** is viable but less flexible for mixed hierarchical model families in one file.
- **JSON/TOML + raw binaries** keeps dependencies lighter but reintroduces brittle synchronization and weak atomicity; acceptable only as transitional interchange, not canonical runtime store.

## 11. Transitional support
FSPS should support:
1. `legacy` loaders (current paths/files)
2. `fsds_hdf5` loader (new standard)

Runtime selection via config/environment, with `fsds_hdf5` preferred when available.
