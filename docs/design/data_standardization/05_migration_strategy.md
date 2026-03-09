# Migration Strategy: Legacy `data/` -> FSDS Standard

## 1. Principles
- No physics change during migration.
- Deterministic conversion with checksums and provenance.
- Side-by-side operation until parity is proven.

## 2. Migration phases

## Phase 0 — Inventory and schema freeze
- Catalog all current assets by role:
  - isochrones, base spectra, auxiliary spectra, dust, nebular, photometry/calibration, xrb.
- Freeze FSDS v1 schema and required metadata fields.

## Phase 1 — Build conversion tooling
Create `tools/fsps_data_convert` (Python recommended for rapid HDF5 tooling):
- Inputs: existing text/unformatted binaries
- Outputs: `fsps_data_v1.h5`
- Emit conversion report:
  - source files + checksums
  - dataset shapes/ranges
  - unit/transform declarations
- Generate a **mock 5D spectral dataset** (`lambda,z,afe,logt,logg`) for early backend/mapper/interpolator validation, even before full production C3K-style grids are available.

## Phase 2 — Golden-data generation
- Convert the official default dataset once.
- Store artifact checksum manifest (`sha256`) in repo or release assets.

## Phase 3 — Dual-path validation
For each library tuple (e.g., mist+miles, pdva+basel, mist+c3k, bpss):
1. Run legacy ingestion
2. Run FSDS ingestion
3. Compare key internal arrays and final outputs under existing tolerances

## Phase 4 — User custom-data pipeline
Provide documented recipe:
1. Prepare source arrays + metadata
2. Run `fsps_data_convert --validate`
3. Run `test_runner` with `data_backend=hdf5`

## Phase 5 — Default switch
- Make FSDS backend default once CI and parity are stable.
- Keep legacy adapter behind explicit compatibility switch for at least one major release.

## 3. Required migration tests
- Schema validation tests
- Round-trip integrity checks
- Regression parity on `.spec`, `.mags`, and selected in-memory grids
- Startup-time benchmark vs legacy path
- Rank-handling tests:
  - 4D input -> 5D internal with `afe=1` (degenerate-axis path)
  - native 5D input with `afe>1` (full interpolation path)
  - sparse-grid tests using `missing_value` and flat-model representation

## 4. Operational details
- Version each standardized data artifact (`fsps_data_v1.h5`, `v1.1`, ...).
- Record exact converter version and git commit in root metadata.
- Maintain changelog of dataset updates independent from code release.

## 5. Exit criteria
- 100% required legacy roles represented in FSDS.
- CI parity green across supported compiler matrix.
- Custom dataset onboarding documented and reproducible.
