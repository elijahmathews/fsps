# Data libraries: survey and standard format

This document describes how the code reads its input libraries (stellar spectra, isochrones,
nebular grids, dust, filters, ...), the standard format that replaces the legacy layout, and the
tools that convert, package and fetch the data.

Contents:

1. [Survey: how the libraries are read today](#1-survey-how-the-libraries-are-read-today)
2. [The standard format (version 1)](#2-the-standard-format-version-1)
3. [Library kinds](#3-library-kinds)
4. [Tools](#4-tools)
5. [How the code selects a format](#5-how-the-code-selects-a-format)
6. [Test library, data releases and CI](#6-test-library-data-releases-and-ci)

## 1. Survey: how the libraries are read today

### Locating the data

`fsps_resolve_paths` (`src/core/fsps_environment.f90`) picks the data root `ctx%sps_home` from,
in order: `SPS_HOME`, `FSPS_DATA_HOME`, `$XDG_DATA_HOME/fsps`, `~/.local/share/fsps`,
`/usr/share/fsps`, `/usr/local/share/fsps`. A candidate is accepted if
`<root>/data/allfilters.dat` or `<root>/data/FILTER_LIST` exists. Every reader builds paths as
`trim(ctx%sps_home) // '/data/...'` with the file names written into the code. `SPS_HOME` wins
over `FSPS_DATA_HOME` when both are set.

### Selecting a library

Libraries are chosen at run time by name, through `fsps_setup(ctx, zin, isoc, spec, dust)`
(defaults `mist`, `miles`, `DL07`); nothing is selected at build time. The names are normalized in
`normalize_library_names` (`src/core/fsps_initialization.f90`): isochrones `mist`, `pdva`
(Padova 2007), `prsc` (PARSEC), `bsti` (BaSTI), `gnva` (Geneva), `bpss` (BPASS, which forces the
spectral library `bpass`); spectral libraries `miles`, `basel`, `c3k_*` (e.g. `c3k_afe+0.0`);
dust emission `DL07` or `THEMIS`. Everything else (nebular grid, hot-star and AGB spectra,
filters, ...) follows from these names or is fixed.

All grid sizes are compile-time constants or `select case` branches on the names
(`set_simulation_dimensions`, `configure_dust_model`, `src/core/fsps_constants.f90`), so a
library with another size or name cannot be used without code changes:

| Library | Sizes set in code | Other metadata set in code |
|---|---|---|
| Isochrones | `nt`, `nz` per name (`NT_MIST` = 107, `NZ_MIST` = 12, ...); `NM` = 2000 points per isochrone | solar metallicity per name (MIST 0.0142, Padova 0.019, PARSEC 0.01524, BaSTI/Geneva/BPASS 0.020) |
| Spectral library | `nspec`, `nzinit` per name (MILES 5994 x 5, BaSeL 1963 x 6, C3K 11149 x 11) | solar metallicity per name (MILES 0.019, BaSeL 0.020, C3K 0.0134); the (log T, log g) grid is always read from the BaSeL files (`basel_logt.dat`, `basel_logg.dat`), `NDIM_LOGT` = 68, `NDIM_LOGG` = 19 |
| Nebular | `NEBNZ` = 11, `NEBNAGE` = 10, `NEBNIP` = 7, `NLAM_NEBCONT` = 1963, `NEMLINE` = 166 | available only for `mist`, `pdva`, `prsc`, `bpss` |
| Dust emission | `ndim_dustem`, `numin`, `nqpah` per name | the q_PAH and U_min grids are arrays in the code, not in the files |
| Filters | `nbands` = 159 | the IR/MIPS calibration (indices and reference wavelengths) is hardcoded by filter index |
| X-ray binaries | `NT_XRB` = 10, `NZ_XRB` = 11, `NSPEC_XRB` = 15000 | ages and [Fe/H] grid in the code |
| Hot stars, AGB, post-AGB, dusty AGB, AGN, attenuation | all sizes are constants in `fsps_io` / `fsps_constants` | WMBasic log g grid (3.5, 4.0, 4.5) in the code |

### The files, their layouts and the implicit metadata

All readers are in `src/io/fsps_io.f90` except `load_dusty_agb_spectra`,
`load_xray_nebular_grid`, `load_xrb_spectra` and the (log T, log g) grid
(`src/core/fsps_initialization.f90`) and the user IMF table (`src/physics/fsps_ssp.f90`). Text
files are read list-directed into FP64 (`WP`). "Processing" is what the reader does after
reading; it stays in the code in the new design (Section 5).

| Library (reader) | Files under `data/` | Layout | Processing on load |
|---|---|---|---|
| Stellar spectra (`read_spectral_binary`, `load_wavelength_grid`, `load_spectral_resolution`, `load_zlegend_file`) | `spectra/MILES/imiles_z<Z>.spectra.bin`, `spectra/BaSeL3.1/basel_wlbc_z<Z>.spectra.bin`, `spectra/C3K/<name>_z<Z>.spectra.bin`; `*.lambda`, `*.res`, `zlegend.dat`; grid `BaSeL3.1/basel_log{t,g}.dat` | one file per metallicity: a Fortran direct-access record of FP32 values, Fortran order `(nspec, 68 logT, 19 logg)`, no header, native byte order; `<Z>` is `F6.4` of the zlegend value | promoted to FP64; interpolated in log flux (floor 1e-70) onto the isochrone metallicities |
| Isochrones (`read_isochrone_database`) | `isochrones/<set>/isoc_z<Z>.dat`, `zlegend.dat` (MIST: `isoc_z{m,p}X.XX.dat`, zlegend holds `m1.50`-style strings converted with Z = Zsun 10^[M/H]) | text; isochrones separated by `#` lines (or any non-numeric line); columns log age, M_ini, M_act, log L, log T, log g, composition (C/O or surface X), phase, and for MIST only log Mdot | log Mdot = -99 when absent; the age of an isochrone is the age on its first row; at most `NM` rows per isochrone |
| BPASS (`read_bpass_data`) | `isochrones/BPASS/bpass.lambda`, `bpass.mass`, `bpass_v2.2_salpeter100.ssp.bin`, `zlegend.dat` | SSP spectra, FP64 direct access, Fortran order `(15000, 43 ages, 12 Z)`; `bpass.mass` rows: log age and the stellar mass at 12 Z | none |
| Nebular (`load_nebular_grid`, `load_xray_nebular_grid`) | `nebular/ZAU_{ND,WD}_<isoc>.{cont,lines}`, `ZAU_WX_{ND,WD}_bpss.*` (`WD`/`ND`: with/without dust in the HII region, `cloudy_dust`; `WX`: X-ray-binary-irradiated) | text; header line, wavelength (or line wavelength) line, then for each (log Z, age, log U) a label line and a spectrum line, loops Z outermost and U innermost | continuum: log10(f + 1e-95) interpolated linearly onto the spectral grid; lines: log10(L + 1e-95); age converted to log10; the line labels must equal the continuum labels |
| Dust emission (`load_dust_emission_table`) | `dust/dustem/{DL07,THEMIS}_MW3.1_<q>.dat`, one file per q_PAH | 2 header lines; rows: wavelength (um) and two spectra per U_min (delta-function U_min, and U_min..U_max power law) | wavelength x 1e4; linear interpolation onto the grid for lambda >= 1 um |
| AGN torus (`load_agn_dust_models`) | `dust/Nenkova08_y010_torusg_n10_q2.0.dat` | 3 header lines, a line of 9 optical depths, then 125 rows: wavelength (A), 9 spectra (f_nu, normalized to unity) | log-log interpolation inside the model's wavelength range, zero outside |
| Filters (`load_filter_definitions`, `apply_legacy_filter_norm`) | `allfilters.dat` (or `alt_filter_file`) | one header line, then per filter a `#` comment line and (wavelength, transmission) rows | repeated wavelengths dropped, negative transmissions clamped; linear interpolation onto the grid; normalized to integral T/lambda = 1; filters 53-56, 95-104 renormalized for a flat nu f_nu and 90-92 for a beta = 2 blackbody (only for the default list) |
| Reference spectra (`load_standard_sed`) | `spectra/A0V_KURUCZ_92.SED` (Vega, F_lambda, 1 header line), `spectra/SUN_STScI.SED` (Sun, f_nu at 10 pc, no header); 1221 rows each | (wavelength, flux) rows | log-log interpolation; Vega multiplied by lambda^2; zero (1e-308) redward of the table |
| Spectral indices (`load_index_definitions`) | `allindices.dat` | 4 header lines, 30 rows of 7 numbers | first 25 rows converted from air to vacuum |
| Attenuation (`load_attenuation_curves`) | `dust/alldirty_{h,c}.dat`, `dust/Gordon03_table4.dat` | alldirty: 3 header lines, 18 blocks (tau_V) of 25 rows: wavelength, tau_V, 6 curves; SMC: 30 rows, wavelength (um) descending | linear interpolation; constant blueward, zero redward |
| WMBasic O stars (`load_wmbasic_spectra`) | `spectra/Hot_spectra/WMBASIC.teff`, `WMBASIC_zlegend.dat`, `WMBASIC_z<Z>.spectra` (12 Z) | rows: wavelength and 3 x 11 spectra (log g outer, T inner); log g grid in the code | onto the grid, then log-flux interpolation onto the isochrone metallicities |
| Wolf-Rayet (`load_wr_spectra`) | `CMFGEN_{WN,WC}.teff`, `CMFGEN_{WN,WC}_Zall.spec` | wavelength line, then 5 Z x 12 T blocks: a (T, Z) label line and a spectrum line | log-flux interpolation in Z, log-log in wavelength |
| Post-AGB (`load_post_agb_spectra`) | `ipagb.teff` (Teff in K), `ipagb_{solar,halo}.spec` | rows: wavelength and 14 spectra | linear onto the grid, floor 1e-308 |
| AGB (`load_agb_spectra`) | `AGB_spectra/{Orich,Crich,Crich_Aringer}.{teff,spec}` | O-rich Teff tabulated per spectral type and 22 [Z/Zsun]; C-rich and `Crich_Aringer` one Teff per spectrum; spectra as rows (wavelength, one column per type) | Teff interpolated to the isochrone Z and logged; spectra linear onto the grid, floor 1e-308 |
| Dusty AGB shells (`load_dusty_agb_spectra`) | `dust/dusty/{Orich,Crich}_dusty.spec` | number of wavelengths, wavelength line, then 6 Teff x 50 tau blocks: (Teff, tau) label and a transfer-function line | linear onto the grid from the first tabulated wavelength |
| X-ray binaries (`load_xrb_spectra`) | `spectra/xrb/xsp.lambda`, `xsp_feh<feh>.spec` (11 [Fe/H]) | 10 lines (ages) of 15000 values each; ages and [Fe/H] in the code | linear onto the grid, floor 1e-308, x L_sun |
| LSF, IMF, SFH (`load_lsf_data`, `fsps_ssp`, `load_tabular_sfh`) | `lsf.dat`, `imf.dat` (or `imf_filename`), `sfh.dat` (or `sfh_filename`) | small user-editable text tables | unchanged by this work |

Unit conventions that are implicit today (they become explicit `units` in the manifests):

- Main stellar library: Eddington flux H_nu (erg s^-1 cm^-2 Hz^-1); the code multiplies by
  16 pi^2 R^2 / L_sun, with R from (M, L, Teff).
- WMBasic, CMFGEN, post-AGB, `Crich_Aringer`: L_nu per unit bolometric luminosity, scaled by L_bol.
- Empirical O- and C-rich AGB (`Orich.spec`, `Crich.spec`): F_lambda normalized to unit integral over Angstroms,
  converted to f_nu by lambda^2 / c in the code.
- Nebular continuum and lines: per unit ionizing photon rate, multiplied by Q_H in the code.
- X-ray binaries: per L_sun (multiplied by L_sun on load).

### Observations that shaped the design

- The large binaries (spectral libraries, BPASS) are already plain arrays: a direct-access
  record has no record markers, so `imiles_z0.0190.spectra.bin` is exactly 5994 x 68 x 19 FP32
  values. Only the metadata (sizes, axes, solar metallicity, byte order) is missing.
- The text libraries carry most of their metadata implicitly: in file names (`F6.4` of Z), in
  `zlegend.dat`, in comment lines, and in code (grids for q_PAH, U_min, WMBasic log g, XRB ages).
  Values are parsed once from text at each setup (the MIST isochrones are 140 MB of text, the
  slowest part of a cold setup).
- Every reader does two things: get the raw table, then process it onto the model grids. Only
  the first depends on the file layout. The new loader therefore replaces only that step: both
  layouts produce the same raw arrays, and the processing code is shared, so results cannot
  differ between the two layouts.
- Some physics branches on the isochrone name (`mist`: WR temperatures from Mdot; `bsti`: log g
  from (M, L, T); `pdva`: horizontal-branch detection; `gnva`: lower IMF bound;
  `bpss`: SSP spectra instead of isochrones). Those stay keyed on the library name; a new
  isochrone library behaves like a generic one (Section 3.2).

## 2. The standard format (version 1)

### Goals

- Self-describing: every size, axis, unit and normalization is in the files, so a new library
  needs no code change as long as it fits one of the library kinds of Section 3.
- Simple to produce: a JSON manifest and raw little-endian arrays. A library can be written from
  Python with the standard `json` module and `numpy.ndarray.tofile`, without reading Fortran.
- Fast and exact to read: arrays are read with Fortran stream I/O at the stored precision (no
  text parsing), and converted libraries reproduce the legacy ones bit for bit.
- Extensible: arrays have named axes, and a library may have axes the code does not use yet
  (e.g. [alpha/Fe], IMF slope); the loader selects one value along each extra axis.

### Why JSON

The Fortran side has to parse the manifest without external dependencies. JSON's grammar is
small (six value types, no indentation rules, one string syntax), so a complete parser is a few
hundred lines of Fortran (`src/io/fsps_json.f90`), and Python reads and writes it with the
standard library. TOML is pleasanter to edit by hand, but its grammar (tables, arrays of tables,
dotted keys, four string forms, dates) needs a much larger parser, and Python's standard library
can only read it. Manifests are written by tools, so hand-editing matters less than a parser that
is small and complete. JSON has no comments; manifests have `title` and `notes` fields instead.

### Layout of a data root

```
<root>/
  fsps-data.json                      root index (format version, release, list of libraries)
  libraries/<kind>/<name>/manifest.json
  libraries/<kind>/<name>/*.bin       arrays referenced by the manifest
  data/                               small text tables in their legacy formats (allindices.dat,
                                      lsf.dat, imf.dat, sfh.dat, ...) and user input files
```

`<kind>` is one of the kinds of Section 3, `<name>` the library name passed to `fsps_setup`
(e.g. `libraries/stellar_spectra/miles`, `libraries/isochrones/mist`). The small text tables and
user files (`sfh_filename`, `imf_filename`, `alt_filter_file`) keep their legacy formats and
their legacy location under `data/`. A data root is self-contained: it holds everything a setup
reads.

### Root index: `fsps-data.json`

```json
{
  "format": "fsps-data",
  "format_version": 1,
  "release": "2026.10.0",
  "title": "Default libraries converted from the legacy layout",
  "libraries": [
    {"kind": "stellar_spectra", "name": "miles", "path": "libraries/stellar_spectra/miles"},
    {"kind": "isochrones", "name": "mist", "path": "libraries/isochrones/mist"}
  ]
}
```

The loader requires `format` = `"fsps-data"` and `format_version` = 1, and finds libraries by
their conventional path; `libraries` is for the tools (packing, verification).

### Library manifest: `manifest.json`

```json
{
  "format": "fsps-data",
  "format_version": 1,
  "kind": "stellar_spectra",
  "name": "miles",
  "title": "MILES empirical library merged with BaSeL",
  "notes": "free text",
  "provenance": {
    "converted_from": "legacy FSPS layout: data/spectra/MILES",
    "converter": "tools/fsps_data/convert_legacy.py",
    "created": "2026-10-19T12:00:00Z"
  },
  "metadata": {"solar_metallicity": 0.019},
  "axes": {
    "Z":          {"units": "mass fraction", "values": [0.0008, 0.0035, 0.0097, 0.019, 0.03]},
    "logT":       {"units": "log10(K)", "values": [3.30103, "..."]},
    "logg":       {"units": "log10(cm s^-2)", "values": [-1.0, "..."]},
    "wavelength": {"units": "angstrom", "medium": "vacuum",
                   "file": "wavelength.bin", "dtype": "<f8", "shape": [5994],
                   "sha256": "..."}
  },
  "arrays": {
    "flux":       {"file": "flux.bin", "dtype": "<f4", "shape": [5, 19, 68, 5994],
                   "axes": ["Z", "logg", "logT", "wavelength"],
                   "units": "erg s^-1 cm^-2 Hz^-1 (Eddington flux H_nu)", "sha256": "..."},
    "resolution": {"file": "resolution.bin", "dtype": "<f8", "shape": [5994],
                   "axes": ["wavelength"], "units": "angstrom (FWHM)", "sha256": "..."}
  }
}
```

Fields:

- `format`, `format_version`, `kind`, `name` (required): `name` equals the directory name.
- `title`, `notes`, `provenance` (optional, free form): who made the library, from what.
- `metadata` (object): scalar properties a kind requires (Section 3), plus anything else.
- `axes` (object): every axis an array refers to. An axis has `units`, and its coordinate
  values either inline (`values`, a JSON array of numbers) or in a file (`file`, `dtype`,
  `shape`, `sha256`, like an array). An axis whose coordinates are not numbers may have
  `labels` (strings) instead, e.g. `["halo", "solar"]`. Optional: `description`.
- `arrays` (object): the data. Each array has
  - `file`: file name relative to the manifest directory;
  - `dtype`: `"<f4"`, `"<f8"`, `"<i4"` or `"<i8"` (little-endian FP32, FP64, int32, int64;
    the NumPy spelling);
  - `shape`: the dimensions, **slowest-varying first** (C order, NumPy's default);
  - `axes`: the axis name of each dimension, in the same order (omitted for 1-D tables whose
    length is set by another array, e.g. isochrone columns);
  - `units` (string), `description` (optional);
  - `sha256`: hex SHA-256 of the file.

A file holds exactly one array: `product(shape) * itemsize` bytes, no header, no padding. With
NumPy, `a.astype("<f4").tofile(path)` writes it and
`np.fromfile(path, dtype="<f4").reshape(shape)` reads it. The Fortran loader sees the same
bytes as an array with the dimensions reversed (Fortran order), so the C-order shape
`[5, 19, 68, 5994]` is the Fortran array `flux(5994, 68, 19, 5)`.

### Conventions

- Wavelengths in Angstrom, increasing; vacuum unless the axis says `"medium": "air"`.
  Within a library, arrays on the same grid share one axis.
- The fastest-varying (last) dimension of a spectral array is wavelength. Each kind fixes the
  order of its required axes (Section 3).
- Extra axes: an array may have more dimensions than its kind requires, all **before** (slower
  than) the required ones, e.g. `["afe", "Z", "logg", "logT", "wavelength"]`. The loader takes
  one index along each extra axis (Section 5); an extra axis of length 1 needs no selection. So
  a grid in [alpha/Fe] or IMF slope can be stored as one library without a format change, and
  the code can learn to interpolate in it later.
- Numbers in manifests are written so that they round-trip exactly (Python's `json` does this).
  Integers that are sizes are JSON integers.
- Units are free-form strings for people; the loader checks units only where the kind gives a
  choice (e.g. `flux_density` of a reference spectrum).

### Checksums and versions

Each file's SHA-256 is in the manifest. The loader checks that every file has exactly the
expected size (catching truncated downloads at no cost) but does not hash at run time; the
`verify` tool checks the hashes, and `fetch` checks them on download (Section 4).

`format_version` changes only for an incompatible change of this specification; adding a kind,
an optional field or an array does not change it. A loader refuses a version it does not know.
`release` (root index) names a set of libraries, i.e. a data release; it changes whenever any
library changes.

## 3. Library kinds

Each kind fixes the names and order of the axes its arrays must end with (any extra axes come
first, Section 3.2), the arrays it needs and the metadata it reads. "Free" sizes may be anything;
"fixed" sizes are set by arrays of the code and are checked on load. Library names are those
passed to `fsps_setup` (lower case) or, for the kinds the code loads under a fixed name, the name
shown. Units are those of the converted libraries; a new library must use the same.

| Kind / name | Arrays (C order) | Axes, metadata | Sizes |
|---|---|---|---|
| `stellar_spectra/<spec>` | `flux` [Z, logg, logT, wavelength] (FP32 or FP64), H_nu in erg s^-1 cm^-2 Hz^-1; `resolution` [wavelength] (FWHM, the `.res` convention) | metallicity axis `Z` (mass fraction) or `feh` ([Fe/H], Z = Z_sun 10^[Fe/H]); `logT` (log K), `logg` (log cgs); `wavelength` (A); metadata `solar_metallicity` | free: wavelength, metallicities; fixed: 68 logT x 19 logg |
| `isochrones/<isoc>` | `npoints` [Z, age_index] int32; `log_age` [Z, age_index] (log yr); columns `mini`, `mact`, `logl`, `logt`, `logg`, `composition`, `phase`, optional `log_mdot` (else -99), one value per point, points ordered by (extra axes,) Z, isochrone, mass | `Z` or `feh`; `age_index` (1..n); `solar_metallicity` | free: metallicities, ages, points (<= 2000 per isochrone) |
| `ssp/bpss` (isochrone name `bpss`) | `spectra` [Z, log_age, wavelength] (L_sun/Hz per M_sun formed); `stellar_mass` [Z, log_age] | `Z`, `log_age`, `wavelength`; `solar_metallicity` | free |
| `nebular/<isoc>_<dust\|nodust>`, `nebular/bpss_xray_<dust\|nodust>` | `continuum` [logZ, age, logU, wavelength] (L_sun/Hz per ionizing photon s^-1); `lines` [logZ, age, logU, line_wavelength] (L_sun per photon s^-1) | `logZ` (log Z/Z_sun), `age` (yr), `logU`, `wavelength`, `line_wavelength` (line names: `data/emlines_info.dat`) | free: continuum wavelength; fixed: 11 x 10 x 7 grid, 166 lines |
| `dust_emission/<DL07\|THEMIS>` | `spectra` [qpah, umin, component, wavelength] (L_nu per H nucleon) | `qpah` (percent), `umin`, `component` = (delta U_min, U_min..U_max) | free (wavelengths >= 1 um are used) |
| `agn_torus/clumpy` | `spectra` [tau, wavelength] (f_nu, unit normalization) | `tau`, `wavelength` | free wavelength; fixed: 9 tau |
| `filters/allfilters` | `npoints` [filter] int32; `wavelength`, `transmission` (concatenated curves) | `filter` (labels); metadata `filters`: one entry per filter with `name`, `description` and optionally `calibration` {`lambda_ref`, `exponent`} | free: number of filters, points |
| `reference_spectra/<vega\|sun>` | `flux` [wavelength] | metadata `flux_density`: `f_lambda` (multiplied by lambda^2 on load) or `f_nu` | free |
| `attenuation/dirty` | `attenuation` [structure, tau_v, model, wavelength] | `structure` (homogeneous, clumpy), `tau_v`, `model` (labels) | fixed: 2 x 18 x 6 x 25 |
| `attenuation/smc_bar` | `extinction` [wavelength] (A_lambda/A_V) | `wavelength` (A, increasing) | free |
| `hot_star_spectra/wmbasic` | `flux` [Z, logg, logT, wavelength] (L_nu / L_bol) | `Z`, `logg`, `logT`, `wavelength` | free wavelength; fixed: 12 Z x 3 logg x 11 logT |
| `hot_star_spectra/cmfgen_<wn\|wc>` | `flux` [Z, logT, wavelength] | `Z`, `logT`, `wavelength` | free wavelength; fixed: 5 Z x 12 logT |
| `hot_star_spectra/post_agb` | `flux` [composition, teff, wavelength] | `composition` (halo, solar), `teff` (K) | free wavelength; fixed: 14 Teff |
| `agb_spectra/<orich\|crich\|crich_aringer>` | `flux` [type, wavelength]; `teff` [type, logZ] (O-rich) or [type] (K) | `type`, `wavelength`, `logZ` (O-rich) | free wavelength; fixed: 9, 5, 9 types, 22 logZ |
| `dusty_agb/<orich\|crich>` | `transfer` [teff, tau, wavelength] | `teff` (K), `tau`, `wavelength` | free wavelength; fixed: 6 Teff x 50 tau |
| `xrb/xsp` | `flux` [feh, log_age, wavelength] (per M_sun, x L_sun on load) | `feh`, `log_age` (log yr), `wavelength` | fixed: 11 [Fe/H] x 10 ages x 15000 |

What "a new library without code changes" covers today:

- A new **spectral library** or **isochrone set** under a new name: all sizes and the solar
  metallicity come from its manifest (`fsps_setup(ctx, zin, 'myisoc', 'myspec')`). The library
  is treated generically: physics that legacy FSPS ties to particular isochrone sets (log Mdot
  for Wolf-Rayet temperatures with `mist`, log g from (M, L, T) for `bsti`, horizontal-branch
  detection for `pdva`, the lower IMF bound for `gnva`) is applied by name only. A nebular grid
  is used if `nebular/<isoc>_<dust|nodust>` exists.
- A replacement of any other library under its fixed name (e.g. new filter curves, a finer
  nebular continuum, other Vega or solar spectra), within the fixed sizes of the table.
- Limits: the (log T, log g) grid of the main spectral library is fixed at 68 x 19 (arrays in
  the context and in the setup-cache tables have that size), and a dust-emission library must
  be named `DL07` or `THEMIS` (the name is a six-character field).

### 3.1 [alpha/Fe]

[alpha/Fe] is a named axis, `afe` (dex), of spectral libraries and isochrones. A library
without it is the [alpha/Fe] = 0 case. With it, `afe` is the first (slowest) axis of the arrays
that depend on it, followed by the metallicity axis:

```json
"axes": {
  "afe": {"units": "[alpha/Fe] (dex)", "values": [-0.2, 0.0, 0.2, 0.4, 0.6], "default": 0.0},
  "feh": {"units": "[Fe/H] (dex); Z = solar_metallicity * 10^[Fe/H]",
          "values": [-2.5, -2.25, -2.0, -1.75, -1.5, -1.25, -1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5]},
  ...
},
"arrays": {
  "flux": {"axes": ["afe", "feh", "logg", "logT", "wavelength"], "shape": [5, 13, 14, 80, 1936], ...}
}
```

- **Metallicity** is given either as `Z` (mass fraction) or as `feh` with the library's
  `solar_metallicity`, whichever the source tabulates; a 2-D grid is the product (afe, Z) or
  (afe, feh) of the two axes.
- **Irregular grids** (not every [alpha/Fe] at every [Fe/H]) keep the rectangular arrays and add
  `available`, an int32 array over [afe, feh] (1 tabulated, 0 absent); absent spectra are NaN,
  absent isochrone blocks have `npoints` = 0. Loading an absent combination stops with a message.
  The upstream alpha-enhanced MIST and C3K grids are regular (13 [Fe/H] x 5 [alpha/Fe]), so they
  need no mask.
- **Selection**: an axis may declare a `default` coordinate. The library is loaded at the
  default unless the name passed to `fsps_setup` selects another value:
  `fsps_setup(ctx, -1, 'mist:afe=0.4', 'c3k_lr:afe=0.4')`. The selection is kept apart from the
  library name (`state%isoc_select`, `state%spec_select`; part of the setup-cache key), so
  name-based physics (`mist`) still applies. The same mechanism selects any other extra axis
  (`name:axis=value,axis=value`, a number or a label).
- **Bit-identical results**: a library with an `afe` axis of one value, or several with the
  default taken, loads exactly the data of the axis-free library (test `unit_standard_library`).

Upstream layouts (converted by `convert_legacy.py --upstream-mist DIR --upstream-c3k DIR`,
tested on synthetic inputs in `tools/fsps_data/test_fsps_data.py`; the fork's data tree does not
contain them):

- MIST: `ISOCHRONES/MIST/zlegend.dat` with lines `(A4,F6.2)` (file tag `m250`, [Fe/H]), and
  `isoc_feh_<tag>_afe_<m2|p0|p2|p4|p6>_vvcrit0.4_full.dat` ([alpha/Fe] = -0.2 ... +0.6),
  solar metallicity 0.0185. Converted to `isochrones/mist_afe`.
- C3K: `SPECTRA/C3K/<name>/` with `<name>.lambda`, `<name>.res`, `<name>_zlegend.dat` ([Fe/H]
  strings such as `-2.50`), `logt.dat` (log T in the first column), `logg.dat` and
  `<name>_feh<feh>_afe<afe>.spec.bin`, each a Fortran (nspec, nlogT, nlogg) FP32 record;
  solar metallicity 0.0185. Note: upstream reads the [Fe/H] strings with format `F4.2` from a
  five-character field, which keeps only the first four characters (`-0.25` becomes -0.2); the
  converter stores the full values.
- Upstream C3K is on its native 80 x 14 (log T, log g) grid, so loading it needs the fixed
  68 x 19 grid of this code made allocatable (see below).

What SSP generation over [alpha/Fe] would need (not done here; the loaders read one slice):

- Spectral library: `state%speclib(nspec, nz, nlogT, nlogg)` gains an afe dimension, as do
  `speclib_agb_smooth`; `speclib_logt`/`speclib_logg` and the matching `fsps_setup_tables_t`
  entries must become allocatable to accept native grids; `NDIM_LOGT`/`NDIM_LOGG` go.
- Isochrones: `mini_isoc`, `mact_isoc`, `logl_isoc`, `logt_isoc`, `logg_isoc`, `ffco_isoc`,
  `phase_isoc`, `lmdot_isoc` (nz, nt, NM), `nmass_isoc` and `timestep_isoc` (nz, nt) gain an afe
  dimension; `zlegend` stays one-dimensional if the [Fe/H] grid is common to all [alpha/Fe]
  (otherwise per afe, with `available`).
- SSP outputs: `spec_ssp_zz(nspec, ntfull, nz)`, `mass_ssp_zz`, `lbol_ssp_zz` (and the
  template-sum weights of `generate_ssp_grid`) gain an afe dimension; the loops over metallicity
  in `fsps_initialization` (spectral library onto the isochrone metallicities) and in the SSP
  driver gain an afe loop. The hot-star, AGB and post-AGB libraries have no [alpha/Fe] axis and
  stay shared.
- Batched CSP with per-model [alpha/Fe]: the SSP basis (`ssp_basis_spec`, `csp_ssp_grid`,
  emission-line grids) over (Z, afe) and a bilinear weight per model in (log Z, afe), in
  `compute_csp_scenario` and `compute_csp_batch` alike; the nebular grids have no [alpha/Fe] axis.

## 4. Tools

All in `tools/` (Python >= 3.8; NumPy for everything except `fetch_data.py`):

| Tool | Purpose |
|---|---|
| `fsps_data/fsps_data.py` | Library writer/reader and verification (`LibraryWriter`, `Library`), root index, and a reader of Fortran list-directed text |
| `fsps_data/convert_legacy.py LEGACY_ROOT OUT_ROOT` | Converts all 37 libraries of the legacy layout (35 s); `--only kind[/name]`; `--upstream-mist`, `--upstream-c3k` for the upstream [alpha/Fe] layouts |
| `fsps_data/make_test_library.py STANDARD_ROOT OUT` | Builds the test library of Section 6 |
| `fsps_data/make_afe_variants.py` | [alpha/Fe] variants of the test library for `unit_standard_library` |
| `fsps_data/provenance.py` | Source, version, citations and redistribution note of each library (written into the manifests and DATA_README.md) |
| `fsps_data/pack_data.py ROOT OUTDIR --release data-vK [--previous M]` | Versioned library archives, release manifest, DATA_README.md (Section 6) |
| `fsps_data/upload_release.sh OUTDIR data-vK` | Publishes a packed release as a GitHub release (run by hand after review) |
| `fetch_data.py DEST [--release data-vK]` | Download from the mirrors, verify and unpack a release |
| `fsps_data/test_fsps_data.py` | Tests of the tools (format round trip, list-directed reader, upstream converters) |

Writing a library from Python:

```python
from fsps_data import LibraryWriter
lib = LibraryWriter(root, "stellar_spectra", "mylib", title="My library")
lib.metadata["solar_metallicity"] = 0.0142
lib.axis("Z", "mass fraction", values=z)                 # 1-D sequences
lib.axis("logT", "log10(K)", values=logt)                # must be 68 values (Section 3)
lib.axis("logg", "log10(cm s^-2)", values=logg)          # 19 values
lib.axis("wavelength", "angstrom", data=wave, medium="vacuum")
lib.array("flux", flux.astype("f4"), axes=["Z", "logg", "logT", "wavelength"], units="H_nu")
lib.array("resolution", res, axes=["wavelength"], units="angstrom (FWHM)")
lib.write()
```

then `write_root_index(root, release)` once all libraries are in place.

## 5. How the code selects a format

- The data root is found as before (`SPS_HOME`, then `FSPS_DATA_HOME`, ...); a directory with
  `fsps-data.json` is accepted as a root.
- `FSPS_DATA_FORMAT` = `auto` (default): each library is read from
  `libraries/<kind>/<name>/manifest.json` if that exists, otherwise from the legacy file
  under `data/`. `standard`: every library from the standard format (a missing manifest is an
  error). `legacy`: the legacy readers only. The format is part of the setup-cache key.
- `src/io/fsps_json.f90` parses the manifests; `src/io/fsps_datalib.f90` validates them and
  reads arrays with stream I/O, checking each file's size; `src/io/fsps_io.f90` holds, for each
  library, a legacy and a standard reader that fill the same raw tables, followed by the shared
  processing (interpolation onto the model grids, unit conversions).
- Exactness: `test_data_parity` (suite `data_parity`, `-Dstandard_data_home=<root>`) sets up
  contexts from both layouts in one process and compares every array the setup loads, an SSP
  grid and a CSP: 539 of 539 comparisons identical over seven library combinations (MIST,
  Padova, PARSEC, BaSTI, Geneva, BPASS; MILES, BaSeL, C3K; DL07, THEMIS; nebular, dust and AGN
  emission, X-ray binaries, dusty AGB), in FP64 and FP32 builds. The unit and regression suites
  pass on either layout, and the 35-case matrix (`tools/compare_builds/run_cases.sh` with
  `DATA_HOME=<root> FSPS_DATA_FORMAT=standard`) gives byte-identical outputs.
- Converting introduced no numerical change: text tables are parsed by the converter with the
  same correctly rounded conversion as Fortran list-directed input, and values derived from
  constants are stored as the compile-time constants of the legacy code evaluate (e.g. the
  X-ray binary ages, correctly rounded log10). The legacy C3K binaries contain NaN fluxes
  (between 1012 and 6072 per metallicity file); they are carried over unchanged.
- Bug fixed on the way: the X-ray-irradiated nebular continuum (BPASS with
  `add_xrb_emission = 1`) was read misaligned (15000 wavelengths read from a 1963-point grid,
  7-value spectra); it is now read like the regular grid.

## 6. Test library, data releases and CI

- **Test library**: `tests/data/fsps-test-data` (6.0 MB, made by `make_test_library.py`):
  everything a default setup (MIST + MILES + DL07, nebular emission) reads, at low resolution
  (two metallicities, every fifth isochrone age and fourth point, every 25th wavelength of
  MILES, 40 points per filter, FP32). `unit_standard_library` runs a setup, an SSP grid and a
  CSP from it alone, plus the [alpha/Fe] checks; with the other unit tests (which build their own
  small files) the unit suite needs no external data. It is not for science.
- **Releases** (hosting decided: GitHub releases of `elijahmathews/opensps-data`, with Zenodo
  snapshots of tagged data releases later; the repository is not created yet):
  - Every library is an immutable, separately versioned archive,
    `<kind>-<name>-v<N>.tar.zst` (`.tar.gz` without zstd), holding `libraries/<kind>/<name>/`;
    `N` is the manifest's `version`, which must be raised whenever a library's files change.
  - A data release (tag `data-vK`) is a small release manifest `data-vK.json` plus a base
    archive `base-data-vK.tar.zst` (root index, small text tables, `DATA_README.md`). The
    manifest lists each archive's library, version, size, SHA-256, the tag it was uploaded to and
    its mirror URLs: `https://github.com/elijahmathews/opensps-data/releases/download/<tag>/<file>`
    first, a Zenodo record second (placeholder, TODO). A release reuses unchanged archives of an
    earlier one (`pack_data.py --previous`) instead of uploading them again.
  - `pack_data.py ROOT OUTDIR --release data-v1` writes the archives, the manifest,
    `DATA_README.md` (what each library is and what to cite, from `tools/fsps_data/provenance.py`;
    a copy is `tools/fsps_data/DATA_README.md`) and `UPLOAD.txt`;
    `tools/fsps_data/upload_release.sh OUTDIR data-v1` publishes it with `gh release create` and
    `gh release upload` after review (never run by the tools).
  - `fetch_data.py DEST --release data-v1` (or `--release-manifest`, `--mirror DIR`) tries the
    mirrors in order, checks archive and file checksums and unpacks; `FSPS_DATA_HOME` then points
    at DEST.
  - Measured (2026-10-20): `data-v1` is 38 archives, 525 MB with zstd (C3K optional); a
    `data-v2` with one library re-versioned reuses the other 36 archives and uploads two files.
    Fetching either through local-directory mirrors gives a root identical to the source.
  - Provenance: every manifest's `provenance` holds source, version, what to cite and the
    redistribution note, taken only from the FSPS sources (README, manual, revision history, data
    headers and READMEs, code comments); these mostly give author and year, so full references,
    versions, URLs and licenses are marked TODO for the maintainers.
- **CI** today: the full-data job converts the checked-out legacy tree and runs `data_parity`;
  all jobs run the test library in the unit suite. Once the data live outside the repository,
  the regression and legacy-parity jobs would fetch the release instead of checking out `data/`:

  ```yaml
  - uses: actions/cache@v4
    with:
      path: ~/.cache/fsps-data
      key: fsps-data-${{ hashFiles('tools/fetch_data.py') }}-<release>
  - run: python3 tools/fetch_data.py fsps-data --release-manifest <URL> --cache ~/.cache/fsps-data
  - run: FSPS_DATA_HOME=$PWD/fsps-data FSPS_DATA_FORMAT=standard meson test -C build --suite regression
  ```

  The cache holds the verified archives (~0.6 GB without C3K), so a run downloads them once per
  release.
- **Size without the data** (measured 2026-10-19): the files tracked outside `data/` are 70 MB
  (68 MB checked out; 49 MB of git objects in a fresh repository after `git gc`). 54 MB of
  that are the regression references `tests/data/sps_ref_*.bin` (FP64 spectra, which barely
  compress), 6 MB the test library. The references are only usable with the full data, so they
  could move into the data release (or be stored as differences from a smaller set), leaving
  about 15 MB of code, tests and test data.
