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
| Attenuation (`load_attenuation_curves`) | `dust/alldirty_{h,c}.dat`, `dust/Gordon03_table4.dat` | WG00: 3 header lines, 18 blocks (tau_V) of 25 rows: wavelength, tau_V, 6 curves; SMC: 30 rows, wavelength (um) descending | linear interpolation; constant blueward, zero redward |
| WMBasic O stars (`load_wmbasic_spectra`) | `spectra/Hot_spectra/WMBASIC.teff`, `WMBASIC_zlegend.dat`, `WMBASIC_z<Z>.spectra` (12 Z) | rows: wavelength and 3 x 11 spectra (log g outer, T inner); log g grid in the code | onto the grid, then log-flux interpolation onto the isochrone metallicities |
| Wolf-Rayet (`load_wr_spectra`) | `CMFGEN_{WN,WC}.teff`, `CMFGEN_{WN,WC}_Zall.spec` | wavelength line, then 5 Z x 12 T blocks: a (T, Z) label line and a spectrum line | log-flux interpolation in Z, log-log in wavelength |
| Post-AGB (`load_post_agb_spectra`) | `ipagb.teff` (Teff in K), `ipagb_{solar,halo}.spec` | rows: wavelength and 14 spectra | linear onto the grid, floor 1e-308 |
| AGB (`load_agb_spectra`) | `AGB_spectra/{Orich,Crich,Crich_Aringer}.{teff,spec}` | O-rich Teff tabulated per spectral type and 22 [Z/Zsun]; C-rich and Aringer one Teff per spectrum; spectra as rows (wavelength, one column per type) | Teff interpolated to the isochrone Z and logged; spectra linear onto the grid, floor 1e-308 |
| Dusty AGB shells (`load_dusty_agb_spectra`) | `dust/dusty/{Orich,Crich}_dusty.spec` | number of wavelengths, wavelength line, then 6 Teff x 50 tau blocks: (Teff, tau) label and a transfer-function line | linear onto the grid from the first tabulated wavelength |
| X-ray binaries (`load_xrb_spectra`) | `spectra/xrb/xsp.lambda`, `xsp_feh<feh>.spec` (11 [Fe/H]) | 10 lines (ages) of 15000 values each; ages and [Fe/H] in the code | linear onto the grid, floor 1e-308, x L_sun |
| LSF, IMF, SFH (`load_lsf_data`, `fsps_ssp`, `load_tabular_sfh`) | `lsf.dat`, `imf.dat` (or `imf_filename`), `sfh.dat` (or `sfh_filename`) | small user-editable text tables | unchanged by this work |

Unit conventions that are implicit today (they become explicit `units` in the manifests):

- Main stellar library: Eddington flux H_nu (erg s^-1 cm^-2 Hz^-1); the code multiplies by
  16 pi^2 R^2 / L_sun, with R from (M, L, Teff).
- WMBasic, CMFGEN, post-AGB, Aringer: L_nu per unit bolometric luminosity, scaled by L_bol.
- Lancon-Mouhcine O- and C-rich AGB: F_lambda normalized to unit integral over Angstroms,
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
