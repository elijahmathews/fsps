# FSPS data layout

This directory contains all model inputs and lookup tables used by FSPS. As of the current layout, the physical libraries have been moved under `data/` for clarity and discoverability.

## Subdirectories

- `data/isochrones/`: isochrone libraries (BaSTI, Padova, MIST, PARSEC, BPASS, Geneva)
- `data/spectra/`: spectral libraries and reference spectra (MILES, BaSeL, C3K, CKC14, Hot_spectra, AGB_spectra, A0V, Sun, XRB)
- `data/nebular/`: nebular emission tables (continuum + line emission)
- `data/dust/`: dust attenuation curves, emission models, and AGN torus templates

## Other data files

The top-level `data/` directory also contains smaller tables such as `FILTER_LIST`, `allfilters.dat`, IMF/SFH tables, and index definitions. These are read directly by the Fortran code via the resolved data root (`DATA_HOME`, derived from `FSPS_DATA_HOME` or `SPS_HOME`).

## Notes

- `OUTPUTS/` remains at the repository root for compatibility with existing scripts and IDL readers.
- When using `FSPS_OUTPUT_HOME`, ensure `OUTPUTS/` exists under that root (e.g., create `$FSPS_OUTPUT_HOME/OUTPUTS`).
- If you move or mirror this repository, ensure `SPS_HOME` (legacy) or `FSPS_DATA_HOME` points to the install root that contains `data/`.
- When `FSPS_DATA_HOME` is not set, FSPS falls back to standard locations such as `/usr/share/fsps`, `/usr/local/share/fsps`, or `~/.local/share/fsps` (if they contain `data/`).
- `FSPS_OUTPUT_HOME` can be set to override the output directory root if you want outputs outside the repo.
