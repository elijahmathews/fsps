# FSPS data layout

This directory contains all model inputs and lookup tables used by FSPS. As of the current layout, the physical libraries have been moved under `data/` for clarity and discoverability.

## Subdirectories

- `data/isochrones/`: isochrone libraries (BaSTI, Padova, MIST, PARSEC, BPASS, Geneva)
- `data/spectra/`: spectral libraries and reference spectra (MILES, BaSeL, C3K, CKC14, Hot_spectra, AGB_spectra, A0V, Sun, XRB)
- `data/nebular/`: nebular emission tables (continuum + line emission)
- `data/dust/`: dust attenuation curves, emission models, and AGN torus templates

## Other data files

The top-level `data/` directory also contains smaller tables such as `FILTER_LIST`, `allfilters.dat`, IMF/SFH tables, and index definitions. These are read directly by the Fortran code via `$SPS_HOME/data/...`.

## Notes

- `OUTPUTS/` remains at the repository root for compatibility with existing scripts and IDL readers.
- If you move or mirror this repository, ensure `SPS_HOME` points to the repository root so relative paths resolve correctly.
