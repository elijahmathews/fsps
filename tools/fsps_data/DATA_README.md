# OpenSPS / FSPS data release data-v1

Libraries for the modernized FSPS code in the standard data format (doc/DATA_FORMAT.md of the code). Unpack the archives into one directory (tools/fetch_data.py does this and verifies every checksum) and point `FSPS_DATA_HOME` at it.

## What to cite

Any use: the FSPS papers, Conroy, Gunn & White 2009, ApJ, 699, 486 (doi:10.1088/0004-637X/699/1/486); Conroy & Gunn 2010, ApJ, 712, 833 (doi:10.1088/0004-637X/712/2/833). In addition, cite the sources of the libraries you use (below). Entries marked TODO are incomplete in the FSPS sources this release was made from.

Redistribution: these files have been distributed publicly with FSPS (https://github.com/cconroy20/fsps). No separate license is stated there; the terms of the original data apply.

| Library | Version | Source | Cite |
|---|---|---|---|
| `agb_spectra/crich` (archive v1) | TODO | Averaged empirical C-rich TP-AGB spectra (data/spectra/AGB_spectra/README) | Lancon & Mouhcine 2002 (AGB_spectra/README; TODO: full reference) |
| `agb_spectra/crich_aringer` (archive v1) | TODO | Synthetic carbon-star spectra used to extrapolate the C-rich TP-AGB spectra redward (FSPS manual, section 1.3) | Aringer et al. 2009 (FSPS manual, doc/REVISION_HISTORY; TODO: full reference) |
| `agb_spectra/orich` (archive v1) | TODO | Averaged empirical O-rich TP-AGB spectra (data/spectra/AGB_spectra/README), extrapolated as described in the FSPS manual (section 1.3) | Lancon & Mouhcine 2002 (AGB_spectra/README; TODO: full reference)<br>Lancon & Wood 2002 (FSPS manual, extrapolation; TODO: full reference) |
| `agn_torus/clumpy` (archive v1) | TODO | Clumpy AGN torus models, 'reformatted for FSPS' (file header of Nenkova08_y010_torusg_n10_q2.0.dat): Y = 10, sigma = 10 (torusg), N0 = 10, q = 2.0; a small subset of the grid | Nenkova et al. 2008 (file header, FSPS manual; TODO: full reference) |
| `attenuation/dirty` (archive v1) | TODO | Radiative-transfer attenuation curves for dusty galactic environments, MW and SMC grains, dusty/shell/cloudy geometries, homogeneous and clumpy media (files alldirty_{h,c}.dat) | Witt & Gordon 2000 (file header, FSPS manual; TODO: full reference) |
| `attenuation/smc_bar` (archive v1) | TODO | SMC bar extinction curve (file Gordon03_table4.dat, table 4 of the source) | Gordon et al. 2003 (FSPS manual; TODO: full reference) |
| `dust_emission/DL07` (archive v1) | MW3.1 | Draine & Li 2007 dust emission models, Milky Way grain mix R_V = 3.1 (files DL07_MW3.1_*; header: L_nu per H nucleon, 3.3 um PAH feature reduced by 50%) | Draine & Li 2007 (FSPS manual; TODO: full reference) |
| `dust_emission/THEMIS` (archive v1) | TODO | THEMIS dust emission model (files THEMIS_MW3.1_*) | Jones et al. 2013 (FSPS manual; TODO: full reference)<br>Jones et al. 2017 (FSPS manual; TODO: full reference) |
| `dusty_agb/crich` (archive v1) | TODO | Circumstellar dust shells of C-rich AGB stars computed with the radiative transfer code DUSTY (FSPS manual, section 1.11) | Villaume et al. 2015 (FSPS manual; TODO: full reference) |
| `dusty_agb/orich` (archive v1) | TODO | Circumstellar dust shells of O-rich AGB stars computed with the radiative transfer code DUSTY (FSPS manual, section 1.11) | Villaume et al. 2015 (FSPS manual; TODO: full reference) |
| `filters/allfilters` (archive v1) | 159 filters | Filter transmission curves of data/allfilters.dat; the source of each curve is in its 'description' (from data/FILTER_LIST) | The references of the individual filter curves (per-filter descriptions) |
| `hot_star_spectra/cmfgen_wc` (archive v1) | TODO | CMFGEN Wolf-Rayet spectra, carbon sequence | Smith et al. 2002 (doc/REVISION_HISTORY, for the OB and WR spectra added in v2.2; TODO: confirm that it applies to these files and give the full reference) |
| `hot_star_spectra/cmfgen_wn` (archive v1) | TODO | CMFGEN Wolf-Rayet spectra, nitrogen sequence | Smith et al. 2002 (doc/REVISION_HISTORY, for the OB and WR spectra added in v2.2; TODO: confirm that it applies to these files and give the full reference) |
| `hot_star_spectra/post_agb` (archive v1) | TODO | Synthetic hydrogen-deficient post-AGB spectra, halo and solar compositions | Rauch 2003 (doc/REVISION_HISTORY; TODO: full reference) |
| `hot_star_spectra/wmbasic` (archive v1) | TODO | WMBasic O-star spectral grid 'courtesy of JJ Eldridge, Mason Ng, and Georgie Taylor' (doc/REVISION_HISTORY) | TODO: reference for the WMBasic grid (not given in the FSPS sources) |
| `isochrones/bsti` (archive v1) | TODO | BaSTI isochrones (solar metallicity 0.020) | TODO: BaSTI references (not given in the FSPS sources) |
| `isochrones/gnva` (archive v1) | TODO | Geneva isochrones, 'older Geneva models' (doc/REVISION_HISTORY; solar metallicity 0.020) | TODO: Geneva references (not given in the FSPS sources) |
| `isochrones/mist` (archive v1) | TODO | MIST isochrones (solar metallicity 0.0142, FSPS manual) | TODO: MIST references (README.md points to CITATION.bib, which lists only the FSPS papers) |
| `isochrones/pdva` (archive v1) | TODO | Padova 2007 isochrones (data/isochrones/Padova/Padova2007; solar metallicity 0.019) | TODO: Padova isochrone references (not given in the FSPS sources) |
| `isochrones/prsc` (archive v1) | TODO | PARSEC isochrones (solar metallicity 0.01524) | TODO: PARSEC references (not given in the FSPS sources) |
| `nebular/bpss_dust` (archive v1) | TODO | Cloudy photoionization tables provided by Nell Byler, ionized by the bpss SSPs, with dust in the ionized gas (FSPS manual, section 1.13) | Byler et al. 2017 (FSPS manual; TODO: full reference) |
| `nebular/bpss_nodust` (archive v1) | TODO | Cloudy photoionization tables provided by Nell Byler, ionized by the bpss SSPs, without dust in the ionized gas (FSPS manual, section 1.13) | Byler et al. 2017 (FSPS manual; TODO: full reference) |
| `nebular/bpss_xray_dust` (archive v1) | TODO | Cloudy tables for BPASS SSPs including X-ray binaries in the ionizing spectrum, with dust in the ionized gas (files ZAU_WX_*) | TODO: reference for the X-ray-binary nebular grids (not given in the FSPS sources) |
| `nebular/bpss_xray_nodust` (archive v1) | TODO | Cloudy tables for BPASS SSPs including X-ray binaries in the ionizing spectrum, without dust in the ionized gas (files ZAU_WX_*) | TODO: reference for the X-ray-binary nebular grids (not given in the FSPS sources) |
| `nebular/mist_dust` (archive v1) | TODO | Cloudy photoionization tables provided by Nell Byler, ionized by the mist SSPs, with dust in the ionized gas (FSPS manual, section 1.13) | Byler et al. 2017 (FSPS manual; TODO: full reference) |
| `nebular/mist_nodust` (archive v1) | TODO | Cloudy photoionization tables provided by Nell Byler, ionized by the mist SSPs, without dust in the ionized gas (FSPS manual, section 1.13) | Byler et al. 2017 (FSPS manual; TODO: full reference) |
| `nebular/pdva_dust` (archive v1) | TODO | Cloudy photoionization tables provided by Nell Byler, ionized by the pdva SSPs, with dust in the ionized gas (FSPS manual, section 1.13) | Byler et al. 2017 (FSPS manual; TODO: full reference) |
| `nebular/pdva_nodust` (archive v1) | TODO | Cloudy photoionization tables provided by Nell Byler, ionized by the pdva SSPs, without dust in the ionized gas (FSPS manual, section 1.13) | Byler et al. 2017 (FSPS manual; TODO: full reference) |
| `nebular/prsc_dust` (archive v1) | TODO | Cloudy photoionization tables provided by Nell Byler, ionized by the prsc SSPs, with dust in the ionized gas (FSPS manual, section 1.13) | Byler et al. 2017 (FSPS manual; TODO: full reference) |
| `nebular/prsc_nodust` (archive v1) | TODO | Cloudy photoionization tables provided by Nell Byler, ionized by the prsc SSPs, without dust in the ionized gas (FSPS manual, section 1.13) | Byler et al. 2017 (FSPS manual; TODO: full reference) |
| `reference_spectra/sun` (archive v1) | TODO | Solar spectrum from STScI (file name SUN_STScI.SED), f_nu at 10 pc | TODO: reference and URL of the STScI solar spectrum |
| `reference_spectra/vega` (archive v1) | TODO | Kurucz 1992 model atmosphere for Vega, Teff = 9400 K, log g = 3.90, [M/H] = -0.50 (file header of A0V_KURUCZ_92.SED) | Kurucz 1992 (file header; TODO: full reference) |
| `ssp/bpss` (archive v1) | 2.2 | BPASS v2.2 binary population SSPs, '-bin-imf135all_100' models: Salpeter IMF, 100 M_sun upper mass; lowest metallicity omitted, ages cut at 15.8 Gyr (FSPS manual, section 1.16) | Eldridge, Stanway et al. 2017 (FSPS manual; TODO: full reference) |
| `stellar_spectra/basel` (archive v1) | 3.1 (WLBC) | BaSeL 3.1 theoretical library, WLBC version (file names basel_wlbc_*); based on Kurucz model atmospheres (FSPS manual, acknowledgements) | TODO: BaSeL 3.1 / WLBC references (not given in the FSPS sources) |
| `stellar_spectra/c3k_afe+0.0` (archive v1) | TODO | C3K theoretical spectra, [alpha/Fe] = 0, downsampled from full C3K resolution and interpolated onto the BaSeL log g - log T grid (data/spectra/C3K/readme.md); based on Kurucz atmospheres (FSPS manual: 'the C3K spectral library under development') | TODO: C3K reference (none given in the FSPS sources) |
| `stellar_spectra/miles` (archive v1) | TODO | MILES empirical stellar library, combined with the BaSeL library outside the MILES wavelength range and stellar-parameter coverage (doc/REVISION_HISTORY) | TODO: MILES library references (not given in the FSPS sources)<br>TODO: BaSeL references (BaSeL acknowledged in the FSPS manual, no reference given) |
| `xrb/xsp` (archive v1) | TODO | X-ray binary spectra of simple populations (11 [Fe/H] x 10 ages) | TODO: reference (none in the FSPS sources) |

## Layout

Each library is one archive, `<kind>-<name>-v<N>.tar.zst` (or `.tar.gz`), holding `libraries/<kind>/<name>/` (a `manifest.json` and raw little-endian arrays); the base archive holds `fsps-data.json` and the small text tables in `data/`. Archives are immutable: a changed library gets a new version and a new archive, and a release manifest `<release>.json` lists, for every library, its version, SHA-256, size and download URLs (mirrors), possibly pointing at archives of an earlier release.

## Open provenance questions (TODO)

- `agb_spectra/crich`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `agb_spectra/crich_aringer`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `agb_spectra/orich`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `agn_torus/clumpy`: URL; license or redistribution terms of the original data (not stated in the FSPS sources)
- `attenuation/dirty`: URL; license or redistribution terms of the original data (not stated in the FSPS sources)
- `attenuation/smc_bar`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `dust_emission/DL07`: URL; license or redistribution terms of the original data (not stated in the FSPS sources)
- `dust_emission/THEMIS`: URL; version; license or redistribution terms of the original data (not stated in the FSPS sources)
- `dusty_agb/crich`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `dusty_agb/orich`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `filters/allfilters`: per-filter references where the description gives none; license or redistribution terms of the original data (not stated in the FSPS sources)
- `hot_star_spectra/cmfgen_wc`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `hot_star_spectra/cmfgen_wn`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `hot_star_spectra/post_agb`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `hot_star_spectra/wmbasic`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `isochrones/bsti`: BaSTI version and URL; license or redistribution terms of the original data (not stated in the FSPS sources)
- `isochrones/gnva`: URL; license or redistribution terms of the original data (not stated in the FSPS sources)
- `isochrones/mist`: MIST version (rotation, vvcrit) and URL; license or redistribution terms of the original data (not stated in the FSPS sources)
- `isochrones/pdva`: URL; license or redistribution terms of the original data (not stated in the FSPS sources)
- `isochrones/prsc`: PARSEC version and URL; license or redistribution terms of the original data (not stated in the FSPS sources)
- `nebular/bpss_dust`: table version; license or redistribution terms of the original data (not stated in the FSPS sources)
- `nebular/bpss_nodust`: table version; license or redistribution terms of the original data (not stated in the FSPS sources)
- `nebular/bpss_xray_dust`: origin of the ZAU_WX tables; license or redistribution terms of the original data (not stated in the FSPS sources)
- `nebular/bpss_xray_nodust`: origin of the ZAU_WX tables; license or redistribution terms of the original data (not stated in the FSPS sources)
- `nebular/mist_dust`: table version; license or redistribution terms of the original data (not stated in the FSPS sources)
- `nebular/mist_nodust`: table version; license or redistribution terms of the original data (not stated in the FSPS sources)
- `nebular/pdva_dust`: table version; license or redistribution terms of the original data (not stated in the FSPS sources)
- `nebular/pdva_nodust`: table version; license or redistribution terms of the original data (not stated in the FSPS sources)
- `nebular/prsc_dust`: table version; license or redistribution terms of the original data (not stated in the FSPS sources)
- `nebular/prsc_nodust`: table version; license or redistribution terms of the original data (not stated in the FSPS sources)
- `reference_spectra/sun`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `reference_spectra/vega`: license or redistribution terms of the original data (not stated in the FSPS sources)
- `ssp/bpss`: URL; license or redistribution terms of the original data (not stated in the FSPS sources)
- `stellar_spectra/basel`: URL; license or redistribution terms of the original data (not stated in the FSPS sources)
- `stellar_spectra/c3k_afe+0.0`: C3K version; license or redistribution terms of the original data (not stated in the FSPS sources)
- `stellar_spectra/miles`: MILES version and URL; license or redistribution terms of the original data (not stated in the FSPS sources)
- `xrb/xsp`: origin of the files in data/spectra/xrb; license or redistribution terms of the original data (not stated in the FSPS sources)
