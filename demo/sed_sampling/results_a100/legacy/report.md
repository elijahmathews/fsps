# Fork vs legacy FSPS on 1000 demo parameter vectors

## Time per vector (all outputs incl. spectra, per 1 Msun formed)

| Code | µs per vector | vs legacy |
|---|---|---|
| legacy aligned_rv | 55,570.5 | 1x |
| legacy aligned | 55,486.3 | 1x |
| legacy stock | 55,360.6 | 1x |
| fork gpu | 27.4 | 2,026x |
| fork cpu | 8,096.6 | 7x |

Legacy aligned_rv: setup 12.1 s, SSPs of all metallicities 215.6 s (once); per vector ztinterp 0.90 ms, compsp 54.57 ms, getindx 0.096 ms.

Legacy aligned: setup 12.1 s, SSPs of all metallicities 215.4 s (once); per vector ztinterp 0.89 ms, compsp 54.50 ms, getindx 0.096 ms.

Legacy stock: setup 12.2 s, SSPs of all metallicities 216.0 s (once); per vector ztinterp 0.88 ms, compsp 54.38 ms, getindx 0.096 ms.

Fork gpu: setup 20.5 s, first call (SSP bases) 1.5 s, then one call of 1000 vectors 27.4 ms.

Fork cpu: setup 12.5 s, first call (SSP bases) 10.5 s, then one call of 1000 vectors 8096.6 ms.

## Differences (median / 99th percentile / max over vectors)

### fork GPU vs legacy aligned_rv

| Quantity | median | p99 | max |
|---|---|---|---|
| mags (all bands) [mag] | 3.22e-04 | 1.75e-03 | 4.89e-03 |
| mags UV (GALEX) [mag] | 1.46e-04 | 4.54e-04 | 5.68e-04 |
| mags optical (SDSS ugriz) [mag] | 2.65e-04 | 5.49e-04 | 6.00e-04 |
| mags near-IR (2MASS JHKs) [mag] | 3.90e-04 | 5.76e-04 | 6.51e-04 |
| mags JWST NIRCam [mag] | 3.74e-04 | 6.05e-04 | 8.77e-04 |
| mags mid-IR (WISE, IRAC) [mag] | 4.26e-04 | 1.20e-03 | 2.51e-03 |
| mags far-IR (MIPS, PACS, SPIRE) [mag] | 3.73e-04 | 2.10e-03 | 2.58e-03 |
| spectrum 912 A - 2.5 um [rel] | 1.15e-05 | 1.58e-04 | 1.22e-02 |
| spectrum 2.5 - 1000 um (dust, AGN) [rel] | 1.81e-04 | 2.14e-03 | 3.02e-03 |
| spectrum < 912 A (ionizing) [rel] | 2.77e-05 | 8.51e-05 | 9.52e-05 |
| indices [A or mag] | 2.58e-05 | 9.51e-04 | 6.42e-03 |
| emission lines [rel] | 9.99e-06 | 1.53e-04 | 4.11e-04 |
| mfrac [rel] | 1.32e-06 | 2.74e-06 | 4.64e-06 |
| log L_bol [dex] | 5.47e-07 | 5.86e-07 | 5.89e-07 |
| dust mass [rel] | 4.39e-06 | 2.28e-05 | 3.18e-05 |

Magnitudes compared: 144795 of 159000 (with flux in both, within 30 mag of the vector's brightest band). Largest per-band max: Bessell_M 0.0049, SPIRE_250 0.0026, SPIRE_350 0.0026, SCUBA_450WB 0.0026, SPIRE_500 0.0026, IRAS_100 0.0026, PACS_100 0.0026, PACS_160 0.0026

### fork GPU vs legacy aligned

| Quantity | median | p99 | max |
|---|---|---|---|
| mags (all bands) [mag] | 2.63e-03 | 8.64e-03 | 1.23e-02 |
| mags UV (GALEX) [mag] | 1.86e-03 | 8.04e-03 | 8.49e-03 |
| mags optical (SDSS ugriz) [mag] | 3.06e-03 | 8.98e-03 | 1.04e-02 |
| mags near-IR (2MASS JHKs) [mag] | 2.91e-03 | 7.61e-03 | 1.05e-02 |
| mags JWST NIRCam [mag] | 2.76e-03 | 8.30e-03 | 1.03e-02 |
| mags mid-IR (WISE, IRAC) [mag] | 1.51e-03 | 6.36e-03 | 8.17e-03 |
| mags far-IR (MIPS, PACS, SPIRE) [mag] | 4.16e-04 | 2.83e-03 | 3.77e-03 |
| spectrum 912 A - 2.5 um [rel] | 2.24e-03 | 4.99e-03 | 9.61e-03 |
| spectrum 2.5 - 1000 um (dust, AGN) [rel] | 5.01e-04 | 3.54e-03 | 7.96e-03 |
| spectrum < 912 A (ionizing) [rel] | 8.11e-05 | 5.12e-04 | 6.33e-04 |
| indices [A or mag] | 3.45e-05 | 9.87e-04 | 6.42e-03 |
| emission lines [rel] | 1.04e-03 | 5.06e-03 | 8.23e-03 |
| mfrac [rel] | 1.32e-06 | 2.74e-06 | 4.64e-06 |
| log L_bol [dex] | 5.47e-07 | 5.86e-07 | 5.89e-07 |
| dust mass [rel] | 2.56e-04 | 1.22e-03 | 1.56e-03 |

Magnitudes compared: 144801 of 159000 (with flux in both, within 30 mag of the vector's brightest band). Largest per-band max: Bessell_M 0.0123, LSST_y 0.0108, VISTA_J 0.0107, WFCAM_J 0.0107, VISTA_Y 0.0106, NEWFIRM_J3 0.0106, WFCAM_Y 0.0106, 2MASS_J 0.0105

### fork GPU vs legacy stock

| Quantity | median | p99 | max |
|---|---|---|---|
| mags (all bands) [mag] | 2.40e-02 | 1.11e+01 | 2.32e+01 |
| mags UV (GALEX) [mag] | 9.69e-02 | 1.64e+01 | 1.94e+01 |
| mags optical (SDSS ugriz) [mag] | 2.57e-02 | 1.41e+01 | 2.22e+01 |
| mags near-IR (2MASS JHKs) [mag] | 2.46e-02 | 2.08e+00 | 9.98e+00 |
| mags JWST NIRCam [mag] | 2.41e-02 | 5.78e+00 | 1.90e+01 |
| mags mid-IR (WISE, IRAC) [mag] | 2.09e-02 | 1.52e-01 | 4.61e-01 |
| mags far-IR (MIPS, PACS, SPIRE) [mag] | 1.28e-02 | 2.09e-02 | 2.21e-02 |
| spectrum 912 A - 2.5 um [rel] | 1.94e-02 | 3.45e-01 | 2.40e+00 |
| spectrum 2.5 - 1000 um (dust, AGN) [rel] | 8.89e-03 | 1.90e-02 | 1.72e-01 |
| spectrum < 912 A (ionizing) [rel] | 2.71e-01 | 2.21e+01 | 2.99e+01 |
| indices [A or mag] | 1.02e-03 | 7.61e-01 | 2.79e+01 |
| emission lines [rel] | 6.46e-03 | 1.79e-02 | 2.78e-02 |
| mfrac [rel] | 1.32e-06 | 2.74e-06 | 4.64e-06 |
| log L_bol [dex] | 5.47e-07 | 5.86e-07 | 5.89e-07 |
| dust mass [rel] | 9.79e-03 | 1.65e-02 | 1.70e-02 |

Magnitudes compared: 149052 of 159000 (with flux in both, within 30 mag of the vector's brightest band). Largest per-band max: Stromgren_v 23.1690, Steidel_Un 22.6701, I2800 22.2984, SDSS_u 22.2007, Stromgren_y 22.0753, Steidel_Rs 21.4009, Stromgren_u 21.0472, LSST_u 20.9780

### fork CPU vs legacy aligned_rv

| Quantity | median | p99 | max |
|---|---|---|---|
| mags (all bands) [mag] | 3.22e-04 | 1.75e-03 | 4.89e-03 |
| mags UV (GALEX) [mag] | 1.46e-04 | 4.54e-04 | 5.68e-04 |
| mags optical (SDSS ugriz) [mag] | 2.65e-04 | 5.49e-04 | 6.00e-04 |
| mags near-IR (2MASS JHKs) [mag] | 3.90e-04 | 5.76e-04 | 6.51e-04 |
| mags JWST NIRCam [mag] | 3.74e-04 | 6.06e-04 | 8.77e-04 |
| mags mid-IR (WISE, IRAC) [mag] | 4.26e-04 | 1.20e-03 | 2.51e-03 |
| mags far-IR (MIPS, PACS, SPIRE) [mag] | 3.73e-04 | 2.10e-03 | 2.58e-03 |
| spectrum 912 A - 2.5 um [rel] | 1.15e-05 | 1.58e-04 | 1.22e-02 |
| spectrum 2.5 - 1000 um (dust, AGN) [rel] | 1.81e-04 | 2.14e-03 | 3.02e-03 |
| spectrum < 912 A (ionizing) [rel] | 2.77e-05 | 8.51e-05 | 9.53e-05 |
| indices [A or mag] | 2.58e-05 | 9.53e-04 | 6.40e-03 |
| emission lines [rel] | 9.98e-06 | 1.53e-04 | 4.11e-04 |
| mfrac [rel] | 1.32e-06 | 2.74e-06 | 4.64e-06 |
| log L_bol [dex] | 5.47e-07 | 5.86e-07 | 5.89e-07 |
| dust mass [rel] | 4.38e-06 | 2.27e-05 | 3.18e-05 |

Magnitudes compared: 144795 of 159000 (with flux in both, within 30 mag of the vector's brightest band). Largest per-band max: Bessell_M 0.0049, SPIRE_250 0.0026, SPIRE_350 0.0026, SCUBA_450WB 0.0026, SPIRE_500 0.0026, IRAS_100 0.0026, PACS_100 0.0026, PACS_160 0.0026

### fork CPU vs legacy aligned

| Quantity | median | p99 | max |
|---|---|---|---|
| mags (all bands) [mag] | 2.63e-03 | 8.64e-03 | 1.23e-02 |
| mags UV (GALEX) [mag] | 1.86e-03 | 8.04e-03 | 8.49e-03 |
| mags optical (SDSS ugriz) [mag] | 3.06e-03 | 8.98e-03 | 1.04e-02 |
| mags near-IR (2MASS JHKs) [mag] | 2.91e-03 | 7.61e-03 | 1.05e-02 |
| mags JWST NIRCam [mag] | 2.76e-03 | 8.30e-03 | 1.03e-02 |
| mags mid-IR (WISE, IRAC) [mag] | 1.51e-03 | 6.36e-03 | 8.17e-03 |
| mags far-IR (MIPS, PACS, SPIRE) [mag] | 4.16e-04 | 2.83e-03 | 3.77e-03 |
| spectrum 912 A - 2.5 um [rel] | 2.24e-03 | 4.99e-03 | 9.61e-03 |
| spectrum 2.5 - 1000 um (dust, AGN) [rel] | 5.01e-04 | 3.54e-03 | 7.96e-03 |
| spectrum < 912 A (ionizing) [rel] | 8.11e-05 | 5.12e-04 | 6.33e-04 |
| indices [A or mag] | 3.45e-05 | 9.89e-04 | 6.40e-03 |
| emission lines [rel] | 1.04e-03 | 5.06e-03 | 8.23e-03 |
| mfrac [rel] | 1.32e-06 | 2.74e-06 | 4.64e-06 |
| log L_bol [dex] | 5.47e-07 | 5.86e-07 | 5.89e-07 |
| dust mass [rel] | 2.56e-04 | 1.22e-03 | 1.56e-03 |

Magnitudes compared: 144801 of 159000 (with flux in both, within 30 mag of the vector's brightest band). Largest per-band max: Bessell_M 0.0123, LSST_y 0.0108, VISTA_J 0.0107, WFCAM_J 0.0107, VISTA_Y 0.0106, NEWFIRM_J3 0.0106, WFCAM_Y 0.0106, 2MASS_J 0.0105

### fork CPU vs legacy stock

| Quantity | median | p99 | max |
|---|---|---|---|
| mags (all bands) [mag] | 2.40e-02 | 1.11e+01 | 2.32e+01 |
| mags UV (GALEX) [mag] | 9.69e-02 | 1.64e+01 | 1.94e+01 |
| mags optical (SDSS ugriz) [mag] | 2.57e-02 | 1.41e+01 | 2.22e+01 |
| mags near-IR (2MASS JHKs) [mag] | 2.46e-02 | 2.08e+00 | 9.98e+00 |
| mags JWST NIRCam [mag] | 2.41e-02 | 5.78e+00 | 1.90e+01 |
| mags mid-IR (WISE, IRAC) [mag] | 2.09e-02 | 1.52e-01 | 4.61e-01 |
| mags far-IR (MIPS, PACS, SPIRE) [mag] | 1.28e-02 | 2.09e-02 | 2.21e-02 |
| spectrum 912 A - 2.5 um [rel] | 1.94e-02 | 3.45e-01 | 2.40e+00 |
| spectrum 2.5 - 1000 um (dust, AGN) [rel] | 8.89e-03 | 1.90e-02 | 1.72e-01 |
| spectrum < 912 A (ionizing) [rel] | 2.71e-01 | 2.21e+01 | 2.99e+01 |
| indices [A or mag] | 1.02e-03 | 7.61e-01 | 2.79e+01 |
| emission lines [rel] | 6.46e-03 | 1.79e-02 | 2.78e-02 |
| mfrac [rel] | 1.32e-06 | 2.74e-06 | 4.64e-06 |
| log L_bol [dex] | 5.47e-07 | 5.86e-07 | 5.89e-07 |
| dust mass [rel] | 9.79e-03 | 1.65e-02 | 1.70e-02 |

Magnitudes compared: 149052 of 159000 (with flux in both, within 30 mag of the vector's brightest band). Largest per-band max: Stromgren_v 23.1690, Steidel_Un 22.6701, I2800 22.2984, SDSS_u 22.2007, Stromgren_y 22.0753, Steidel_Rs 21.4009, Stromgren_u 21.0472, LSST_u 20.9780

### fork GPU (FP32) vs fork CPU (FP64)

| Quantity | median | p99 | max |
|---|---|---|---|
| mags (all bands) [mag] | 3.74e-08 | 1.84e-07 | 9.51e-07 |
| mags UV (GALEX) [mag] | 4.14e-08 | 1.97e-07 | 3.74e-07 |
| mags optical (SDSS ugriz) [mag] | 3.54e-08 | 1.70e-07 | 3.34e-07 |
| mags near-IR (2MASS JHKs) [mag] | 3.92e-08 | 1.58e-07 | 2.42e-07 |
| mags JWST NIRCam [mag] | 3.62e-08 | 1.50e-07 | 2.95e-07 |
| mags mid-IR (WISE, IRAC) [mag] | 4.09e-08 | 1.63e-07 | 3.36e-07 |
| mags far-IR (MIPS, PACS, SPIRE) [mag] | 5.52e-08 | 2.21e-07 | 3.38e-07 |
| spectrum 912 A - 2.5 um [rel] | 3.22e-08 | 1.71e-07 | 8.95e-07 |
| spectrum 2.5 - 1000 um (dust, AGN) [rel] | 4.66e-08 | 2.03e-07 | 4.73e-07 |
| spectrum < 912 A (ionizing) [rel] | 2.75e-08 | 1.77e-07 | 3.52e-07 |
| indices [A or mag] | 1.89e-07 | 4.01e-06 | 4.48e-05 |
| emission lines [rel] | 6.03e-08 | 2.37e-07 | 4.48e-07 |
| mfrac [rel] | 2.22e-16 | 8.88e-16 | 1.11e-15 |
| log L_bol [dex] | 4.44e-16 | 1.24e-14 | 5.78e-14 |
| dust mass [rel] | 4.21e-08 | 1.54e-07 | 2.03e-07 |

Magnitudes compared: 144795 of 159000 (with flux in both, within 30 mag of the vector's brightest band). Largest per-band max: WFC3_UVIS_F275W 0.0000, Stromgren_b 0.0000, Stromgren_u 0.0000, WFPC2_F336W 0.0000, U 0.0000, WFC3_UVIS_F225W 0.0000, Stromgren_v 0.0000, GALEX_FUV 0.0000
