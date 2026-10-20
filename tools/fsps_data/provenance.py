"""Provenance of the libraries converted from the FSPS data tree.

Everything here is taken from the repository's own sources: README.md and CITATION.bib, the FSPS
manual (doc/MANUAL.pdf), doc/REVISION_HISTORY, the data files' headers and READMEs, and the code
comments of the readers. Those sources mostly give author and year only; full references,
versions, URLs and licenses they do not state are marked TODO for the maintainers to complete.
Nothing here is a license grant: the redistribution note records only that FSPS has distributed
these files publicly.

Each entry: source (who produced the data / project), version, cite (list of strings),
redistribution, notes, todo (list of open questions). Used by convert_legacy.py (written into
each manifest's "provenance") and by write_data_readme() (DATA_README.md of a data release).
"""

FSPS_CITE = [
    "Conroy, Gunn & White 2009, ApJ, 699, 486 (doi:10.1088/0004-637X/699/1/486)",
    "Conroy & Gunn 2010, ApJ, 712, 833 (doi:10.1088/0004-637X/712/2/833)",
]
FSPS_URL = "https://github.com/cconroy20/fsps"
REDISTRIBUTION = ("Distributed publicly with FSPS (" + FSPS_URL + "), from whose data tree this library was "
                  "converted. No separate license is stated in the FSPS distribution; the terms of the original "
                  "data apply.")
TODO_LICENSE = "license or redistribution terms of the original data (not stated in the FSPS sources)"

P = {}


def _e(key, source, cite, version="TODO", notes="", todo=(), url="TODO"):
    P[key] = {"source": source, "url": url, "version": version, "cite": list(cite),
              "cite_also": FSPS_CITE, "redistribution": REDISTRIBUTION, "notes": notes,
              "todo": list(todo) + [TODO_LICENSE]}


# --- stellar spectra
_e("stellar_spectra/miles", "MILES empirical stellar library, combined with the BaSeL library outside the MILES "
   "wavelength range and stellar-parameter coverage (doc/REVISION_HISTORY)",
   ["TODO: MILES library references (not given in the FSPS sources)",
    "TODO: BaSeL references (BaSeL acknowledged in the FSPS manual, no reference given)"],
   todo=["MILES version and URL"])
_e("stellar_spectra/basel", "BaSeL 3.1 theoretical library, WLBC version (file names basel_wlbc_*); based on "
   "Kurucz model atmospheres (FSPS manual, acknowledgements)",
   ["TODO: BaSeL 3.1 / WLBC references (not given in the FSPS sources)"],
   version="3.1 (WLBC)", todo=["URL"])
_e("stellar_spectra/c3k_afe+0.0", "C3K theoretical spectra, [alpha/Fe] = 0, downsampled from full C3K resolution and "
   "interpolated onto the BaSeL log g - log T grid (data/spectra/C3K/readme.md); based on Kurucz atmospheres "
   "(FSPS manual: 'the C3K spectral library under development')",
   ["TODO: C3K reference (none given in the FSPS sources)"],
   todo=["C3K version"])
# --- isochrones
_e("isochrones/mist", "MIST isochrones (solar metallicity 0.0142, FSPS manual)",
   ["TODO: MIST references (README.md points to CITATION.bib, which lists only the FSPS papers)"],
   todo=["MIST version (rotation, vvcrit) and URL"])
_e("isochrones/pdva", "Padova 2007 isochrones (data/isochrones/Padova/Padova2007; solar metallicity 0.019)",
   ["TODO: Padova isochrone references (not given in the FSPS sources)"], todo=["URL"])
_e("isochrones/prsc", "PARSEC isochrones (solar metallicity 0.01524)",
   ["TODO: PARSEC references (not given in the FSPS sources)"], todo=["PARSEC version and URL"])
_e("isochrones/bsti", "BaSTI isochrones (solar metallicity 0.020)",
   ["TODO: BaSTI references (not given in the FSPS sources)"], todo=["BaSTI version and URL"])
_e("isochrones/gnva", "Geneva isochrones, 'older Geneva models' (doc/REVISION_HISTORY; solar metallicity 0.020)",
   ["TODO: Geneva references (not given in the FSPS sources)"], todo=["URL"])
_e("ssp/bpss", "BPASS v2.2 binary population SSPs, '-bin-imf135all_100' models: Salpeter IMF, 100 M_sun upper "
   "mass; lowest metallicity omitted, ages cut at 15.8 Gyr (FSPS manual, section 1.16)",
   ["Eldridge, Stanway et al. 2017 (FSPS manual; TODO: full reference)"], version="2.2", todo=["URL"])
# --- nebular
for _iso in ("mist", "pdva", "prsc", "bpss"):
    for _d, _dn in (("dust", "with"), ("nodust", "without")):
        _e(f"nebular/{_iso}_{_d}", f"Cloudy photoionization tables provided by Nell Byler, ionized by the {_iso} SSPs, "
           f"{_dn} dust in the ionized gas (FSPS manual, section 1.13)",
           ["Byler et al. 2017 (FSPS manual; TODO: full reference)"], todo=["table version"])
    if _iso == "bpss":
        for _d, _dn in (("dust", "with"), ("nodust", "without")):
            _e(f"nebular/bpss_xray_{_d}", "Cloudy tables for BPASS SSPs including X-ray binaries in the ionizing "
               f"spectrum, {_dn} dust in the ionized gas (files ZAU_WX_*)",
               ["TODO: reference for the X-ray-binary nebular grids (not given in the FSPS sources)"],
               todo=["origin of the ZAU_WX tables"])
# --- dust
_e("dust_emission/DL07", "Draine & Li 2007 dust emission models, Milky Way grain mix R_V = 3.1 (files DL07_MW3.1_*; "
   "header: L_nu per H nucleon, 3.3 um PAH feature reduced by 50%)",
   ["Draine & Li 2007 (FSPS manual; TODO: full reference)"], version="MW3.1", todo=["URL"])
_e("dust_emission/THEMIS", "THEMIS dust emission model (files THEMIS_MW3.1_*)",
   ["Jones et al. 2013 (FSPS manual; TODO: full reference)", "Jones et al. 2017 (FSPS manual; TODO: full reference)"],
   todo=["URL", "version"])
_e("agn_torus/clumpy", "Clumpy AGN torus models, 'reformatted for FSPS' (file header of "
   "Nenkova08_y010_torusg_n10_q2.0.dat): Y = 10, sigma = 10 (torusg), N0 = 10, q = 2.0; a small subset of the grid",
   ["Nenkova et al. 2008 (file header, FSPS manual; TODO: full reference)"], todo=["URL"])
_e("attenuation/dirty", "Radiative-transfer attenuation curves for dusty galactic environments, MW and SMC grains, "
   "dusty/shell/cloudy geometries, homogeneous and clumpy media (files alldirty_{h,c}.dat)",
   ["Witt & Gordon 2000 (file header, FSPS manual; TODO: full reference)"], todo=["URL"])
_e("attenuation/smc_bar", "SMC bar extinction curve (file Gordon03_table4.dat, table 4 of the source)",
   ["Gordon et al. 2003 (FSPS manual; TODO: full reference)"])
_e("dusty_agb/orich", "Circumstellar dust shells of O-rich AGB stars computed with the radiative transfer code DUSTY "
   "(FSPS manual, section 1.11)", ["Villaume et al. 2015 (FSPS manual; TODO: full reference)"])
_e("dusty_agb/crich", "Circumstellar dust shells of C-rich AGB stars computed with the radiative transfer code DUSTY "
   "(FSPS manual, section 1.11)", ["Villaume et al. 2015 (FSPS manual; TODO: full reference)"])
# --- hot stars, AGB
_e("hot_star_spectra/wmbasic", "WMBasic O-star spectral grid 'courtesy of JJ Eldridge, Mason Ng, and Georgie Taylor' "
   "(doc/REVISION_HISTORY)", ["TODO: reference for the WMBasic grid (not given in the FSPS sources)"])
_e("hot_star_spectra/cmfgen_wn", "CMFGEN Wolf-Rayet spectra, nitrogen sequence",
   ["Smith et al. 2002 (doc/REVISION_HISTORY, for the OB and WR spectra added in v2.2; TODO: confirm that it applies "
    "to these files and give the full reference)"])
_e("hot_star_spectra/cmfgen_wc", "CMFGEN Wolf-Rayet spectra, carbon sequence",
   ["Smith et al. 2002 (doc/REVISION_HISTORY, for the OB and WR spectra added in v2.2; TODO: confirm that it applies "
    "to these files and give the full reference)"])
_e("hot_star_spectra/post_agb", "Synthetic hydrogen-deficient post-AGB spectra, halo and solar compositions",
   ["Rauch 2003 (doc/REVISION_HISTORY; TODO: full reference)"])
_e("agb_spectra/orich", "Averaged empirical O-rich TP-AGB spectra (data/spectra/AGB_spectra/README), extrapolated as "
   "described in the FSPS manual (section 1.3)",
   ["Lancon & Mouhcine 2002 (AGB_spectra/README; TODO: full reference)",
    "Lancon & Wood 2002 (FSPS manual, extrapolation; TODO: full reference)"])
_e("agb_spectra/crich", "Averaged empirical C-rich TP-AGB spectra (data/spectra/AGB_spectra/README)",
   ["Lancon & Mouhcine 2002 (AGB_spectra/README; TODO: full reference)"])
_e("agb_spectra/crich_aringer", "Synthetic carbon-star spectra used to extrapolate the C-rich TP-AGB spectra redward "
   "(FSPS manual, section 1.3)", ["Aringer et al. 2009 (FSPS manual, doc/REVISION_HISTORY; TODO: full reference)"])
_e("xrb/xsp", "X-ray binary spectra of simple populations (11 [Fe/H] x 10 ages)",
   ["TODO: reference (none in the FSPS sources)"], todo=["origin of the files in data/spectra/xrb"])
# --- calibration
_e("filters/allfilters", "Filter transmission curves of data/allfilters.dat; the source of each curve is in its "
   "'description' (from data/FILTER_LIST)",
   ["The references of the individual filter curves (per-filter descriptions)"], version="159 filters",
   todo=["per-filter references where the description gives none"])
_e("reference_spectra/vega", "Kurucz 1992 model atmosphere for Vega, Teff = 9400 K, log g = 3.90, [M/H] = -0.50 "
   "(file header of A0V_KURUCZ_92.SED)", ["Kurucz 1992 (file header; TODO: full reference)"])
_e("reference_spectra/sun", "Solar spectrum from STScI (file name SUN_STScI.SED), f_nu at 10 pc",
   ["TODO: reference and URL of the STScI solar spectrum"])

# The upstream alpha-enhanced layouts (converted only on request)
_e("isochrones/mist_afe", "MIST isochrones with [alpha/Fe] variations (upstream FSPS layout isoc_feh_*_afe_*)",
   ["TODO: MIST references"], todo=["MIST version"])


def provenance(kind, name):
    key = f"{kind}/{name}"
    if key in P:
        return dict(P[key])
    if kind == "stellar_spectra" and name.startswith("c3k"):
        d = dict(P["stellar_spectra/c3k_afe+0.0"])
        d["source"] = f"C3K theoretical spectra ({name}), upstream FSPS layout"
        return d
    return {"source": "TODO", "url": "TODO", "version": "TODO", "cite": ["TODO"], "cite_also": FSPS_CITE,
            "redistribution": REDISTRIBUTION, "notes": "", "todo": ["everything"]}


def write_data_readme(path, release, libraries, layout_note=""):
    """Writes the release's DATA_README.md: what each library is and what to cite."""
    lines = [
        f"# OpenSPS / FSPS data release {release}",
        "",
        "Libraries for the modernized FSPS code in the standard data format (doc/DATA_FORMAT.md of the code). "
        "Unpack the archives into one directory (tools/fetch_data.py does this and verifies every checksum) and point "
        "`FSPS_DATA_HOME` at it.",
        "",
        "## What to cite",
        "",
        "Any use: the FSPS papers, " + "; ".join(FSPS_CITE) + ". In addition, cite the sources of the libraries "
        "you use (below). Entries marked TODO are incomplete in the FSPS sources this release was made from.",
        "",
        "Redistribution: these files have been distributed publicly with FSPS (" + FSPS_URL + "). No separate "
        "license is stated there; the terms of the original data apply.",
        "",
        "| Library | Version | Source | Cite |",
        "|---|---|---|---|",
    ]
    for kind, name, ver in libraries:
        p = provenance(kind, name)
        lines.append(f"| `{kind}/{name}` (archive v{ver}) | {p['version']} | {p['source']} | "
                     + "<br>".join(p["cite"]) + " |")
    lines += ["", "## Layout", "", layout_note or (
        "Each library is one archive, `<kind>-<name>-v<N>.tar.zst` (or `.tar.gz`), holding "
        "`libraries/<kind>/<name>/` (a `manifest.json` and raw little-endian arrays); the base archive holds "
        "`fsps-data.json` and the small text tables in `data/`. Archives are immutable: a changed library gets a new "
        "version and a new archive, and a release manifest `<release>.json` lists, for every library, its version, "
        "SHA-256, size and download URLs (mirrors), possibly pointing at archives of an earlier release."),
        "", "## Open provenance questions (TODO)", ""]
    for kind, name, _ in libraries:
        p = provenance(kind, name)
        lines.append(f"- `{kind}/{name}`: " + "; ".join(p["todo"]))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
