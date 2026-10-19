#!/usr/bin/env python3
"""Batched CSPs for SED-fitting samplers: N independent parameter vectors -> N SEDs per call.

Every parameter vector in a batch has its own values of all 14 fitting parameters (stellar mass,
metallicity, a 7-bin non-parametric SFH, dust, dust emission, nebular, AGN and redshift), drawn
independently from Prospector-style priors, as a sampler proposes them. The timed path is the
whole step a likelihood needs: parameter arrays -> SFH weights -> FSPS batched CSP (SSP-basis
interpolation in metallicity, nebular emission, dust attenuation and emission, IGM, AGN) ->
redshifting, distance and 159 filter convolutions, 30 spectral indices, 166 line
luminosities, surviving stellar mass, L_bol and dust mass, all back in host (numpy) memory and
scaled to each vector's mass. What is shared across the batch is printed by --describe and in
README.md, with what each output contains.

Modes:
  demo_batch_sed.py --lib build/libfsps.so [--sizes 1,10,100,1000,10000] [--spec] [--smooth 0]
      describe the setup, check per-vector correctness, time the batch, write figures/results.
  demo_batch_sed.py --lib build_cpu/libfsps.so --replay out/vectors.npz
      recompute saved vectors with another build (e.g. CPU FP64) and compare.
"""

import argparse
import json
import os
import platform
import subprocess
import sys
import time

import numpy as np

import fsps_batch

# ---------------------------------------------------------------------------------------------
# Model definition: what is fixed for the whole batch, what each vector varies
# ---------------------------------------------------------------------------------------------

# Libraries chosen at setup (fixed: they define the SSP basis)
LIBRARIES = {"isochrones": "MIST", "stellar spectra": "MILES (FSPS 'miles' spectral library)",
             "dust emission templates": "Draine & Li 2007 (DL07)"}

# Context settings shared by every vector (FSPS names). Integers are switches.
INT_SETTINGS = {
    "imf_type": 1,             # Chabrier 2003 IMF
    "dust_type": 4,            # Kriek & Conroy 2013 attenuation curve (dust_index = slope offset)
    "add_dust_emission": 1,    # energy-balance dust emission (DL07)
    "add_neb_emission": 1,     # nebular emission from the CLOUDY grid, computed per vector
    "add_neb_continuum": 1,    # nebular continuum
    "nebemlineinspec": 1,      # emission lines in the spectrum (and so in the photometry)
    "add_agn_dust": 1,         # AGN torus emission (Nenkova et al. 2008 templates)
    "add_igm_absorption": 1,   # Madau 1995 IGM attenuation at each vector's redshift
    "add_stellar_remnants": 1, # remnants counted in the surviving mass
    "compute_indices": 1,      # Lick / other spectral indices
    "compute_vega_mags": 0,    # AB magnitudes
    "vactoair_flag": 0,        # vacuum wavelengths
}
FLOAT_SETTINGS = {
    "H0": 70.0, "om0": 0.3, "ol0": 0.7,   # flat LambdaCDM: distances, t_univ(z)
    "dust_tesc": 7.0,          # log age (yr) below which stars sit in birth clouds (dust1)
    "frac_nodust": 0.0, "frac_obrun": 0.0, "igm_factor": 1.0,
}
COSMO = (FLOAT_SETTINGS["H0"], FLOAT_SETTINGS["om0"])

# Priors (Prospector-alpha-like, Leja et al. 2017, 2019). Each vector draws all of them.
N_SFH_BINS = 7
PRIORS = {
    "logmass":     ("uniform", 8.0, 12.0, "log10 stellar mass formed [Msun] (scales the outputs)"),
    "logzsol":     ("uniform", -1.98, 0.19, "stellar log(Z/Zsun); SSP basis interpolated per vector"),
    "dust2":       ("uniform", 0.0, 4.0, "diffuse dust optical depth at 5500 A"),
    "dust_ratio":  ("uniform", 0.0, 2.0, "dust1 / dust2 (birth-cloud dust), dust1 passed per vector"),
    "dust_index":  ("uniform", -1.0, 0.4, "attenuation slope offset (Kriek & Conroy)"),
    "fagn":        ("loguniform", 1e-5, 3.0, "AGN torus luminosity / stellar L_bol"),
    "agn_tau":     ("loguniform", 5.0, 150.0, "AGN torus optical depth"),
    "gas_logz":    ("uniform", -2.0, 0.5, "gas-phase log(Z/Zsun)"),
    "gas_logu":    ("uniform", -4.0, -1.0, "log ionization parameter"),
    "duste_umin":  ("uniform", 0.1, 25.0, "DL07 minimum radiation field"),
    "duste_gamma": ("loguniform", 1e-3, 0.15, "DL07 fraction of dust in high radiation field"),
    "duste_qpah":  ("uniform", 0.5, 7.0, "DL07 PAH mass fraction [%]"),
    "zred":        ("uniform", 0.05, 5.0, "redshift (observed frame, IGM, distance, SFH age limit)"),
    "logsfr_ratios": ("student_t", 2.0, 0.3, f"{N_SFH_BINS - 1} log SFR ratios of adjacent bins "
                      "(continuity prior, Student-t df=2, scale 0.3)"),
}
# Parameters handed to FSPS per vector (logmass and the SFH act through the weight matrix)
FSPS_ARRAYS = ["logzsol", "dust2", "dust1", "dust_index", "fagn", "agn_tau", "gas_logz", "gas_logu",
               "duste_umin", "duste_gamma", "duste_qpah", "zred"]


def draw_vectors(rng, n):
    """n independent draws from PRIORS (a sampler's proposals)."""
    p = {}
    for k, (kind, a, b, _) in PRIORS.items():
        if kind == "uniform":
            p[k] = rng.uniform(a, b, n)
        elif kind == "loguniform":
            p[k] = np.exp(rng.uniform(np.log(a), np.log(b), n))
        else:
            p[k] = b * rng.standard_t(a, (n, N_SFH_BINS - 1))
    return p


# ---------------------------------------------------------------------------------------------
# Non-parametric SFH -> SSP weights (part of the timed path)
# ---------------------------------------------------------------------------------------------

try:
    import numexpr as _ne
except ImportError:
    _ne = None
_MASS_EXPR = " + ".join(f"s{b} * where(e > l{b}, where(e < l{b + 1}, e - l{b}, l{b + 1} - l{b}), 0.0)"
                        for b in range(N_SFH_BINS))


class SfhWeights:
    """Continuity SFH (Leja et al. 2019): 7 bins in lookback time, 0-30 Myr, 30-100 Myr, then
    log-spaced up to 0.85 t_univ(z) and a last bin to t_univ(z); SFR constant in each bin,
    log SFR ratios between adjacent bins. Weight of SSP age k = mass formed in the lookback-time
    interval that age represents (edges at the log-midpoints of the SSP grid), normalized to
    1 Msun formed per vector; ages older than t_univ(z) get none."""

    def __init__(self, log_age):
        t = 10.0 ** log_age
        mid = 10.0 ** (0.5 * (log_age[1:] + log_age[:-1]))
        self.lo = np.concatenate([[0.0], mid])          # SSP interval edges, yr
        self.hi = np.concatenate([mid, [t[-1]]])
        self.edges = np.concatenate([self.lo, self.hi[-1:]])
        # t_univ(z) table (flat LambdaCDM of FLOAT_SETTINGS), interpolated per vector
        h0, om = COSMO
        z = np.linspace(0.0, 20.0, 4001)
        a = 1.0 / (1.0 + z)
        aa = np.linspace(1e-6, 1.0, 200001)
        integrand = 1.0 / (aa * np.sqrt(om / aa**3 + (1.0 - om)))
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(aa))])
        hubble_time_yr = 977.8e9 / h0
        self._z, self._tuniv = z, hubble_time_yr * np.interp(a, aa, cum)

    def t_univ(self, zred):
        return np.interp(zred, self._z, self._tuniv)

    def __call__(self, zred, logsfr_ratios):
        n = zred.size
        tu = np.minimum(self.t_univ(zred), self.hi[-1])
        edges = np.empty((n, N_SFH_BINS + 1))
        edges[:, 0], edges[:, 1], edges[:, 2] = 0.0, 3e7, 1e8
        edges[:, 3:N_SFH_BINS] = np.exp(np.linspace(np.log(1e8), np.log(0.85 * tu), N_SFH_BINS - 2)[1:].T)
        edges[:, N_SFH_BINS] = tu
        # SFR per bin (youngest first): sfr_{i+1} = sfr_i * 10**(-ratio_i)
        logsfr = np.concatenate([np.zeros((n, 1)), -np.cumsum(logsfr_ratios, axis=1)], axis=1)
        sfr = 10.0 ** (logsfr - logsfr.max(axis=1, keepdims=True))
        # Mass formed within lookback time t is piecewise linear (SFR constant per bin): evaluate
        # it at the SSP interval edges (ages beyond t_univ fall in no bin); weights are its
        # increments. sum_b sfr_b * clip(t - edge_b, 0, width_b), one multithreaded pass (numexpr)
        loc = {"e": self.edges[None, :]}
        for b in range(N_SFH_BINS + 1):
            loc[f"l{b}"] = edges[:, b:b + 1]
        for b in range(N_SFH_BINS):
            loc[f"s{b}"] = sfr[:, b:b + 1]
        if _ne is not None:
            mass_le = _ne.evaluate(_MASS_EXPR, local_dict=loc)
        else:
            mass_le = sum(sfr[:, b:b + 1] * np.clip(self.edges[None, :] - edges[:, b:b + 1], 0.0,
                                                    edges[:, b + 1:b + 2] - edges[:, b:b + 1])
                          for b in range(N_SFH_BINS))
        w = np.diff(mass_le, axis=1)
        w /= w.sum(axis=1, keepdims=True)
        return np.ascontiguousarray(w)


# ---------------------------------------------------------------------------------------------
# The model: context set-up and one evaluation of a batch
# ---------------------------------------------------------------------------------------------

class BatchModel:
    def __init__(self, lib, smooth=0.0):
        t0 = time.perf_counter()
        self.ctx = fsps_batch.Context(lib, "mist", "miles", "DL07")
        for k, v in INT_SETTINGS.items():
            self.ctx.set_int(k, v)
        for k, v in FLOAT_SETTINGS.items():
            self.ctx.set_float(k, v)
        if smooth > 0:
            self.ctx.set_int("smooth_velocity", 1)
            self.ctx.set_float("sigma_smooth", smooth)
        self.setup_s = time.perf_counter() - t0
        self.dims = {k: self.ctx.dim(k) for k in ("nspec", "ntfull", "nbands", "nindx", "nz", "nemline")}
        self.log_age = self.ctx.log_age_grid()
        self.wave = self.ctx.wavelengths()
        self.leff = self.ctx.filter_leff()
        self.sfh = SfhWeights(self.log_age)

    def outputs(self, n, spec=False):
        d = self.dims
        o = {"mags": np.empty((n, d["nbands"])), "indx": np.empty((n, d["nindx"])),
             "emlines": np.empty((n, d["nemline"])), "mfrac": np.empty(n), "lbol": np.empty(n),
             "mdust": np.empty(n)}
        if spec:
            o["spec"] = np.empty((n, d["nspec"]))
        return o

    def evaluate(self, p, out, timing=None):
        """Parameter vectors p (dict of arrays) -> SEDs in out, scaled to each vector's mass.
        Everything here is the per-sampler-step cost."""
        t0 = time.perf_counter()
        w = self.sfh(p["zred"], p["logsfr_ratios"])
        arrays = {k: p[k] for k in FSPS_ARRAYS if k != "dust1"}
        arrays["dust1"] = p["dust_ratio"] * p["dust2"]
        t1 = time.perf_counter()
        self.ctx.compute_csp_batch(arrays, w, out["mags"], spec=out.get("spec"), mass=out["mfrac"],
                                   lbol=out["lbol"], mdust=out["mdust"], emlines=out["emlines"],
                                   indx=out["indx"])
        t2 = time.perf_counter()
        # Outputs are per 1 Msun formed: scale to 10**logmass (mfrac stays a fraction)
        m = p["logmass"]
        out["mags"] -= 2.5 * m[:, None]
        out["lbol"] += m                       # log10 Lsun
        out["mdust"] *= 10.0 ** m
        out["emlines"] *= 10.0 ** m[:, None]
        if "spec" in out:
            out["spec"] *= 10.0 ** m[:, None]
        t3 = time.perf_counter()
        if timing is not None:
            timing["sfh_weights"] = t1 - t0
            timing["fsps_call"] = t2 - t1
            timing["mass_scaling"] = t3 - t2
        return w


# ---------------------------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------------------------

def compare(a, b, mags_window=30.0):
    """Worst differences: magnitudes and indices absolute (magnitudes: finite bands within 30 mag
    of the brightest), L_bol in dex, other arrays relative where the reference is not negligible."""
    res = {}
    ma, mb = a["mags"], b["mags"]
    ok = np.isfinite(mb) & (mb < np.nanmin(mb, axis=1, keepdims=True) + mags_window)
    res["mags_abs"] = float(np.max(np.abs(ma - mb)[ok]))
    for k in ("indx", "emlines", "mfrac", "lbol", "mdust", "spec"):
        if k not in a or k not in b:
            continue
        x, y = np.atleast_2d(a[k]), np.atleast_2d(b[k])
        if k == "lbol":
            res["lbol_abs_dex"] = float(np.max(np.abs(x - y)))
            continue
        if k == "indx":   # magnitudes or equivalent widths (Angstrom), some near zero: absolute
            res["indx_abs"] = float(np.nanmax(np.abs(x - y)))
            continue
        else:
            d = np.abs(x - y) / np.maximum(np.abs(y), 1e-6 * np.max(np.abs(y), axis=1, keepdims=True))
        res[k + "_rel"] = float(np.nanmax(d))
    return res


def one_at_a_time(model, p, idx, spec):
    """The vectors idx computed one per call with their parameters set as context settings, the
    path of a non-batched code (only logzsol, which selects the SSP bases, stays an array)."""
    out = model.outputs(len(idx), spec)
    for j, i in enumerate(idx):
        q = {k: (v[i:i + 1] if v.ndim == 1 else v[i:i + 1]) for k, v in p.items()}
        dust1 = q["dust_ratio"] * q["dust2"]
        for k in FSPS_ARRAYS:
            if k != "logzsol":
                model.ctx.set_float(k, (dust1 if k == "dust1" else q[k])[0])
        o1 = model.outputs(1, spec)
        w = model.sfh(q["zred"], q["logsfr_ratios"])
        model.ctx.compute_csp_batch({"logzsol": q["logzsol"]}, w, o1["mags"], spec=o1.get("spec"),
                                    mass=o1["mfrac"], lbol=o1["lbol"], mdust=o1["mdust"],
                                    emlines=o1["emlines"], indx=o1["indx"])
        m = q["logmass"][0]
        o1["mags"] -= 2.5 * m
        o1["lbol"] += m
        o1["mdust"] *= 10.0 ** m
        o1["emlines"] *= 10.0 ** m
        if spec:
            o1["spec"] *= 10.0 ** m
        for k in out:
            out[k][j] = o1[k][0]
    return out


# ---------------------------------------------------------------------------------------------

def describe(model):
    d = model.dims
    lines = ["FIXED FOR THE WHOLE BATCH (shared by every parameter vector)"]
    for k, v in LIBRARIES.items():
        lines.append(f"  {k:26s} {v}")
    for k, v in {**INT_SETTINGS, **FLOAT_SETTINGS}.items():
        lines.append(f"  {k:26s} {v}")
    lines.append(f"  SSP basis                  {d['nz']} metallicities x {d['ntfull']} ages x {d['nspec']} "
                 "wavelengths, built once and kept on the GPU")
    lines.append(f"  filters                    all {d['nbands']} FSPS filters (FILTER_LIST), observed frame")
    lines.append("VARIED PER VECTOR (independent draws, no two vectors share a value)")
    for k, (kind, a, b, doc) in PRIORS.items():
        lines.append(f"  {k:14s} {kind:10s} [{a:g}, {b:g}]  {doc}")
    lines.append("RETURNED TO HOST MEMORY PER VECTOR (inside the timing)")
    lines.append(f"  mags     {d['nbands']} AB magnitudes (redshifted, IGM, distance modulus, filter integrals)")
    lines.append(f"  indx     {d['nindx']} spectral indices;  emlines {d['nemline']} line luminosities")
    lines.append("  mfrac    surviving stellar mass (incl. remnants) per mass formed;  lbol, mdust")
    lines.append(f"  spec     (with --spec) {d['nspec']}-pixel rest-frame spectrum")
    return "\n".join(lines)


def gpu_others():
    """Compute processes on the GPU other than this one (a shared node slows the kernels)."""
    try:
        q = subprocess.run(["nvidia-smi", "--query-compute-apps=pid,used_memory", "--format=csv,noheader"],
                           capture_output=True, text=True, timeout=20).stdout
        return [l.strip() for l in q.splitlines() if l.strip() and l.split(",")[0].strip() != str(os.getpid())]
    except Exception:
        return []


def gpu_name():
    try:
        return subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"],
                              capture_output=True, text=True, timeout=20).stdout.strip()
    except Exception:
        return "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lib", required=True, help="path to libfsps.so (e.g. an -Daccel=gpu build)")
    ap.add_argument("--sizes", default="1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000",
                    help="batch sizes to time")
    ap.add_argument("--reps", type=int, default=20, help="timed calls per size (fresh draws each)")
    ap.add_argument("--spec", action="store_true", help="also time with spectra copied to the host")
    ap.add_argument("--smooth", type=float, default=0.0, help="velocity dispersion [km/s] (0: none)")
    ap.add_argument("--check", type=int, default=24, help="vectors recomputed one at a time")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="out")
    ap.add_argument("--replay", help="npz written by a run: recompute its vectors and compare")
    ap.add_argument("--describe", action="store_true", help="print the setup and exit")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    model = BatchModel(args.lib, args.smooth)
    if args.describe:
        print(describe(model))
        return

    if args.replay:
        z = np.load(args.replay)
        p = {k[2:]: z[k] for k in z.files if k.startswith("p_")}
        ref = {k[2:]: z[k] for k in z.files if k.startswith("o_")}
        out = model.outputs(len(p["zred"]), "spec" in ref)
        model.evaluate(p, out)
        res = compare(out, ref)
        print(f"replay of {len(p['zred'])} vectors with {args.lib}: worst differences vs the saved run")
        for k, v in res.items():
            print(f"  {k:16s} {v:.2e}")
        with open(os.path.join(args.out, "replay.json"), "w") as f:
            json.dump({"lib": args.lib, "n": int(len(p["zred"])), "differences": res}, f, indent=1)
        return

    print(describe(model))
    print(f"\nGPU: {gpu_name()};  host: {platform.node()};  library: {args.lib}")
    rng = np.random.default_rng(args.seed)
    sizes = [int(s) for s in args.sizes.split(",")]
    result = {"gpu": gpu_name(), "lib": args.lib, "smooth_kms": args.smooth, "dims": model.dims,
              "setup_s": model.setup_s, "timing": []}

    # One-time cost: the first call builds the SSP bases of all metallicities and moves the
    # context to the GPU (a sampler pays this once per run)
    p = draw_vectors(rng, 16)
    t0 = time.perf_counter()
    model.evaluate(p, model.outputs(16))
    result["first_call_s"] = time.perf_counter() - t0
    print(f"\nOne-time: setup {model.setup_s:.1f} s; first call (12 SSP bases, context to GPU) "
          f"{result['first_call_s']:.1f} s")

    # Timing: fresh, independent vectors every call
    print("\nTiming (wall clock in Python, per call: SFH weights + FSPS batch + mass scaling; "
          "outputs in numpy arrays)")
    print(f"{'N':>7s} {'output':>12s} {'us/vector median':>17s} {'min':>8s}   "
          f"{'SFH weights':>11s} {'FSPS call':>10s} {'scaling':>8s}  (us/vector, median)")
    for spec in ([False, True] if args.spec else [False]):
        for n in sizes:
            out = model.outputs(n, spec)
            model.evaluate(draw_vectors(rng, n), out)      # warm-up at this size
            reps = max(3, min(args.reps, int(3e5 // n)))
            tot, parts = [], {"sfh_weights": [], "fsps_call": [], "mass_scaling": []}
            for _ in range(reps):
                p = draw_vectors(rng, n)                   # the sampler's proposals (not timed)
                tm = {}
                t0 = time.perf_counter()
                model.evaluate(p, out, tm)
                tot.append(time.perf_counter() - t0)
                for k in parts:
                    parts[k].append(tm[k])
            others = gpu_others()
            us = 1e6 * np.array(tot) / n
            med = {k: 1e6 * float(np.median(v)) / n for k, v in parts.items()}
            label = "+spectra" if spec else "phot+props"
            print(f"{n:7d} {label:>12s} {np.median(us):17.2f} {us.min():8.2f}   {med['sfh_weights']:11.2f} "
                  f"{med['fsps_call']:10.2f} {med['mass_scaling']:8.2f}"
                  + (f"   (GPU shared with {len(others)} other process(es))" if others else ""))
            result["timing"].append({"n": n, "spec": spec, "reps": reps, "us_median": float(np.median(us)),
                                     "us_min": float(us.min()), "parts_us_median": med,
                                     "other_gpu_processes": others})
            # the batch really varies every parameter per vector
            if n == min(max(sizes), 10000) and not spec:
                result["n_unique"] = {k: int(np.unique(v).size) for k, v in p.items() if v.ndim == 1}

    # Per-vector correctness: vectors of a big batch against the same vectors one per call
    n = min(max(sizes), 10000)
    p = draw_vectors(rng, n)
    out = model.outputs(n, True)
    w = model.evaluate(p, out)
    idx = np.sort(rng.choice(n, size=min(args.check, n), replace=False))
    single = one_at_a_time(model, p, idx, True)
    sub = {k: v[idx] for k, v in out.items()}
    res = compare(sub, single)
    result["check_batch_vs_single"] = res
    print(f"\nCheck: {len(idx)} vectors of the N = {n} batch vs the same vectors one per call, "
          "parameters set as context scalars (worst differences):")
    for k, v in res.items():
        print(f"  {k:16s} {v:.2e}")
    spread = {k: (float(np.nanmin(sub["mags"][:, b])), float(np.nanmax(sub["mags"][:, b])))
              for k, b in (("SDSS r", 14), ("WISE W3", 84))}
    print("  (outputs differ between vectors as they should: e.g. "
          + ", ".join(f"{k} {a:.1f} ... {b:.1f} mag" for k, (a, b) in spread.items()) + ")")

    np.savez(os.path.join(args.out, "vectors.npz"), **{"p_" + k: v[idx] for k, v in p.items()},
             **{"o_" + k: v for k, v in sub.items()})
    with open(os.path.join(args.out, "results.json"), "w") as f:
        json.dump(result, f, indent=1)
    try:
        import plots
        plots.make(model, p, out, w, result, args.out)
        print(f"\nFigures in {args.out}/")
    except ImportError as e:
        print(f"(no figures: {e})")


if __name__ == "__main__":
    sys.exit(main())
