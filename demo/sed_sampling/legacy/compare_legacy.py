#!/usr/bin/env python3
"""The demo's parameter vectors through this fork's batch and through legacy (upstream) FSPS:
speed and numerical differences on identical inputs. Driven by run_legacy_compare.sh.

  compare_legacy.py fork --lib LIB --tag TAG --dir DIR [--n N] [--seed S]
      (first call) draws N vectors from the demo's priors, writes DIR/inputs.npz, the legacy
      driver's DIR/vectors.bin and DIR/settings.txt; then evaluates them with the fork library
      LIB and writes DIR/fork_TAG.npz (outputs and timing)
  compare_legacy.py verify-data --legacy SPS_HOME --ref REV --data DATA --dir DIR
      checks that every file of the fork's data tree (as the legacy layout names it) is the
      legacy revision's own file; writes DIR/data_provenance.txt
  compare_legacy.py report --dir DIR
      compares DIR/fork_*.npz with DIR/legacy_*.bin, writes DIR/report.md and figures

Every vector is evaluated per 1 Msun formed (logmass is a caller-side scaling, identical in both
codes, so it is left out); the SSP weights are computed once here and given to both codes.
"""

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import demo_batch_sed as demo  # noqa: E402

# Settings both codes get (legacy: settings.txt before SPS_SETUP; fork: the context), the
# demo's plus the defaults that shape the SSPs, set explicitly on both sides
EXTRA_INT = {"add_agb_dust_model": 1, "use_wr_spectra": 1, "tpagb_norm_type": 2}
EXTRA_FLOAT = {"dust1_index": -1.0}

# Index windows of FILTER_LIST (1-based numbers) for the band groups of the report
BAND_GROUPS = [("UV (GALEX)", [62, 63]), ("optical (SDSS ugriz)", [13, 14, 15, 16, 17]),
               ("near-IR (2MASS JHKs)", [10, 11, 12]), ("JWST NIRCam", list(range(115, 123))),
               ("mid-IR (WISE, IRAC)", [83, 84, 85, 86, 53, 54, 55, 56]),
               ("far-IR (MIPS, PACS, SPIRE)", [90, 91, 92, 95, 96, 97, 98, 99, 100])]


def filter_names():
    names = []
    with open(os.path.join(os.environ.get("FSPS_DATA_HOME", "."), "data", "FILTER_LIST")) as f:
        for line in f:
            p = line.split()
            if len(p) >= 2 and p[0].isdigit():
                names.append(p[1])
    return names


def _legacy_name(f):
    """A fork data file's path in the legacy SPS_HOME layout."""
    top, _, rest = f.partition("/")
    return {"isochrones": "ISOCHRONES/" + rest, "spectra": "SPECTRA/" + rest}.get(
        top, f if top in ("dust", "nebular") else "data/" + f)


def _git_blob_sha(path):
    import hashlib
    h = hashlib.sha1()
    h.update(b"blob %d\0" % os.path.getsize(path))
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_data(args):
    import subprocess
    fork = {}
    for dp, _, fs in os.walk(args.data):
        for f in fs:
            rel = os.path.relpath(os.path.join(dp, f), args.data)
            fork[_legacy_name(rel)] = _git_blob_sha(os.path.join(dp, f))
    is_git = subprocess.run(["git", "-C", args.legacy, "rev-parse"], capture_output=True).returncode == 0

    def blobs(ref):
        out = subprocess.run(["git", "-C", args.legacy, "ls-tree", "-r", ref], capture_output=True,
                             text=True, check=True).stdout
        return {l.split("\t")[1]: l.split()[2] for l in out.splitlines()}

    if is_git:
        leg, head = blobs(args.ref), blobs("HEAD")
    else:
        leg = {k: _git_blob_sha(os.path.join(args.legacy, k)) for k in fork
               if os.path.exists(os.path.join(args.legacy, k))}
        head = {}
    same = sorted(k for k in fork if leg.get(k) == fork[k])
    diff = sorted(k for k in fork if k in leg and leg[k] != fork[k])
    absent = sorted(k for k in fork if k not in leg)
    lines = [f"Fork data tree {args.data} against legacy {args.legacy} at {args.ref}",
             f"  identical to the legacy revision's file: {len(same)} of {len(fork)}"]
    for k in diff:
        note = " (= legacy HEAD's file)" if head.get(k) == fork[k] else ""
        lines.append(f"  DIFFERENT: {k}{note}")
    for k in absent:
        lines.append(f"  not in legacy: {k}")
    txt = "\n".join(lines) + "\n"
    with open(os.path.join(args.dir, "data_provenance.txt"), "w") as fh:
        fh.write(txt)
    print(txt, end="")


def prepare(model, n, seed, d):
    rng = np.random.default_rng(seed)
    p = demo.draw_vectors(rng, n)
    p["logmass"][:] = 0.0
    w = model.sfh(p["zred"], p["logsfr_ratios"])
    np.savez(os.path.join(d, "inputs.npz"), weights=w, log_age=model.log_age, **{"p_" + k: v for k, v in p.items()})
    par = np.stack([p["dust_ratio"] * p["dust2"] if k == "dust1" else p[k] for k in demo.FSPS_ARRAYS])
    with open(os.path.join(d, "vectors.bin"), "wb") as f:
        np.array([n, w.shape[1]], dtype=np.int32).tofile(f)
        # tofile writes C order: (n, 12) rows are Fortran's params(12, n) columns
        np.ascontiguousarray(par.T, dtype=np.float64).tofile(f)
        np.ascontiguousarray(w, dtype=np.float64).tofile(f)        # weights(ntfull, n)
    with open(os.path.join(d, "settings.txt"), "w") as f:
        for k, v in {**demo.INT_SETTINGS, **EXTRA_INT}.items():
            f.write(f"{k}={v}\n")
        for k, v in {**demo.FLOAT_SETTINGS, **EXTRA_FLOAT}.items():
            f.write(f"{k}={v!r}\n")


def run_fork(args):
    model = demo.BatchModel(args.lib, extra_ints=EXTRA_INT, extra_floats=EXTRA_FLOAT)
    if not os.path.exists(os.path.join(args.dir, "inputs.npz")):
        prepare(model, args.n, args.seed, args.dir)
    z = np.load(os.path.join(args.dir, "inputs.npz"))
    if not np.array_equal(z["log_age"], model.log_age):
        raise SystemExit("SSP age grid differs from the one the inputs were made for")
    p = {k[2:]: z[k] for k in z.files if k.startswith("p_")}
    w = z["weights"]
    n = w.shape[0]
    arrays = {k: p[k] for k in demo.FSPS_ARRAYS if k != "dust1"}
    arrays["dust1"] = p["dust_ratio"] * p["dust2"]
    out = model.outputs(n, spec=True)

    def call():
        model.ctx.compute_csp_batch(arrays, w, out["mags"], spec=out["spec"], mass=out["mfrac"],
                                    lbol=out["lbol"], mdust=out["mdust"], emlines=out["emlines"],
                                    indx=out["indx"])

    t0 = time.perf_counter()
    call()                                   # first call: the 12 SSP bases (and the device copy)
    first = time.perf_counter() - t0
    reps = []
    for _ in range(args.reps):
        t0 = time.perf_counter()
        call()
        reps.append(time.perf_counter() - t0)
    np.savez(os.path.join(args.dir, f"fork_{args.tag}.npz"), lib=args.lib, setup_s=model.setup_s,
             first_call_s=first, call_s=np.array(reps), wave=model.wave, **out)
    print(f"fork {args.tag}: {n} vectors, {1e6 * np.median(reps) / n:.1f} us per vector "
          f"(one call, spectra included)")


def read_legacy(path):
    with open(path, "rb") as f:
        n, nspec, nb, ni, ne, nt = np.fromfile(f, np.int32, 6)
        tf = np.fromfile(f, np.float64, nt)
        lam = np.fromfile(f, np.float64, nspec)

        def arr(k, m=None):
            a = np.fromfile(f, np.float64, k * n)
            return a.reshape(n, k) if m is None else a

        o = {"mags": arr(nb), "spec": arr(nspec), "indx": arr(ni), "emlines": arr(ne)}
        o["mfrac"], o["lbol"], o["mdust"] = (np.fromfile(f, np.float64, n) for _ in range(3))
        t = np.fromfile(f, np.float64, 5)
    o["time_full"], o["wave"] = tf, lam
    o["timing"] = dict(zip(["setup_s", "ssp_all_z_s", "ztinterp_s", "compsp_s", "getindx_s"], t))
    return o


# ---------------------------------------------------------------------------------------------

def stats(x):
    x = np.asarray(x)[np.isfinite(x)]
    if x.size == 0:
        return (np.nan, np.nan, np.nan)
    return (float(np.median(x)), float(np.percentile(x, 99)), float(np.max(x)))


def differences(a, b, wave, names):
    """a vs reference b: per-quantity (median, 99th percentile, max) over vectors (and pixels,
    bands, indices, lines)."""
    r = {}
    # bands with flux in both (legacy writes 99 for none), within 30 mag of the vector's
    # brightest band (per 1 Msun formed at z ~ 1-5 magnitudes are ~45-55)
    ma, mb = a["mags"], b["mags"]
    ok = np.isfinite(ma) & np.isfinite(mb) & (mb < 90) & (ma < 90)
    ok &= mb < np.min(np.where(ok, mb, np.inf), axis=1, keepdims=True) + 30.0
    dm = np.where(ok, np.abs(ma - mb), np.nan)
    r["n_mags_compared"] = int(ok.sum())
    r["n_mags_total"] = int(ok.size)
    r["mags (all bands) [mag]"] = stats(dm)
    r["per_band_max"] = {names[i] if i < len(names) else str(i + 1): float(np.nanmax(dm[:, i]))
                         for i in range(dm.shape[1]) if np.isfinite(dm[:, i]).any()}
    for label, idx in BAND_GROUPS:
        r[f"mags {label} [mag]"] = stats(dm[:, np.array(idx) - 1])
    sa, sb = a["spec"], b["spec"]
    peak = np.max(sb, axis=1, keepdims=True)
    rel = np.abs(sa / np.where(sb > 0, sb, np.nan) - 1.0)
    sig = sb > 1e-6 * peak
    for label, lo, hi in [("912 A - 2.5 um", 912, 2.5e4), ("2.5 - 1000 um (dust, AGN)", 2.5e4, 1e7),
                          ("< 912 A (ionizing)", 0, 912)]:
        sel = sig & (wave[None, :] >= lo) & (wave[None, :] < hi)
        r[f"spectrum {label} [rel]"] = stats(rel[sel])
    r["_rel_spec"] = np.where(sig, rel, np.nan)
    r["indices [A or mag]"] = stats(np.abs(a["indx"] - b["indx"]))
    lb = b["emlines"]
    lsel = lb > 1e-6 * np.max(lb, axis=1, keepdims=True)
    r["emission lines [rel]"] = stats(np.abs(a["emlines"] / np.where(lsel, lb, np.nan) - 1)[lsel])
    r["mfrac [rel]"] = stats(np.abs(a["mfrac"] / b["mfrac"] - 1))
    r["log L_bol [dex]"] = stats(np.abs(a["lbol"] - b["lbol"]))
    r["dust mass [rel]"] = stats(np.abs(a["mdust"] / b["mdust"] - 1))
    return r


def report(args):
    d = args.dir
    names = filter_names()
    forks = {t: dict(np.load(os.path.join(d, f"fork_{t}.npz"))) for t in ("gpu", "cpu")
             if os.path.exists(os.path.join(d, f"fork_{t}.npz"))}
    legs = {t: read_legacy(os.path.join(d, f"legacy_{t}.bin")) for t in ("aligned_rv", "aligned", "stock")
            if os.path.exists(os.path.join(d, f"legacy_{t}.bin"))}
    inp = np.load(os.path.join(d, "inputs.npz"))
    n = inp["weights"].shape[0]
    wave = next(iter(legs.values()))["wave"]
    for f in forks.values():
        if not np.allclose(f["wave"], wave, rtol=1e-12):
            raise SystemExit("wavelength grids differ")
    for leg in legs.values():
        if not np.allclose(leg["time_full"], inp["log_age"], atol=1e-12):
            raise SystemExit("SSP age grids differ")

    comps = []
    for ft in forks:
        for lt in legs:
            comps.append((f"fork {ft.upper()} vs legacy {lt}", differences(forks[ft], legs[lt], wave, names)))
    if "gpu" in forks and "cpu" in forks:
        comps.append(("fork GPU (FP32) vs fork CPU (FP64)", differences(forks["gpu"], forks["cpu"], wave, names)))

    # Timing (per vector, all outputs incl. spectra on both sides)
    tim = {}
    for lt, leg in legs.items():
        t = leg["timing"]
        tim[f"legacy {lt}"] = 1e6 * (t["ztinterp_s"] + t["compsp_s"] + t["getindx_s"]) / n
    for ft, f in forks.items():
        tim[f"fork {ft}"] = 1e6 * float(np.median(f["call_s"])) / n

    with open(os.path.join(d, "report.md"), "w") as fo:
        fo.write(f"# Fork vs legacy FSPS on {n} demo parameter vectors\n\n")
        fo.write("## Time per vector (all outputs incl. spectra, per 1 Msun formed)\n\n")
        fo.write("| Code | µs per vector | vs legacy |\n|---|---|---|\n")
        ref = tim.get("legacy aligned", tim.get("legacy stock"))
        for k, v in tim.items():
            fo.write(f"| {k} | {v:,.1f} | {ref / v:,.0f}x |\n")
        for lt, leg in legs.items():
            t = leg["timing"]
            fo.write(f"\nLegacy {lt}: setup {t['setup_s']:.1f} s, SSPs of all metallicities "
                     f"{t['ssp_all_z_s']:.1f} s (once); per vector ztinterp "
                     f"{1e3 * t['ztinterp_s'] / n:.2f} ms, compsp {1e3 * t['compsp_s'] / n:.2f} ms, "
                     f"getindx {1e3 * t['getindx_s'] / n:.3f} ms.\n")
        for ft, f in forks.items():
            fo.write(f"\nFork {ft}: setup {float(f['setup_s']):.1f} s, first call (SSP bases) "
                     f"{float(f['first_call_s']):.1f} s, then one call of {n} vectors "
                     f"{1e3 * float(np.median(f['call_s'])):.1f} ms.\n")
        fo.write("\n## Differences (median / 99th percentile / max over vectors)\n")
        for title, r in comps:
            fo.write(f"\n### {title}\n\n| Quantity | median | p99 | max |\n|---|---|---|---|\n")
            for k, v in r.items():
                if k.startswith("_") or k in ("per_band_max", "n_mags_compared", "n_mags_total"):
                    continue
                fo.write(f"| {k} | {v[0]:.2e} | {v[1]:.2e} | {v[2]:.2e} |\n")
            worst = sorted(r["per_band_max"].items(), key=lambda kv: -kv[1])[:8]
            fo.write(f"\nMagnitudes compared: {r['n_mags_compared']} of {r['n_mags_total']} "
                     "(with flux in both, within 30 mag of the vector's brightest band). Largest per-band max: "
                     + ", ".join(f"{k} {v:.4f}" for k, v in worst) + "\n")
    summary = {"n": n, "time_us_per_vector": tim,
               "differences": {t: {k: v for k, v in r.items() if not k.startswith("_")} for t, r in comps}}
    with open(os.path.join(d, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    figure(comps, wave, os.path.join(d, "legacy_differences.png"))
    print(open(os.path.join(d, "report.md")).read())


def figure(comps, wave, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plots import COLORS, INK, _style
    fig, ax = plt.subplots(figsize=(9, 4.6), dpi=150)
    shown = [(t, r) for t, r in comps if t.startswith("fork GPU")]   # GPU vs each legacy build, vs CPU
    labels = {"legacy stock": "legacy stock (upstream as is)",
              "legacy aligned_rv": "legacy + constants, AGN fix, Calzetti R_V",
              "legacy aligned": "legacy + constants, AGN fix",
              "fork CPU (FP64)": "fork CPU FP64 (precision of the GPU kernels)"}
    for c, (title, r) in enumerate(shown[:5]):
        label = next((v for k, v in labels.items() if title.endswith(k)), title)
        rel = r["_rel_spec"]
        with np.errstate(all="ignore"):
            med = np.nanmedian(rel, axis=0)
        ax.loglog(wave / 1e4, med, color=COLORS[c], lw=1.2, label="fork GPU vs " + label)
    ax.set_xlim(0.05, 1000)
    ax.set_xlabel("rest-frame wavelength [µm]", color=INK)
    ax.set_ylabel("|fork / reference - 1|, median over vectors", color=INK)
    ax.set_title(f"Spectra of the same parameter vectors: per-pixel relative difference", color=INK,
                 fontsize=10, loc="left")
    _style(ax)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["fork", "report", "verify-data"])
    ap.add_argument("--legacy")
    ap.add_argument("--ref")
    ap.add_argument("--data")
    ap.add_argument("--lib")
    ap.add_argument("--tag")
    ap.add_argument("--dir", required=True)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--reps", type=int, default=3)
    args = ap.parse_args()
    os.makedirs(args.dir, exist_ok=True)
    if args.mode == "fork":
        run_fork(args)
    elif args.mode == "verify-data":
        verify_data(args)
    else:
        report(args)


if __name__ == "__main__":
    main()
