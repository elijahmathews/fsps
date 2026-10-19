"""Figures for demo_batch_sed.py: time per vector vs batch size, and example SEDs of one batch."""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Categorical order (blue, orange, aqua, yellow, magenta), fixed; text stays in ink colors
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK, MUTED, GRID = "#1f2328", "#59636e", "#d8dee4"

# A broad, typical photometric set (FILTER_LIST numbers, 1-based)
SHOW_BANDS = [62, 63, 13, 14, 15, 16, 17, 10, 11, 12, 115, 117, 118, 119, 120, 121, 122,
              83, 84, 85, 86, 90, 95, 96, 97, 98, 99, 100]


def _style(ax):
    ax.grid(True, which="major", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)


def timing_figure(result, path, throughput=None):
    rows = result["timing"]
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=150)
    series = [("end to end (Python: SFH weights + FSPS batch + mass scaling)", False, "us_median", 0),
              ("FSPS batch call alone", False, "fsps", 1),
              ("end to end, rest-frame spectra also copied to numpy", True, "us_median", 2)]
    for label, spec, key, c in series:
        r = [x for x in rows if x["spec"] == spec]
        if not r:
            continue
        n = np.array([x["n"] for x in r])
        y = np.array([x["parts_us_median"]["fsps_call"] if key == "fsps" else x[key] for x in r])
        ax.plot(n, y, "-o", color=COLORS[c], lw=2, ms=4.5, label=label)
        dy = {0: 5, 1: -6}.get(c, 0)
        ax.annotate(f"{y[-1]:.1f} µs", (n[-1], y[-1]), xytext=(7, dy), textcoords="offset points",
                    va="center", color=INK, fontsize=8)
    if throughput:
        n2, y2 = throughput
        ax.plot([n2], [y2], "D", color=COLORS[3], ms=7, mec="white", mew=1,
                label=f"2 processes sharing the GPU, aggregate (N = {n2} each)")
        ax.annotate(f"{y2:.1f} µs", (n2, y2), xytext=(0, -11), textcoords="offset points", ha="center",
                    color=INK, fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("batch size N (parameter vectors per call)", color=INK)
    ax.set_ylabel("wall time per parameter vector [µs]", color=INK)
    ax.set_title(f"Batched CSP + photometry, {result['gpu'].split(',')[0]}", color=INK, fontsize=11, loc="left")
    _style(ax)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def sed_figure(model, p, out, path, k=5):
    # vectors spread in redshift so the panels show the per-vector observed frame
    order = np.argsort(p["zred"])
    idx = order[(np.linspace(0.05, 0.95, k) * (len(order) - 1)).astype(int)]
    lam = model.wave
    nu = 2.998e18 / lam
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.0), dpi=150)
    bands = np.array(SHOW_BANDS) - 1
    for c, i in enumerate(idx):
        lab = (f"z={p['zred'][i]:.2f}, logM={p['logmass'][i]:.1f}, [Z]={p['logzsol'][i]:+.2f}, "
               f"τ₂={p['dust2'][i]:.2f}, f_AGN={p['fagn'][i]:.0e}")
        a1.loglog(lam / 1e4, nu * out["spec"][i], color=COLORS[c], lw=1.0, label=lab)
        m = out["mags"][i, bands]
        ok = np.isfinite(m) & (m < 40)
        a2.semilogx(model.leff[bands][ok] / 1e4, m[ok], "o", color=COLORS[c], ms=5, mec="white", mew=0.8)
    a1.set_xlim(0.09, 1000)
    top = np.nanmax(nu[None, :] * out["spec"][idx])
    a1.set_ylim(top * 1e-6, top * 3)
    a1.set_xlabel("rest-frame wavelength [µm]", color=INK)
    a1.set_ylabel(r"$\nu L_\nu$ [L$_\odot$]", color=INK)
    a1.set_title("Rest-frame spectra (5 vectors of one batch)", color=INK, fontsize=10, loc="left")
    a2.invert_yaxis()
    a2.set_xlabel("observed filter effective wavelength [µm]", color=INK)
    a2.set_ylabel("apparent AB magnitude", color=INK)
    a2.set_title("Their photometry (28 of the 159 filters: GALEX to SPIRE)", color=INK,
                 fontsize=10, loc="left")
    for a in (a1, a2):
        _style(a)
    fig.legend(*a1.get_legend_handles_labels(), frameon=False, fontsize=8, labelcolor=INK,
               loc="lower center", ncol=3)
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(path)
    plt.close(fig)
    return idx


def _throughput(outdir):
    """(N, us per vector) of the 2-process line of throughput_processes.txt, if it was run."""
    try:
        for line in open(os.path.join(outdir, "throughput_processes.txt")):
            if line.startswith("2 process"):
                return int(line.split("N = ")[1].split()[0]), float(line.split(":")[1].split()[0])
    except (OSError, IndexError, ValueError):
        pass
    return None


def make(model, p, out, w, result, outdir):
    timing_figure(result, os.path.join(outdir, "timing_vs_batch_size.png"), _throughput(outdir))
    if "spec" in out:
        sed_figure(model, p, out, os.path.join(outdir, "example_seds.png"))
