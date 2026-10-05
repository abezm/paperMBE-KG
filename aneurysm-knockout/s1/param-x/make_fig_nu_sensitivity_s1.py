#!/usr/bin/env python3
"""Sensitivity of the equilibrated (MBE-KG) response to the transition decay nu of the
degradation profile, Eq. (34), at fixed axial extent s_z = 3 mm (s1 copy, 3 Oct 2026).

nu and s_z are not calibrated against data, so the dependence of the predictions on them
has to be reported rather than assumed. Two quantities are shown:

  (A) the relative volume J at a FIXED nominal lesion (the grid point nearest to a 20%
      loss of elastic-fiber integrity and a 10% loss of collagen cross-linking), and
  (B) the size of the admissible region, as the fraction of the 20x20 (d_e, d_c) sweep
      that yields a bounded equilibrated solution.

Admissible is the paper's rule (Section 5.5): converged and apex lambda_r/lambda_ro <= 2.5.

Reads MASTER_PARAMETRIC_STUDY.csv of s1/param-x. Writes nu-sensitivity-s1.pdf next to main.tex.
"""
import csv
import numpy as np
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

HERE = os.path.dirname(os.path.abspath(__file__))
TEX_DIR = os.path.normpath(os.path.join(HERE, "..", "..", "..", "paper"))
GRID = 400          # 20 x 20 sweep per configuration
LAM_MAX = 2.5      # admissible: converged and apex lambda_r/lambda_ro <= 2.5
SZ = "3.0"          # the axial extent of the nu sweep (author, 3 Oct 2026)
NUS = [1, 2, 3, 4, 5, 6]

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.8,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})

# The MASTER csv also holds runs of OTHER materials (the "unit-frac" folder carries unit
# mass fractions). Filtering on type/nu/zod/alph alone swallows them and corrupts nu=2,
# so the folder must be filtered too: this figure is one material only.
rows = [r for r in csv.DictReader(open(os.path.join(HERE, "MASTER_PARAMETRIC_STUDY.csv")))
        if r["alph"] == "ON" and r["z4_zod"] == SZ and "unit-frac" not in r["folder"]]

seen = {(r["type"], r["z3_nu"], r["x_elastin_deg"], r["y_collagen_deg"]) for r in rows}
assert len(seen) == len(rows), "duplicate (d_e,d_c) points would be double-counted"

xs = sorted({float(r["x_elastin_deg"]) for r in rows})
ys = sorted({float(r["y_collagen_deg"]) for r in rows})
DE = min(xs, key=lambda v: abs(v - 0.20))   # nominal lesion: ~20% elastin loss
DC = min(ys, key=lambda v: abs(v - 0.10))   #                 ~10% collagen cross-link loss

# both curves are MBE-KG results of the same model: one colour, told apart by line and marker
STYLE = {"axis": ("fusiform", "o", "-", "black"),
         "asym": ("saccular", "s", "--", "black")}

def tidy_axis(a, lo, hi, which="y", nbins=5):
    """First and last major tick exactly on the frame (house rule of Figs. 2, 3, 6, 7)."""
    from matplotlib.ticker import MaxNLocator, FormatStrFormatter
    t = MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
    step = t[1] - t[0]
    d = next(k for k in range(6) if abs(round(step, k) - step) < 1e-9 * max(1.0, abs(step)))
    ax_ = a.yaxis if which == "y" else a.xaxis
    (a.set_yticks if which == "y" else a.set_xticks)(t)
    (a.set_ylim if which == "y" else a.set_xlim)(t[0], t[-1])
    ax_.set_major_formatter(FormatStrFormatter(f"%.{d}f"))


fig, ax = plt.subplots(1, 2, figsize=(6.6, 2.9))
ylo, yhi = [np.inf, np.inf], [-np.inf, -np.inf]
for t, (label, mk, ls, col) in STYLE.items():
    Jnom, frac = [], []
    for nu in NUS:
        g = [r for r in rows if r["type"] == t and r["z3_nu"] == f"{nu}.0"]
        adm = [r for r in g if float(r["z1_lamb_r"]) <= LAM_MAX]
        hit = [r for r in g
               if abs(float(r["x_elastin_deg"]) - DE) < 1e-9
               and abs(float(r["y_collagen_deg"]) - DC) < 1e-9]
        Jnom.append(float(hit[0]["z2_J"]) if hit else float("nan"))
        frac.append(100.0 * len(adm) / GRID)
    ax[0].plot(NUS, Jnom, ls=ls, marker=mk, ms=4, mfc="none", color=col, lw=1.1, label=label, clip_on=False)
    ax[1].plot(NUS, frac, ls=ls, marker=mk, ms=4, mfc="none", color=col, lw=1.1, clip_on=False)
    for k, v in enumerate((Jnom, frac)):
        ylo[k], yhi[k] = min(ylo[k], np.nanmin(v)), max(yhi[k], np.nanmax(v))
    print(f"{t:5s} J(nominal) = {[round(v,3) for v in Jnom]}   admissible % = {[round(v,1) for v in frac]}")

ax[0].set_ylabel(r"$J$ [-]")
ax[1].set_ylabel("admissible states [%]")
for i, a in enumerate(ax):
    a.set_xlabel(r"$\nu$ [-]")
    a.set_xlim(NUS[0], NUS[-1])                  # the axis ends on labelled ticks
    a.set_xticks(NUS)
    tidy_axis(a, 0.0 if i == 1 else ylo[i], yhi[i])   # (b) is a percentage: from 0 % (author, 4 Oct 2026)
    a.yaxis.set_minor_locator(AutoMinorLocator(2))
    a.tick_params(which="minor", length=2)
    # (b) top right: top left holds the saccular point at nu = 1 (58.6 %)
    a.text(0.04 if i == 0 else 0.96, 0.90, f"({'ab'[i]})", transform=a.transAxes,
           ha="left" if i == 0 else "right")
ax[0].legend(fontsize=8, frameon=False, loc="center right")
fig.tight_layout()
out = os.path.join(TEX_DIR, "nu-sensitivity-s1.pdf")
if os.path.exists(out):
    raise SystemExit(f"refusing to overwrite {out}")
fig.savefig(out)
print("\nnominal lesion used: d_e = %.2f, d_c = %.2f (s_z = %s mm)" % (DE, DC, SZ))
print("wrote", out)
