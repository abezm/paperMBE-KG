#!/usr/bin/env python3
"""Figure 2: uniaxial tension, full-KG (rate-dependent) against MBE-KG (equilibrated).

Same layout as the aneurysm figures (make_fig_aneurysm_4fib.py): three panels side by
side against the normalized G&R time s/s_load,

  (A) axial Cauchy stress [kPa]
  (B) principal stretch ratio lambda/lambda_o
  (C) growth stretch lambda_g = J^(1/3)

Lines  = finite element (FEBio)    (prep_uniax_consistent.py   -> <case>_c.txt)
Symbols= semi-analytical solution  (semi-analytical/uniax_4fib.py       -> <case>_a.txt)
Red = MBE-KG (equilibrated).  Black = full-KG at s_load = 1, 10, 100 days.

The pull is along the material axial direction: the specimen sits off the global z axis,
so the fibre basis is uniform over it and alpha_0 is measured from the pull direction.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, FormatStrFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
TEX = os.path.normpath(os.path.join(HERE, "..", "..", "paper"))

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})

RED = "red"
CASES = [("sp1",   "black", r"$s_{\mathrm{load}} = 1$ day"),
         ("sp10",  "black", r"10 days"),
         ("sp100", "black", r"100 days"),
         ("MBE",   RED,     "MBE-KG")]

KEYS = [1, 2, 3]
YLAB = [r"$\sigma$ [kPa]", r"$\lambda/\lambda_{o}$ [-]",
        r"$\lambda_{\mathrm{g}} = J^{1/3}$ [-]"]
# per panel: (case, anchor x on its curve, label x, label y) in data coordinates
ANN = {1: [("sp1", 3.5, 3.8, 690), ("sp10", 2.1, 3.0, 430), ("sp100", 0.5, 1.4, 225)],
       2: [("sp1", 3.5, 4.0, 1.08), ("sp10", 3.5, 4.1, 1.45), ("sp100", 0.75, 2.6, 1.33)],
       3: [("sp1", 3.5, 4.0, 1.06), ("sp10", 3.5, 4.1, 1.45), ("sp100", 0.75, 2.6, 1.33)]}


def tidy_yaxis(a, lo, hi, nbins=5):
    """Major ticks on round values, with the first and last one exactly on the frame,
    so the extreme labels sit at the corners of the box. The number of decimals is the
    least that keeps two neighbouring ticks distinct."""
    t = MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
    step = t[1] - t[0]
    d = next(k for k in range(6) if abs(round(step, k) - step) < 1e-9 * max(1.0, abs(step)))
    a.set_yticks(t)
    a.set_ylim(t[0], t[-1])
    a.yaxis.set_major_formatter(FormatStrFormatter(f"%.{d}f"))
    a.yaxis.set_minor_locator(AutoMinorLocator(2))


ana = {t: np.loadtxt(os.path.join(HERE, f"{t}_a.txt")) for t, _, _ in CASES}
fem = {t: np.loadtxt(os.path.join(HERE, f"{t}_c.txt")) for t, _, _ in CASES}
fem = {t: c[c[:, 0] >= 0] for t, c in fem.items()}

fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.9))

for i, k in enumerate(KEYS):
    a = ax[i]
    lo, hi = np.inf, -np.inf
    for tag, col, _ in CASES:
        z = 3 if tag == "MBE" else 2
        s, c = ana[tag], fem[tag]
        a.plot(c[:, 0], c[:, k], color=col, lw=1.0, zorder=z)        # FE: solid line
        # semi-analytical symbols at the same spacing in s/s_load (0.2) for every curve
        grid = np.arange(0.0, s[-1, 0] + 1e-9, 0.2)
        idx = np.unique([np.argmin(abs(s[:, 0] - g)) for g in grid])
        a.plot(s[idx, 0], s[idx, k], ls="none", marker="o", ms=3.2,
               mfc="none", mec=col, mew=0.7, zorder=z)
        lo = min(lo, s[:, k].min(), c[:, k].min())
        hi = max(hi, s[:, k].max(), c[:, k].max())

    for tag, xa, xl, yl in ANN[k]:
        s = fem[tag]
        ya = s[s[:, 0] <= xa][-1, k]
        txt = next(lab for t, _, lab in CASES if t == tag)
        a.annotate(txt, xy=(xa, ya), xytext=(xl, yl), textcoords="data",
                   ha="center", va="center", fontsize=8,
                   arrowprops=dict(arrowstyle="->", lw=0.6, shrinkA=0, shrinkB=2))

    a.set_ylabel(YLAB[i])
    a.set_xlabel(r"$s/s_{\mathrm{load}}$ [-]")
    a.set_xlim(0, 5)
    tidy_yaxis(a, lo, hi)
    a.xaxis.set_major_locator(MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
    a.xaxis.set_minor_locator(AutoMinorLocator(2))
    a.tick_params(which="minor", length=2)
    a.tick_params(which="major", length=4)
    a.text(0.03, 0.88, f"({'abc'[i]})", transform=a.transAxes,
           bbox=dict(fc="white", ec="none", pad=1.5))

h = [plt.Line2D([], [], color=RED, lw=1.0, marker="o", ms=3.2,
                mfc="none", mec=RED, label="MBE-KG"),
     plt.Line2D([], [], color="black", lw=1.0, marker="o", ms=3.2,
                mfc="none", mec="black", label="Full KG"),
     plt.Line2D([], [], color="0.35", lw=1.0, label="finite element"),
     plt.Line2D([], [], color="0.35", ls="none", marker="o", ms=3.2,
                mfc="none", mec="0.35", label="semi-analytical")]
fig.legend(handles=h, loc="upper center", ncol=4, fontsize=8, frameon=False,
           bbox_to_anchor=(0.5, 1.0), columnspacing=1.2, handletextpad=0.5)

fig.tight_layout(rect=(0, 0, 1, 0.92))
out = os.path.join(TEX, "tension-4fib-all-again.pdf")
fig.savefig(out)
print(f"wrote {out}")
for tag, _, _ in CASES:
    e = ana[tag][-1]
    print(f"  {tag:6s} end: sigma={e[1]:7.1f} kPa  lam/lam_o={e[2]:.3f}  lam_g={e[3]:.3f}")
