#!/usr/bin/env python3
"""Figure 3: stress homeostasis through growth, on the uniaxial sample.

Passive (growth suppressed, black) against equilibrated MBE-KG (red), as the axial force
is raised from the homeostatic value F_o. The house style of Figs. 2, 6 and 7 is used:
red = MBE-KG, black = the other model, solid line = finite element, circles = the
analytical result.

Three panels side by side against the applied force F/F_o:
  (A) axial Cauchy stress: the passive branch climbs, the equilibrated one does not move.
  (B) axial stretch: the passive branch saturates (collagen locks), the equilibrated one
      keeps stretching because it is growing.
  (C) relative volume J: where the load went instead, namely into mass. The circles are
      the CLOSED FORM J = (F/F_o)^{3/2}, which the finite-element solution reproduces.

The closed form: if growth is isotropic and the Cauchy stress is held at sigma_o, the
current cross-section scales as J^{2/3}, so F = sigma_o A_o J^{2/3}, hence J=(F/F_o)^{3/2}.

The axis is cut at 3 F_o. Equilibrated solutions were verified to exist up to 8 F_o, but
that is an existence bound, not a physiological one, and is reported in the text only.

Reads the FEBio logs of uniax-MBE8.feb (growth on) and uniax-PASS8.feb (gr_start beyond
the end of the run, so growth never activates).
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

RED = "red"          # MBE-KG (equilibrated), as in Figs. 2, 6 and 7
BLACK = "black"      # passive (growth suppressed)
FMAX = 3.0           # the plotted branch: F_o to 3 F_o


def load(tag):
    """Return F/F_o, lambda, sigma_zz [kPa] and J along the loading branch."""
    def block(path, ncol):
        t, d = None, {}
        for ln in open(os.path.join(HERE, path)):
            ln = ln.strip()
            if ln.startswith("*Time"):
                t = float(ln.split("=")[1])
            elif ln and ln[0].isdigit():
                p = ln.split(",")
                if len(p) >= ncol + 1:
                    d[t] = [float(x) for x in p[1:ncol + 1]]
        return d
    node = block(f"node7_{tag}.txt", 3)          # ux, uy, uz  (node 7 = the loaded corner)
    elem = block(f"elem1_{tag}.txt", 5)          # sx, sy, sz, J, Fgxx
    out = {k: [] for k in ("F", "lam", "sig", "J")}
    for t in sorted(set(node) & set(elem)):
        if t < 1.0:                              # before F_o: the pre-load ramp
            continue
        uz = node[t][2]                          # pulled along the axial direction
        _, _, sz, J, _ = elem[t]
        out["F"].append(1.0 + 7.0 * (t - 1.0) / 49.0)   # the deck's load curve
        out["lam"].append(1.0 + uz)              # one face fixed -> lambda = 1 + u
        out["sig"].append(sz)
        out["J"].append(J)
    return {k: np.array(v) for k, v in out.items()}


def tidy_yaxis(a, lo, hi, nbins=5):
    """First and last major tick exactly on the frame, so the extreme labels sit at the
    corners of the box (as in Figs. 2, 6 and 7)."""
    t = MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
    step = t[1] - t[0]
    d = next(k for k in range(6) if abs(round(step, k) - step) < 1e-9 * max(1.0, abs(step)))
    a.set_yticks(t)
    a.set_ylim(t[0], t[-1])
    a.yaxis.set_major_formatter(FormatStrFormatter(f"%.{d}f"))
    a.yaxis.set_minor_locator(AutoMinorLocator(2))


m, p = load("MBE8"), load("PASS8")
SIG_O = m["sig"][0]
cut = lambda d: {k: v[d["F"] <= FMAX + 1e-9] for k, v in d.items()}
m, p = cut(m), cut(p)

fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.7))

# (A) the stress: held at sigma_o by growth, climbing without it
ax[0].plot(p["F"], p["sig"], color=BLACK, lw=1.0, zorder=2)
ax[0].plot(m["F"], m["sig"], color=RED, lw=1.0, zorder=3)
ax[0].set_ylabel(r"$\sigma$ [kPa]")
# curve labels: no prefix; each in empty space, in the colour of its curve, with a thin
# arrow to the curve it names
lab = dict(ha="center", va="bottom", fontsize=8)
arr = lambda c: dict(arrowstyle="-|>", color=c, lw=0.6, mutation_scale=7, shrinkA=1, shrinkB=0)
ax[0].annotate(r"$\sigma_{o}=%.1f$ kPa" % SIG_O, xy=(2.7, SIG_O), xytext=(2.52, 400),
               color=RED, arrowprops=arr(RED), **lab)
tidy_yaxis(ax[0], min(p["sig"].min(), m["sig"].min()), p["sig"].max())

# (B) the stretch: the passive branch locks, the equilibrated one keeps going
ax[1].plot(p["F"], p["lam"], color=BLACK, lw=1.0, zorder=2)
ax[1].plot(m["F"], m["lam"], color=RED, lw=1.0, zorder=3)
ax[1].set_ylabel(r"$\lambda$ [-]")
tidy_yaxis(ax[1], min(p["lam"].min(), m["lam"].min()),
           max(p["lam"].max(), m["lam"].max()))

# (C) where the load went: into mass. Circles = the closed form J = (F/F_o)^{3/2}
Fc = np.linspace(1.0, FMAX, 21)
ax[2].plot(p["F"], p["J"], color=BLACK, lw=1.0, zorder=2)
ax[2].plot(m["F"], m["J"], color=RED, lw=1.0, zorder=3)
ax[2].plot(Fc, Fc ** 1.5, ls="none", marker="o", ms=3.6, mfc="none",
           mec=RED, mew=0.9, zorder=4)
ax[2].set_ylabel(r"$J$ [-]")
# in the empty upper-left of the panel, clear of every curve
ax[2].annotate(r"$(F/F_{o})^{3/2}$", xy=(1.9, 1.9 ** 1.5), xytext=(1.45, 4.3),
               color=RED, arrowprops=arr(RED), **lab)
ax[2].annotate(r"$J=1$", xy=(2.5, 1.0), xytext=(2.45, 1.75),
               color=BLACK, arrowprops=arr(BLACK), **lab)
tidy_yaxis(ax[2], 0.0, m["J"].max())   # start at 0 so the passive J=1 line is not on the frame

for i, a in enumerate(ax):
    a.set_xlabel(r"applied force  $F/F_{o}$ [-]")
    a.set_xlim(1, FMAX)
    a.xaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
    a.xaxis.set_minor_locator(AutoMinorLocator(2))
    a.tick_params(which="minor", length=2)
    a.tick_params(which="major", length=4)
    a.text(0.03, 0.88, f"({'abc'[i]})", transform=a.transAxes,
           bbox=dict(fc="white", ec="none", pad=1.5))

h = [plt.Line2D([], [], color=RED, lw=1.0, label="MBE-KG (equilibrated)"),
     plt.Line2D([], [], color=BLACK, lw=1.0, label="passive (no growth)"),
     plt.Line2D([], [], color=RED, ls="none", marker="o", ms=3.6, mfc="none",
                mec=RED, label="closed form")]
fig.legend(handles=h, loc="upper center", ncol=3, fontsize=8, frameon=False,
           bbox_to_anchor=(0.5, 1.0))

fig.tight_layout(rect=(0, 0, 1, 0.92))
out = os.path.join(TEX, "uniax-mbe-homeostasis.pdf")
fig.savefig(out)

err = 100 * np.abs(m["J"] - m["F"] ** 1.5) / m["F"] ** 1.5
print(f"wrote {out}")
print(f"  sigma_o = {SIG_O:.1f} kPa")
print(f"  at {FMAX:.0f} F_o:  MBE  sigma = {m['sig'][-1]:.1f} kPa, lambda = {m['lam'][-1]:.3f}, "
      f"J = {m['J'][-1]:.3f}")
print(f"              passive  sigma = {p['sig'][-1]:.1f} kPa, lambda = {p['lam'][-1]:.3f}, "
      f"J = {p['J'][-1]:.4f}")
print(f"  closed form J=(F/Fo)^(3/2): max error over the branch = {err.max():.4f}%")
