#!/usr/bin/env python3
"""Figure 3: uniaxial coupon over the load range F_o <= F <= 8 F_o (red = MBE-KG,
black = hyperelastic, line = finite element, open circles = closed form J = (F/F_o)^{3/2}).
Usage:

    python3 make_fig_homeostasis_s1.py [A|B|C|D] [--overwrite-v2]

  A  linear axes everywhere (default; the figure of the paper).
  B  log F/F_o in all panels; log sigma in (a); log-log in (c), where the closed form is
     a straight line of slope 3/2; linear lambda in (b).
  C  as B, with the part of the equilibrated curve beyond J = 4 (F/F_o > 4^{2/3}) in grey.
  D  as B, with the three stretch ratios lambda_i/lambda_oi of both responses in (b).

Writes uniax-mbe-homeostasis-s1.pdf into the manuscript folder and a PNG next to this script;
refuses to overwrite either unless --overwrite-v2 is given.
Reads node7_MBE8.txt and elem1_MBE8.txt from this folder and node7_PASS8.txt and elem1_PASS8.txt
from ../../uniax-fractions/ (ux, uy, uz of node 7; sx, sy, sz, J, Fgxx of element 1),
(sx, sy, sz, J, Fgxx), as make_numbers.py does: x = radial, y = circumferential,
z = axial, lambda_i = 1 + u_i (one face of the unit cube fixed per direction).
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import (AutoMinorLocator, MaxNLocator, FormatStrFormatter,
                               FixedLocator, NullFormatter, LogLocator, NullLocator)

HERE = os.path.dirname(os.path.abspath(__file__))

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8,
                     "xtick.direction": "in", "ytick.direction": "in",
                     "xtick.top": True, "ytick.right": True})

RED = "red"
BLACK = "black"
GREY = "0.6"
FMAX = 8.0
J_ADM = 4.0
F_ADM = J_ADM ** (2.0 / 3.0)        # 2.52: closed form J = (F/F_o)^{3/2} = 4

VAR = sys.argv[1].upper() if len(sys.argv) > 1 else "A"
assert VAR in "ABCD" and len(VAR) == 1, "variant must be A, B, C or D"
LOG = VAR in "BCD"


def load(tag):
    """F/F_o, sigma_zz [kPa], J and the stretches (axial, circ, radial) on the branch."""
    def block(path, ncol):
        t, d = None, {}
        # hyperelastic run (no growth, unaffected by lo_opt) is read from uniax-fractions
        src = HERE if "MBE8" in path else os.path.join(HERE, "..", "..", "uniax-fractions")
        for ln in open(os.path.join(src, path)):
            ln = ln.strip()
            if ln.startswith("*Time"):
                t = float(ln.split("=")[1])
            elif ln and ln[0].isdigit():
                p = ln.split(",")
                if len(p) >= ncol + 1:
                    d[t] = [float(x) for x in p[1:ncol + 1]]
        return d
    node = block(f"node7_{tag}.txt", 3)
    elem = block(f"elem1_{tag}.txt", 5)
    out = {k: [] for k in ("F", "sig", "J", "lz", "ly", "lx")}
    for t in sorted(set(node) & set(elem)):
        if t < 1.0:
            continue
        ux, uy, uz = node[t]
        _, _, sz, J, _ = elem[t]
        out["F"].append(1.0 + 7.0 * (t - 1.0) / 49.0)
        out["sig"].append(sz)
        out["J"].append(J)
        out["lz"].append(1.0 + uz)      # axial
        out["ly"].append(1.0 + uy)      # circumferential
        out["lx"].append(1.0 + ux)      # radial
    d = {k: np.array(v) for k, v in out.items()}
    d["lam"] = d["lz"]
    return d


def tidy_yaxis(a, lo, hi, nbins=5):
    t = MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
    step = t[1] - t[0]
    d = next(k for k in range(6) if abs(round(step, k) - step) < 1e-9 * max(1.0, abs(step)))
    a.set_yticks(t)
    a.set_ylim(t[0], t[-1])
    a.yaxis.set_major_formatter(FormatStrFormatter(f"%.{d}f"))
    a.yaxis.set_minor_locator(AutoMinorLocator(2))


def log_yaxis(a, lo, hi, ticks):
    a.set_yscale("log")
    a.set_ylim(lo, hi)
    a.yaxis.set_major_locator(FixedLocator(ticks))
    a.yaxis.set_major_formatter(FormatStrFormatter("%g"))
    a.yaxis.set_minor_locator(NullLocator())     # v2: no unlabelled minor clutter


def plot_mbe(a, x, y, **kw):
    """Equilibrated curve; in variant C the part beyond J = 4 is grey."""
    if VAR == "C":
        k = np.searchsorted(x, F_ADM)
        yb = np.interp(F_ADM, x, y)
        a.plot(np.r_[x[:k], F_ADM], np.r_[y[:k], yb], color=RED, lw=1.0, zorder=3, **kw)
        a.plot(np.r_[F_ADM, x[k:]], np.r_[yb, y[k:]], color=GREY, lw=1.0, zorder=3, **kw)
    else:
        a.plot(x, y, color=RED, lw=1.0, zorder=3, **kw)


m, p = load("MBE8"), load("PASS8")
SIG_O = m["sig"][0]

fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.7))
lab = dict(ha="center", va="bottom", fontsize=8)
arr = lambda c: dict(arrowstyle="-|>", color=c, lw=0.6, mutation_scale=7, shrinkA=1, shrinkB=0)

# ---------------------------------------------------------------- (a) stress
ax[0].plot(p["F"], p["sig"], color=BLACK, lw=1.0, zorder=2)
plot_mbe(ax[0], m["F"], m["sig"])
ax[0].set_ylabel(r"$\sigma$ [kPa]")
if LOG:
    log_yaxis(ax[0], 200, 3200, [200, 400, 800, 1600, 3200])   # v2: ends on labelled ticks
    # in C the arrow lands on the red (J <= 4) part of the line
    xa, xt = ((2.1, 4.2) if VAR == "C" else (5.0, 4.2))
    ax[0].annotate(r"$\sigma_{o}=%.1f$ kPa" % SIG_O, xy=(xa, SIG_O), xytext=(xt, 420),
                   color=RED, arrowprops=arr(RED), **lab)
else:
    tidy_yaxis(ax[0], 0.0, p["sig"].max())
    ax[0].annotate(r"$\sigma_{o}=%.1f$ kPa" % SIG_O, xy=(6.5, SIG_O), xytext=(5.6, 700),
                   color=RED, arrowprops=arr(RED), **lab)

# ---------------------------------------------------------------- (b) stretch
if VAR == "D":
    for d in (m, p):
        for k in ("lz", "ly", "lx"):
            d["r" + k] = d[k] / d[k][0]
    # hyperelastic: axial solid, circumferential dashed, radial dash-dotted
    ax[1].plot(p["F"], p["rlz"], color=BLACK, lw=1.0, zorder=2)
    ax[1].plot(p["F"], p["rly"], color=BLACK, lw=1.0, ls="--", zorder=2)
    ax[1].plot(p["F"], p["rlx"], color=BLACK, lw=1.0, ls="-.", zorder=2)
    # equilibrated: axial line, the two transverse ones as sparse markers on it
    ax[1].plot(m["F"], m["rlz"], color=RED, lw=1.0, zorder=3)
    Fk = np.geomspace(1.0, FMAX, 12)
    ax[1].plot(Fk, np.interp(Fk, m["F"], m["rly"]), ls="none", marker="s", ms=3.4,
               mfc="none", mec=RED, mew=0.8, zorder=4)
    Fk2 = np.sqrt(Fk[:-1] * Fk[1:])                  # staggered, so the shapes do not stack
    ax[1].plot(Fk2, np.interp(Fk2, m["F"], m["rlx"]), ls="none", marker="^", ms=3.6,
               mfc="none", mec=RED, mew=0.8, zorder=4)
    ax[1].set_ylabel(r"$\lambda_i/\lambda_{oi}$ [-]")
    lo = min(p[k].min() for k in ("rly", "rlx"))
    tidy_yaxis(ax[1], np.floor(lo * 10) / 10, m["rlz"].max())
    hb = [plt.Line2D([], [], color="0.25", lw=1.0, label="axial"),
          plt.Line2D([], [], color="0.25", lw=1.0, ls="--", marker="s", ms=3.4,
                     mfc="none", mew=0.8, label="circ."),
          plt.Line2D([], [], color="0.25", lw=1.0, ls="-.", marker="^", ms=3.6,
                     mfc="none", mew=0.8, label="radial")]
    ax[1].legend(handles=hb, loc="upper left", bbox_to_anchor=(0.0, 0.86), fontsize=7,
                 frameon=False, handlelength=2.6, borderaxespad=0.4, labelspacing=0.25)
else:
    ax[1].plot(p["F"], p["lam"], color=BLACK, lw=1.0, zorder=2)
    plot_mbe(ax[1], m["F"], m["lam"])
    ax[1].set_ylabel(r"$\lambda$ [-]")
    tidy_yaxis(ax[1], min(p["lam"].min(), m["lam"].min()),
               max(p["lam"].max(), m["lam"].max()))

# ---------------------------------------------------------------- (c) volume
ax[2].plot(p["F"], p["J"], color=BLACK, lw=1.0, zorder=2)
plot_mbe(ax[2], m["F"], m["J"])
Fc = np.geomspace(1.0, FMAX, 12) if LOG else np.linspace(1.0, FMAX, 12)
ax[2].plot(Fc, Fc ** 1.5, ls="none", marker="o", ms=3.6, mfc="none",
           mec=RED, mew=0.9, zorder=4)
ax[2].set_ylabel(r"$J$ [-]")
if LOG:
    log_yaxis(ax[2], 0.5, 32, [0.5, 1, 2, 4, 8, 16, 32])   # v2: ends on labelled ticks
    ax[2].annotate(r"$(F/F_{o})^{3/2}$", xy=(2.0, 2.0 ** 1.5), xytext=(1.5, 8.0),
                   color=RED, arrowprops=arr(RED), **lab)
    ax[2].annotate(r"$J=1$", xy=(5.0, 1.0), xytext=(5.0, 1.9),
                   color=BLACK, arrowprops=arr(BLACK), **lab)
else:
    tidy_yaxis(ax[2], 0.0, m["J"].max())
    ax[2].annotate(r"$(F/F_{o})^{3/2}$", xy=(4.0, 8.0), xytext=(2.6, 14.0),
                   color=RED, arrowprops=arr(RED), **lab)
    ax[2].annotate(r"$J=1$", xy=(6.5, 1.0), xytext=(6.5, 4.5),
                   color=BLACK, arrowprops=arr(BLACK), **lab)

if VAR == "C":
    for a in ax:
        a.axvline(F_ADM, color="0.35", lw=0.6, ls=":", zorder=1)
    ax[2].text(F_ADM * 1.04, 22.0, r"$J=4$", fontsize=7, color="0.35", ha="left", va="center")

for i, a in enumerate(ax):
    a.set_xlabel(r"applied force  $F/F_{o}$ [-]")
    if LOG:
        a.set_xscale("log")
        a.set_xlim(1, FMAX)
        a.xaxis.set_major_locator(FixedLocator([1, 2, 3, 4, 6, 8]))
        a.xaxis.set_major_formatter(FormatStrFormatter("%g"))
        a.xaxis.set_minor_locator(FixedLocator([1.5, 2.5, 5, 7]))
        a.xaxis.set_minor_formatter(NullFormatter())
    else:
        a.set_xlim(1, FMAX)
        # rABC-6 (3 Oct): one major tick per integer force ratio (1, 2, 4, 6, 8 crowded 1 and 2)
        a.xaxis.set_major_locator(FixedLocator(list(range(1, int(FMAX) + 1))))
        a.xaxis.set_minor_locator(NullLocator())
    a.tick_params(which="minor", length=2)
    a.tick_params(which="major", length=4)
    a.text(0.03, 0.88, f"({'abc'[i]})", transform=a.transAxes,
           bbox=dict(fc="white", ec="none", pad=1.5))

h = [plt.Line2D([], [], color=RED, lw=1.0, label="MBE-KG (equilibrated)"),
     plt.Line2D([], [], color=BLACK, lw=1.0, label="hyperelastic (growth disabled)"),
     plt.Line2D([], [], color=RED, ls="none", marker="o", ms=3.6, mfc="none",
                mec=RED, label="closed form")]
if VAR == "C":
    h.insert(1, plt.Line2D([], [], color=GREY, lw=1.0, label=r"MBE-KG, $J>4$"))
fig.legend(handles=h, loc="upper center", ncol=len(h), fontsize=8, frameon=False,
           bbox_to_anchor=(0.5, 1.0))

fig.tight_layout(rect=(0, 0, 1, 0.92))
# manuscript PDF and an inspection PNG next to the script
pdf = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", "paper", "uniax-mbe-homeostasis-s1.pdf"))
png = os.path.join(HERE, "uniax-mbe-homeostasis-s1.png")
for q in (pdf, png):
    assert not os.path.exists(q) or "--overwrite-v2" in sys.argv, f"{q} exists"
fig.savefig(pdf)
fig.savefig(png, dpi=150)
print(f"wrote {pdf}\n      {png}")

# ---------------------------------------------------------------- numbers to cross-check
for F0 in (3.0, 8.0):
    for nm, d in (("MBE ", m), ("PASS", p)):
        print(f"  F={F0:.0f}Fo {nm}: sigma={np.interp(F0, d['F'], d['sig']):9.2f} kPa  "
              f"lam_z={np.interp(F0, d['F'], d['lz']):.4f}  J={np.interp(F0, d['F'], d['J']):.4f}")
print(f"  sigma_o = {SIG_O:.2f} kPa; F at last record: MBE {m['F'][-1]:.4f}, PASS {p['F'][-1]:.4f}")
print(f"  max |J_hyperelastic - 1| = {np.max(abs(p['J'] - 1)):.3e}")
J3 = m["J"] ** (1.0 / 3.0)
for k in ("lz", "ly", "lx"):
    r = m[k] / m[k][0]
    print(f"  MBE {k}/{k}_o vs J^(1/3): max rel diff = {np.max(abs(r - J3) / J3):.3e}")
print(f"  closed form: max rel err = {np.max(abs(m['J'] - m['F']**1.5) / m['F']**1.5):.3e}")
print(f"  MBE max |sigma-sigma_o|/sigma_o = {np.max(abs(m['sig'] - SIG_O)) / SIG_O:.3e}")
