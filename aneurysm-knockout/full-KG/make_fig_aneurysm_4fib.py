#!/usr/bin/env python3
"""Figures 5 and 6 of the paper: full-KG (rate-dependent, s_load = 1, 10, 100 days)
against MBE-KG (rate-independent), for the axisymmetric (fusiform) and asymmetric
(saccular) lesions. Styled after the uniaxial figure (tension-4fib-all-again.pdf),
except that no semi-analytical solution exists here: every curve is a FE result.

Reads the FEBio node_data (x;ux;y;uy) and element_data (sx;sy;sz;J;sed) logs written
by the decks in this directory. The s_load = 1 day curves come from ds-aneu/, the FIXED-step
(Delta s = 0.001 day) runs: the vstep/ runs they replaced used an adaptive step (<dtmax lc="2">),
which is not the converged increment the text claims for every rate-dependent result. Note the file naming: the *-axisym-* series are written
by rate-fusiform-*.feb and the *-axis-* series by rate-saccular-*.feb.

The MBE solution is rate-independent, so its abscissa carries no time content; it is
drawn against the normalized fraction of the imposed degradation and then held.

    python3 make_fig_aneurysm_4fib.py           # writes both PDFs next to main.tex
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, FormatStrFormatter

TEX_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "..", "paper")
ONSET = 1.0  # growth onset (end of the pressurization stage)

plt.rcParams.update({
    "font.size": 9, "axes.linewidth": 0.8,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
})


def read_node(nf):
    lam, t = {}, None
    for ln in open(nf):
        ln = ln.strip()
        if ln.startswith("*Time"):
            t = float(ln.split("=")[1])
        elif ln and ln[0].isdigit() and len(ln.split(",")) >= 5:
            x, ux, y, uy = map(float, ln.split(",")[1:5])
            r = (x ** 2 + y ** 2) ** .5
            r0 = ((x - ux) ** 2 + (y - uy) ** 2) ** .5
            lam[round(t, 4)] = r / r0
    return lam


def read_elem(ef):
    sv, J, t = {}, {}, None
    for ln in open(ef):
        ln = ln.strip()
        if ln.startswith("*Time"):
            t = float(ln.split("=")[1])
        elif ln and ln[0].isdigit() and len(ln.split(",")) >= 5:
            p = ln.split(",")
            sx, sy, sz, j = (float(p[1]), float(p[2]), float(p[3]), float(p[4]))
            sv[round(t, 4)] = (sx + sy + sz) / 3.0
            J[round(t, 4)] = j
    return sv, J


def series(node_f, elem_f, s_load, mbe=False):
    lam = read_node(node_f)
    sv, J = read_elem(elem_f)
    ts = [t for t in sorted(set(lam) & set(sv)) if t >= ONSET]
    lro = lam[ts[0]]
    if mbe:
        # Rate-independent: the abscissa carries no time content. The MBE solve ramps
        # the degradation in over part of its pseudo-time and then holds. Map that ramp
        # onto [0, 1] so that its completion aligns with s/s_load = 1, as for the
        # transients, and hold the equilibrated value thereafter.
        Jend = J[ts[-1]]
        tol = 1e-4 * abs(Jend - 1.0)
        t_ramp = next(t for t in ts if abs(J[t] - Jend) <= tol)
        span = max(t_ramp - ONSET, 1e-9)
        ts = [t for t in ts if t <= t_ramp]
        x = [(t - ONSET) / span for t in ts]
    else:
        x = [(t - ONSET) / s_load for t in ts]
    d = dict(x=x, sv=[sv[t] for t in ts],
             lr=[lam[t] / lro for t in ts], J=[J[t] for t in ts])
    k = [i for i, xx in enumerate(d["x"]) if xx <= 5.0]      # el eje del paper llega a 5
    return {kk: [d[kk][i] for i in k] for kk in ("x", "sv", "lr", "J")}


def resample(d, key, n=26):
    """Circle markers at ~n evenly spaced abscissae, as in the uniaxial figure."""
    xs, ys = d["x"], d[key]
    if len(xs) <= n:
        return xs, ys
    step = (xs[-1] - xs[0]) / (n - 1)
    out_x, out_y, nxt, j = [], [], xs[0], 0
    for x, y in zip(xs, ys):
        if x >= nxt - 1e-12:
            out_x.append(x); out_y.append(y); nxt += step
    return out_x, out_y


CASES = {
    "aneurysm-4fib-axis.pdf": dict(
        mbe=("apex/node-axisym-mbe.txt", "apex/elem-axisym-mbe.txt"),
        full=[(1, "apex/ds-aneu/node-ds001.txt", "apex/ds-aneu/elem-ds001.txt"),
              (10, "apex/node-axisym-sp10.txt", "apex/elem-axisym-sp10.txt"),
              (100, "apex/node-axisym-sp100.txt", "apex/elem-axisym-sp100.txt")],
        # per panel: (s_load, anchor x on its curve, label x, label y) in data coordinates
        ann={"sv": [(1, 3.5, 3.8, 217.5), (10, 1.6, 2.5, 211.5), (100, 0.4, 2.0, 205.75)],
             "lr": [(1, 3.0, 3.5, 1.037), (10, 1.3, 2.8, 1.008), (100, 1.0, 1.6, 0.996)],
             "J":  [(1, 4.0, 3.8, 1.007), (10, 1.5, 2.8, 1.0175), (100, 0.9, 2.2, 1.0012)]}),
    "aneurysm-4fib-asym.pdf": dict(
        mbe=("apex/node-axis-FINALmbe.txt", "apex/elem-axis-FINALmbe.txt"),
        full=[(1, "apex/ds-aneu/node-axis-FINALfk-SDS.txt", "apex/ds-aneu/elem-axis-FINALfk-SDS.txt"),
              (10, "apex/node-axis-sp10.txt", "apex/elem-axis-sp10.txt"),
              (100, "apex/node-axis-sp100.txt", "apex/elem-axis-sp100.txt")],
        ann={"sv": [(1, 3.5, 3.8, 221.0), (10, 1.6, 2.5, 212.5), (100, 0.4, 2.0, 205.75)],
             "lr": [(1, 3.0, 3.5, 1.068), (10, 1.3, 2.8, 1.015), (100, 1.0, 1.6, 0.990)],
             "J":  [(1, 4.0, 3.8, 1.012), (10, 1.5, 2.8, 1.024), (100, 0.9, 2.2, 1.0018)]}),
}

KEYS = ["sv", "lr", "J"]
YLAB = [r"$\sigma_{\mathrm{v}}$ [kPa]", r"$\lambda_{r}/\lambda_{ro}$ [-]", r"$J$ [-]"]

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

for out, cfg in CASES.items():
    m = series(*cfg["mbe"], s_load=None, mbe=True)
    fulls = {sl: series(nf, ef, sl) for sl, nf, ef in cfg["full"]}

    fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.9))
    for i, k in enumerate(KEYS):
        a = ax[i]
        # MBE-KG: ramp over the imposed degradation, then held at the equilibrated value
        mx = list(m["x"]) + [5.0]
        my = list(m[k]) + [m[k][-1]]
        # every curve here is a finite-element result, so all are drawn as solid lines:
        # red = MBE-KG (equilibrated), black = full KG. Circles are reserved, in Fig. 2,
        # for the semi-analytical solution, which has no counterpart in this example.
        a.plot(mx, my, color="red", lw=1.0, zorder=3)
        lo, hi = min(my), max(my)
        for sl, d in fulls.items():
            a.plot(d["x"], d[k], color="black", lw=1.0, zorder=2)
            lo, hi = min(lo, min(d[k])), max(hi, max(d[k]))
        for sl, xa, xl, yl in cfg["ann"][k]:
            d = fulls[sl]
            ya = [v for x, v in zip(d["x"], d[k]) if x <= xa][-1]
            txt = r"$s_{\mathrm{load}} = 1$ day" if sl == 1 else f"{sl} days"
            a.annotate(txt, xy=(xa, ya), xytext=(xl, yl), textcoords="data",
                       ha="center", va="center", fontsize=8,
                       arrowprops=dict(arrowstyle="->", lw=0.6,
                                       shrinkA=0, shrinkB=2))
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

    h = [plt.Line2D([], [], color="red", lw=1.0, label="MBE-KG"),
         plt.Line2D([], [], color="black", lw=1.0, label="Full KG")]
    fig.legend(handles=h, loc="upper center", ncol=2, fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path = os.path.normpath(os.path.join(TEX_DIR, out))
    fig.savefig(path)
    print("wrote", path, "| MBE J =", round(m["J"][-1], 5),
          "| full-KG J(end) =", {sl: round(d["J"][-1], 5) for sl, d in fulls.items()})
