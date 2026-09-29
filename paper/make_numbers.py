"""Single generator of every result quoted in main.tex. Writes numbers.tex and the bodies
of the two convergence tables (macros convUniRows and convAneuRows).

Sources, all under ../aneurysm-knockout/ (AK):
  AK/uniax-fractions/        uniaxial coupon: FEBio logs (elem1_*, node7_*), the plot-ready
                             series {MBE,sp1,sp10,sp100}_c.txt (FE) and _a.txt
                             (semi-analytical, ../semi-analytical/uniax_4fib.py), ds-study/
  AK/kJ-study/               volumetric-penalty sweep
  AK/mesh-study/             radial mesh refinement of the two lesions
  AK/param-x/                the four lesion fields (matrix-*), and the parametric sweep
                             MASTER_PARAMETRIC_STUDY.csv
  AK/full-KG/                rate comparison (MBE and full-KG), ds-aneu/, ds-study/
  AK/timing/                 wall time, steps and Newton iterations, and the logs whose
                             residuals fill the convergence tables (timing.py)

Every file opened for reading is written to numbers-inputs.txt, so that build.sh can
refuse to compile when an input is newer than numbers.tex.

Run: python3 make_numbers.py
"""
import builtins
import csv
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
AK = os.path.join(ROOT, "aneurysm-knockout")
MACROS = {}

_open, INPUTS = builtins.open, set()


def _recording_open(file, mode="r", *a, **k):
    if isinstance(file, str) and "r" in mode and "+" not in mode:
        p = os.path.abspath(file)
        if p.startswith(ROOT + os.sep) and not p.endswith(".py"):
            INPUTS.add(os.path.relpath(p, HERE))
    return _open(file, mode, *a, **k)


builtins.open = _recording_open


def put(name, value):
    if name in MACROS:
        raise SystemExit(f"macro {name} defined twice")
    if not re.fullmatch(r"[A-Za-z]+", name):
        raise SystemExit(f"macro name {name} must contain letters only")
    MACROS[name] = value


def pct(x, d=1):
    return f"{100.0 * x:.{d}f}"


def bound(x):
    """Smallest of 0.001, 0.01, 0.1, 1 (in percent) above the fraction x, for 'less than'."""
    return next(f"{b:g}" for b in (0.001, 0.01, 0.1, 1.0, 10.0) if 100.0 * x < b)


def sci(x, d=2):
    """Table format of the residuals: 3.81E+02."""
    return f"{x:.{d}E}"


def blocks(path):
    """FEBio logfile: '*Time = t' headers, then 'id,v1,v2,...' rows.
    Returns a sorted list of (t, {id: [v1, v2, ...]})."""
    out, t = [], None
    for ln in open(path):
        s = ln.strip()
        if s.startswith("*Time"):
            t = float(s.split("=")[1])
            out.append((t, {}))
        elif s and s[0].isdigit() and t is not None:
            p = s.split(",")
            out[-1][1][int(p[0])] = [float(v) for v in p[1:]]
    return out


def at(bl, t):
    """Block whose time is closest to t."""
    return min(bl, key=lambda b: abs(b[0] - t))[1]


def one(bl):
    """Series of a single-entity log: times and an array of values."""
    ts = np.array([b[0] for b in bl])
    vs = np.array([next(iter(b[1].values())) for b in bl])
    return ts, vs


def ak(*p):
    return os.path.join(AK, *p)


# =================================================================== uniaxial coupon
U = ak("uniax-fractions")
T_O = 1.0                                   # end of the pre-load: F = F_o
te, e = one(blocks(os.path.join(U, "elem1_MBE8.txt")))      # sx, sy, sz, J, Fg
tn, n = one(blocks(os.path.join(U, "node7_MBE8.txt")))      # ux, uy, uz
io = int(np.argmin(abs(te - T_O)))
jo = int(np.argmin(abs(tn - T_O)))
sig_o = e[io, 2]
lam_o = 1.0 + n[jo, 2]                      # axial (pull along z = N_z)
lam_ci = 1.0 + n[jo, 1]                     # circumferential (y = N_theta)
lam_ra = 1.0 + n[jo, 0]                     # radial (x = N_r)
put("uniSigo", f"{sig_o:.1f}")
xs_ = [float(m.group(1)) for m in re.finditer(r'<node id="\d+">([-\d.]+),', open(os.path.join(U, "uniax-MBE8.feb")).read())]
xc, g = 0.5 * (min(xs_) + max(xs_)), 0.5 / np.sqrt(3.0)
put("uniOffset", f"{xc:.0f}")
put("uniFibreSpread", f"{2 * np.degrees(np.arctan(g / (xc - g))):.2f}")
put("uniLamo", f"{lam_o:.2f}")
put("uniWoz", f"{lam_o:.2f}")
put("uniWoc", f"{lam_ci:.2f}")
put("uniWor", f"{lam_ra:.2f}")


def branch(tag):
    """F/F_o, lambda, sigma_zz, J and the three stretch ratios on the loading branch.
    The load curve of uniax-{MBE8,PASS8}.feb raises F linearly from F_o at t = 1 to
    8 F_o at t = 50."""
    tb, eb = one(blocks(os.path.join(U, f"elem1_{tag}.txt")))
    tnb, nb = one(blocks(os.path.join(U, f"node7_{tag}.txt")))
    common = sorted(set(np.round(tb, 9)) & set(np.round(tnb, 9)))
    common = [t for t in common if t >= T_O - 1e-9]
    E = {round(t, 9): v for t, v in zip(tb, eb)}
    N = {round(t, 9): v for t, v in zip(tnb, nb)}
    F = np.array([1.0 + 7.0 * (t - 1.0) / 49.0 for t in common])
    sig = np.array([E[t][2] for t in common])
    J = np.array([E[t][3] for t in common])
    st = np.array([[1.0 + N[t][2], 1.0 + N[t][1], 1.0 + N[t][0]] for t in common])
    return F, sig, J, st


Fp, sp, Jp, stp = branch("PASS8")
Fm, sm, Jm, stm = branch("MBE8")
put("uniPassSigThree", f"{np.interp(3.0, Fp, sp):.1f}")
put("uniPassLamThree", f"{np.interp(3.0, Fp, stp[:, 0]):.2f}")
put("uniPassLamEight", f"{np.interp(8.0, Fp, stp[:, 0]):.2f}")
put("uniPassSigEight", f"{np.interp(8.0, Fp, sp) / 1000.0:.1f}")
put("uniMbeLamThree", f"{np.interp(3.0, Fm, stm[:, 0]):.2f}")
put("uniMbeJThree", f"{np.interp(3.0, Fm, Jm):.2f}")
put("uniMbeJEight", f"{np.interp(8.0, Fm, Jm):.1f}")
put("uniMbeLamEight", f"{np.interp(8.0, Fm, stm[:, 0]):.2f}")
put("uniMbeSigSpread", bound(np.max(abs(sm - sig_o)) / sig_o))
closed = np.max(abs(Jm - Fm ** 1.5) / Fm ** 1.5)
put("uniClosedErr", pct(closed, 3))
iso = np.max(abs(stm / stm[0] - Jm[:, None] ** (1.0 / 3.0)) / Jm[:, None] ** (1.0 / 3.0))
put("uniIsoErr", pct(iso, 3))


def series(tag, kind):
    d = np.loadtxt(os.path.join(U, f"{tag}_{kind}.txt"))   # x, sigma, lam/lam_o, lam_g
    return d[d[:, 0] >= 0.0]


fe = {t: series(t, "c") for t in ("MBE", "sp1", "sp10", "sp100")}
sa = {t: series(t, "a") for t in ("MBE", "sp1", "sp10", "sp100")}
SLOAD = {"MBE": 10.0, "sp1": 1.0, "sp10": 10.0, "sp100": 100.0}
dsig, dlam = 0.0, 0.0
for t in fe:
    a, f = sa[t], fe[t]
    for k, col in ((1, "sig"), (2, "lam"), (3, "lam")):
        fi = np.interp(a[:, 0], f[:, 0], f[:, k])
        m = abs(a[:, k]) > 1e-9
        r = np.max(abs(fi[m] - a[m, k]) / abs(a[m, k]))
        if col == "sig":
            dsig = max(dsig, r)
        else:
            dlam = max(dlam, r)
put("uniFeSemiSig", pct(dsig, 1))
put("uniFeSemiLam", pct(dlam, 1))


def end(tag, k):
    f = fe[tag]
    return np.interp(5.0, f[:, 0], f[:, k])


J_end = {t: end(t, 3) ** 3 for t in fe}
put("uniVolInc", f"{100.0 * (J_end['MBE'] - 1.0):.0f}")
put("uniShortOne", f"{100.0 * (1.0 - J_end['sp1'] / J_end['MBE']):.0f}")
put("uniShortTen", pct(1.0 - J_end["sp10"] / J_end["MBE"], 1))
put("uniShortHundred", bound(abs(1.0 - J_end["sp100"] / J_end["MBE"])))
put("uniStressOne", f"{100.0 * abs(end('sp1', 1) - end('MBE', 1)) / end('MBE', 1):.0f}")


def settle(tag, k, ref, tol=0.05):
    """Days after the onset of loading from which the quantity stays within tol of ref."""
    f = fe[tag]
    bad = np.nonzero(abs(f[:, k] - ref) / abs(ref) >= tol)[0]
    i = bad[-1] + 1 if len(bad) else 0
    return f[min(i, len(f) - 1), 0] * SLOAD[tag]


lg_mbe = end("MBE", 3)
put("uniSettleSigTen", f"{settle('sp10', 1, sig_o):.0f}")
put("uniSettleSigHundred", f"{settle('sp100', 1, sig_o):.0f}")
put("uniSettleGrowTen", f"{settle('sp10', 3, lg_mbe):.0f}")
put("uniSettleGrowHundred", f"{settle('sp100', 3, lg_mbe):.0f}")

# time-step refinement of the rate-dependent coupon (s_load = 10 days)
DS = ["0.005", "0.01", "0.02", "0.05", "0.1", "0.2", "0.5", "1.0"]
Jds = {}
for ds in DS:
    tb, eb = one(blocks(os.path.join(U, "ds-study", f"elem1_ds{ds}.txt")))
    Jds[ds] = eb[-1, 3]
put("uniDsChange", pct(abs(Jds["0.5"] - Jds["0.1"]) / Jds["0.1"], 1))

# ================================================================= volumetric penalty
KJ = ["1e3", "1e4", "1e5", "1e6", "1e7", "1e8"]


def kj_case(tag):
    bl = blocks(tag)
    s1 = at(bl, 1.0)
    v = next(iter(s1.values()))
    svo = (v[0] + v[1] + v[2]) / 3.0
    return v[3], svo, next(iter(bl[-1][1].values()))[3]


kj = {}
for k in KJ:
    kj[("uniax", k)] = kj_case(ak("kJ-study", f"uniax_kJ{k}_elem1_MBE.txt"))
    # element 10 (eA) is the apex, the point of maximum dilatation of the fusiform lesion;
    # element 1 (eB) sits at the contracted end of the segment
    kj[("aneu", k)] = kj_case(ak("kJ-study", f"aneu_kJ{k}_matrix-fus-el-eA.txt"))
    # Stage-I constraint: the larger deviation |J - 1| of the two probed elements
    _jb = kj_case(ak("kJ-study", f"aneu_kJ{k}_matrix-fus-el-eB.txt"))[0]
    if abs(_jb - 1.0) > abs(kj[("aneu", k)][0] - 1.0):
        kj[("aneu", k)] = (_jb,) + tuple(kj[("aneu", k)][1:])
plateau = ["1e5", "1e6", "1e7"]
put("kJConstraint", "10^{-4}" if max(abs(kj[(c, k)][0] - 1.0)
                                    for c in ("uniax", "aneu") for k in ("1e6", "1e7", "1e8")) < 1e-4
    else "FAIL")
_v = max(abs(kj[(c, "1e5")][0] - 1.0) for c in ("uniax", "aneu"))
_e = int(np.floor(np.log10(_v)))
put("kJFiveConstraint", rf"{_v / 10 ** _e:.1f}\times10^{{{_e}}}")
put("kJSvoVar", pct(max(abs(kj[(c, k)][1] / kj[(c, "1e6")][1] - 1.0)
                        for c in ("uniax", "aneu") for k in plateau), 2))
put("kJGrowthVar", pct(max(abs(kj[(c, k)][2] / kj[(c, "1e6")][2] - 1.0)
                           for c in ("uniax", "aneu") for k in plateau), 2))
put("kJLowConstraint", pct(max(abs(kj[(c, "1e3")][0] - 1.0) for c in ("uniax", "aneu")), 1))
put("kJLowGrowth", f"{100.0 * max(abs(kj[(c, '1e3')][2] / kj[(c, '1e6')][2] - 1.0) for c in ('uniax', 'aneu')):.0f}")

# ======================================================================= mesh study
for m, name in (("fus", "Fus"), ("sac", "Sac")):
    Jmax = []
    for nr in (1, 2, 3):
        bl = blocks(ak("mesh-study", f"{m}_nr{nr}_matrix-{m}-el-ALL.txt"))
        Jmax.append(max(v[0] for v in bl[-1][1].values()))
    for nr, word in zip((1, 2, 3), ("One", "Two", "Three")):
        put(f"mesh{name}{word}", f"{Jmax[nr - 1]:.3f}")
    put(f"mesh{name}Change", pct((Jmax[2] - Jmax[0]) / Jmax[0], 1))

# ===================================================================== lesion fields
P = ak("param-x")
CASES = {"FusEl": "matrix-fus-el", "FusCc": "matrix-fus-elcc",
         "SacEl": "matrix-sac-el", "SacCc": "matrix-sac-elcc"}


def ref_radii(feb):
    rad, inn = {}, False
    for line in open(feb):
        s = line.strip()
        if s.startswith("<Nodes"):
            inn = True
            continue
        if s.startswith("</Nodes>"):
            break
        if inn and s.startswith("<node"):
            nid = int(s.split('id="')[1].split('"')[0])
            x, y, _ = [float(v) for v in s.split(">", 1)[1].split("<")[0].split(",")]
            rad[nid] = np.hypot(x, y)
    return rad


def outer_diam(rad, pts):
    ids = np.array(sorted(rad))
    Rr = np.array([rad[i] for i in ids])
    outer = Rr > Rr.max() - 0.01
    d = [2.0 * np.hypot(*pts[i][:2]) for i in ids[outer]]
    return min(d), max(d)


def node_lamr(path):
    """Radial stretch of node A (the apex), normalized by its homeostatic value, exactly
    as read_all_fused.py computes it for the parametric sweep."""
    ts, v = one(blocks(path))                      # x, ux, y, uy
    r = np.hypot(v[:, 0], v[:, 2])
    r0 = np.hypot(v[:, 0] - v[:, 1], v[:, 2] - v[:, 3])
    lr = r / r0
    return lr[-1] / lr[int(np.argmin(abs(ts - 1.0)))]


def inner_lamr(rad, nb):
    """Largest radial stretch of the inner (luminal) wall at the end of the insult,
    relative to the homeostatic state (t = 1)."""
    ids = sorted(rad)
    Rr = np.array([rad[i] for i in ids])
    inner = [i for i, r in zip(ids, Rr) if r < Rr.min() + 1e-3]
    h, f = at(nb, 1.0), nb[-1][1]
    return max(np.hypot(*f[i][:2]) / np.hypot(*h[i][:2]) for i in inner)


def node_ratio(nb, nid):
    """Radial position of node nid at the end of the insult over its homeostatic value."""
    h, f = at(nb, 1.0), nb[-1][1]
    return np.hypot(*f[nid][:2]) / np.hypot(*h[nid][:2])


field = {}
for key, tag in CASES.items():
    bl = blocks(os.path.join(P, f"{tag}-ALL.txt"))                 # J, W_R
    fin, hom = bl[-1][1], at(bl, 1.0)
    J = {i: v[0] for i, v in fin.items()}
    W = {i: v[1] for i, v in fin.items()}
    imax, imin = max(J, key=J.get), min(J, key=J.get)
    nb = blocks(os.path.join(P, f"{tag}-NODES.txt"))
    rad = ref_radii(os.path.join(P, f"{tag}.feb"))
    field[key] = dict(
        Jmax=J[imax], Jmin=J[imin], WRmax=max(W.values()), WRends=W[imin],
        Wapex=W[imax] / J[imax], WRhom=np.mean([v[1] for v in hom.values()]),
        Whom=np.mean([v[1] / v[0] for v in hom.values()]),
        lamr=node_ratio(nb, 20),
        lamrMaxInner=inner_lamr(rad, nb),
        lamrA=node_lamr(os.path.join(P, f"{tag}-nA.txt")),
        dmax=outer_diam(rad, {i: v for i, v in nb[-1][1].items()})[1],
        dhom=outer_diam(rad, {i: v for i, v in at(nb, 1.0).items()})[0])
for key, f in field.items():
    put(f"field{key}Jmax", f"{f['Jmax']:.2f}")
    put(f"field{key}Jmin", f"{f['Jmin']:.2f}")
    put(f"field{key}Lamr", f"{f['lamr']:.2f}")
    put(f"field{key}Diam", f"{f['dmax']:.1f}")
    put(f"field{key}WRmax", f"{f['WRmax']:.0f}")
    put(f"field{key}WRends", f"{f['WRends']:.0f}")
    put(f"field{key}Wapex", f"{f['Wapex']:.0f}")
# in-vivo wall stresses of the calibrated aorta at its original homeostatic state (t = 1),
# element 10 (mid-length, insult axis): the axial one justifies the uniaxial load F_o (C66)
_h = at(blocks(os.path.join(P, "matrix-fus-el-eA.txt")), 1.0)[10]
put("aneuSttO", f"{_h[0]:.0f}")
put("aneuSzzo", f"{_h[2]:.0f}")
put("fieldWhom", f"{field['FusEl']['Whom']:.0f}")
put("fieldWRhom", f"{field['FusEl']['WRhom']:.0f}")
put("fieldDiamHom", f"{field['FusEl']['dhom']:.2f}")
put("fieldFusElVolInc", f"{100 * (field['FusEl']['Jmax'] - 1):.0f}")
put("fieldSacElVolInc", f"{100 * (field['SacEl']['Jmax'] - 1):.0f}")
put("fieldFusElContr", f"{100 * (1 - field['FusEl']['Jmin']):.0f}")
put("fieldSacElContr", f"{100 * (1 - field['SacEl']['Jmin']):.0f}")
put("fieldFusWdrop", f"{100 * (1 - field['FusEl']['Wapex'] / field['FusEl']['Whom']):.0f}")
put("fieldSacWdrop", f"{100 * (1 - field['SacEl']['Wapex'] / field['SacEl']['Whom']):.0f}")
put("fieldFusCcAdd", f"{field['FusCc']['WRmax'] - field['FusEl']['WRmax']:.0f}")
put("fieldSacCcAdd", f"{field['SacCc']['WRmax'] - field['SacEl']['WRmax']:.0f}")

def elem_centroids(feb):
    t = open(feb).read()
    nodes = {int(m.group(1)): tuple(map(float, m.group(2).split(",")))
             for m in re.finditer(r'<node id="\s*(\d+)">([^<]*)<', t)}
    els = {int(m.group(1)): [int(v) for v in m.group(2).split(",")]
           for m in re.finditer(r'<elem id="\s*(\d+)">([^<]*)<', t)}
    return {e: np.mean([nodes[n] for n in c], axis=0) for e, c in els.items()}


# (the saccular rate comparison is read at the apex since 28 sep: no probe angle or
# off-apex J to report; elem_centroids is kept for checks)

# homeostatic geometry of the fusiform segment (t = 1, end of Stage I)
nb = blocks(os.path.join(P, "matrix-fus-el-NODES.txt"))
rad = ref_radii(os.path.join(P, "matrix-fus-el.feb"))
h = at(nb, 1.0)
ids = np.array(sorted(rad))
Rr = np.array([rad[i] for i in ids])
rc = np.array([np.hypot(*h[i][:2]) for i in ids])
a_o = rc[Rr < Rr.min() + 0.001].mean()
b_o = rc[Rr > Rr.max() - 0.001].mean()
put("aneuAo", f"{a_o:.3f}")
put("aneuWo", f"{b_o - a_o:.2f}")

# ================================================================== rate comparison
F = ak("full-KG")
ONSET = 1.0


def rate_series(node_f, elem_f, s_load):
    nb, eb = blocks(node_f), blocks(elem_f)
    lam = {}
    for t, d in nb:
        x, ux, y, uy = next(iter(d.values()))[:4]
        lam[round(t, 4)] = np.hypot(x, y) / np.hypot(x - ux, y - uy)
    sv, J = {}, {}
    for t, d in eb:
        v = next(iter(d.values()))
        sv[round(t, 4)] = (v[0] + v[1] + v[2]) / 3.0
        J[round(t, 4)] = v[3]
    ts = [t for t in sorted(set(lam) & set(sv)) if t >= ONSET]
    lro = lam[ts[0]]
    x = np.array([(t - ONSET) / s_load for t in ts])
    k = x <= 5.0 + 1e-9
    return dict(x=x[k], sv=np.array([sv[t] for t in ts])[k],
                lr=np.array([lam[t] / lro for t in ts])[k],
                J=np.array([J[t] for t in ts])[k], svo=sv[ts[0]])


# rate comparison read at the apex (node 20, element 10) for both lesions: full-KG/apex/,
# copies of the rate decks with only the two primary loggers moved (28 sep); the earlier
# outputs (node 1209 / element 90, 81 deg from the saccular insult centre) stay in full-KG/
RATE = {
    "Fus": dict(mbe=("apex/node-axisym-mbe.txt", "apex/elem-axisym-mbe.txt"),
                one=("apex/ds-aneu/node-ds001.txt", "apex/ds-aneu/elem-ds001.txt"),
                hundred=("apex/node-axisym-sp100.txt", "apex/elem-axisym-sp100.txt")),
    "Sac": dict(mbe=("apex/node-axis-FINALmbe.txt", "apex/elem-axis-FINALmbe.txt"),
                one=("apex/ds-aneu/node-axis-FINALfk-SDS.txt", "apex/ds-aneu/elem-axis-FINALfk-SDS.txt"),
                hundred=("apex/node-axis-sp100.txt", "apex/elem-axis-sp100.txt")),
}
for m, cfg in RATE.items():
    mb = rate_series(os.path.join(F, cfg["mbe"][0]), os.path.join(F, cfg["mbe"][1]), 1.0)
    s1 = rate_series(os.path.join(F, cfg["one"][0]), os.path.join(F, cfg["one"][1]), 1.0)
    s100 = rate_series(os.path.join(F, cfg["hundred"][0]),
                       os.path.join(F, cfg["hundred"][1]), 100.0)
    Jm, lrm, svo = mb["J"][-1], mb["lr"][-1], mb["svo"]
    put(f"rate{m}J", f"{Jm:.3f}")
    put(f"rate{m}Lamr", f"{lrm:.3f}")
    put(f"rate{m}OneFrac", f"{100 * (s1['J'][-1] - 1) / (Jm - 1):.0f}")
    put(f"rate{m}OneOvershoot", pct(np.max(s1["sv"]) / svo - 1.0, 1))
    put(f"rate{m}OnePeakLamr", f"{np.max(s1['lr']):.3f}")
    put(f"rate{m}HundredSvDev", pct(np.max(abs(s100["sv"] / s100["svo"] - 1.0)), 1))
    bad = np.nonzero(abs(s100["J"] - Jm) / (Jm - 1.0) >= 0.05)[0]
    i = bad[-1] + 1 if len(bad) else 0
    put(f"rate{m}HundredSettle", f"{s100['x'][min(i, len(s100['x']) - 1)] * 100:.0f}")
    if m == "Fus":
        put("rateSvo", f"{svo:.0f}")
        J_mbe_fus = Jm
# sigma_vo of the calibrated material (fusiform deck at k_J = 1e6), element average. The
# original homeostatic state is uniform along the segment: element 1 (one end) and element 10
# (apex) carry the same value
put("aneuSvo", f"{kj[('aneu', '1e6')][1]:.0f}")
assert abs(kj_case(ak("kJ-study", "aneu_kJ1e6_matrix-fus-el-eB.txt"))[1] - kj[("aneu", "1e6")][1]) < 1e-6

# time-step refinement of the rate-dependent fusiform lesion, s_load = 1 day
frac = {}
for ds in ("01", "005", "001"):
    ts, v = one(blocks(os.path.join(F, "ds-aneu", f"elem-ds{ds}.txt")))
    frac[ds] = (np.interp(ONSET + 5.0, ts, v[:, 3]) - 1.0) / (J_mbe_fus - 1.0)
put("dsAneuCoarse", f"{100 * frac['01']:.0f}")
put("dsAneuMid", f"{100 * frac['005']:.0f}")
put("dsAneuFine", f"{100 * frac['001']:.0f}")
put("dsAneuRatio", f"{(frac['01'] - frac['001']) / (frac['005'] - frac['001']):.1f}")

# terminal state: equilibrated against rate-dependent (130 fixed steps), and the
# equilibrated solution at 20 steps (convergence-table run) against 522 steps
_, v = one(blocks(os.path.join(F, "ds-study", "elem-axisym-sp100_a_ds1.0.txt")))
J_full = v[-1, 3]
_, v = one(blocks(ak("timing", "conv-aneu-MBE_elem-axisym-mbe.txt")))
J_mbe20 = v[-1, 3]
put("termJ", f"{J_mbe_fus:.4f}")
put("termMbeFull", pct(abs(J_full - J_mbe_fus) / J_mbe_fus, 4))   # same lesion, apex element 10
put("termMbeSteps", pct(abs(J_mbe20 - J_mbe_fus) / J_mbe_fus, 3))
steps522 = int(re.search(r"Number of time steps completed[\. ]*: (\d+)",
                         open(os.path.join(F, "rate-fusiform-MBE.log")).read()).group(1))
put("termMbeManySteps", f"{steps522}")

# ============================================================ timing and tables
tim = {r["case"]: r for r in csv.DictReader(open(ak("timing", "timing.csv")))}
for c, name in (("uniax-MBE", "UniMbe"), ("uniax-full", "UniFull"), ("uniax-HYP", "UniHyp"),
                ("aneu-MBE", "AneuMbe"), ("aneu-full", "AneuFull"), ("aneu-HYP", "AneuHyp")):
    w = float(tim[c]["wall_s"])
    put(f"time{name}", f"{w:.3f}" if w < 1 else f"{w:.0f}")
    put(f"steps{name}", tim[c]["steps"])
    put(f"newton{name}", tim[c]["newton"])
a_m, a_f = tim["aneu-MBE"], tim["aneu-full"]
put("speedTime", f"{float(a_f['wall_s']) / float(a_m['wall_s']):.1f}")
put("speedSteps", f"{int(a_f['steps']) / int(a_m['steps']):.1f}")
put("speedNewton", f"{int(a_f['newton']) / int(a_m['newton']):.1f}")
put("itPerStepMbe", f"{int(a_m['newton']) / int(a_m['steps']):.1f}")
put("itPerStepFull", f"{int(a_f['newton']) / int(a_f['steps']):.1f}")
u_m, u_f = tim["uniax-MBE"], tim["uniax-full"]
put("speedUniTime", f"{float(u_f['wall_s']) / float(u_m['wall_s']):.1f}")

BLOCK = re.compile(
    r"Nonlinear solution status: time=\s*([\d.eE+-]+).*?"
    r"convergence norms :\s+INITIAL\s+CURRENT\s+REQUIRED\s*\n"
    r"\s*residual\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+[\d.eE+-]+\s*\n"
    r"\s*energy\s+([\d.eE+-]+)\s+([\d.eE+-]+)", re.S)


def representative(log, tmin, tmax):
    """Residual and energy norms of the step with the median Newton count among the
    steps in (tmin, tmax], the interval over which the load or the insult is applied;
    iteration 0 is the INITIAL norm of the step, iteration k the CURRENT norm after k
    iterations."""
    steps, cur, ct = [], [], None
    for m in BLOCK.finditer(open(log).read()):
        t = float(m.group(1))
        rec = tuple(float(m.group(i)) for i in (2, 3, 4, 5))
        if ct is None or t != ct:
            if cur:
                steps.append((ct, cur))
            cur, ct = [rec], t
        else:
            cur.append(rec)
    if cur:
        steps.append((ct, cur))
    st = sorted([s for s in steps if tmin < s[0] <= tmax + 1e-9], key=lambda s: len(s[1]))
    _, b = st[len(st) // 2]
    return [b[0][0]] + [x[1] for x in b], [b[0][2]] + [x[3] for x in b]


def table(runs, name):
    cols = [representative(ak("timing", log), t0, t1) for log, t0, t1 in runs]
    nrow = max(len(c[0]) for c in cols)
    lines = []
    for i in range(nrow):
        cells = []
        for fr, en in cols:
            cells += [sci(fr[i]), sci(en[i])] if i < len(fr) else ["--", "--"]
        # the layout of the original tables: label, then three pairs separated by spacers
        lines.append(f"Iteration {i} & {cells[0]} & {cells[1]} &  & {cells[2]} & {cells[3]} &  & "
                     f"{cells[4]} & {cells[5]} & \\\\")
    put(name, "\n".join(lines))
    return [len(c[0]) - 1 for c in cols]


# loading interval: uniaxial force ramp F_o -> 3 F_o over s_load = 10 days from t = 1;
# aneurysm insult over (1, 2] (MBE, degradation fraction) and (1, 101] (full-KG,
# s_load = 100 days); the hyperelastic aneurysm run is the pressurization (0, 1]
it_u = table([("conv-uniax-MBE.log", 1.0, 11.0), ("u_ds0.5.log", 1.0, 11.0),
              ("conv-uniax-HYP.log", 1.0, 11.0)], "convUniRows")
it_a = table([("conv-aneu-MBE.log", 1.0, 2.0), ("a_ds1.0.log", 1.0, 101.0),
              ("conv-aneu-HYP.log", 0.0, 1.0)], "convAneuRows")
for name, v in zip(("UniMbe", "UniFull", "UniHyp", "AneuMbe", "AneuFull", "AneuHyp"),
                   it_u + it_a):
    put(f"repIt{name}", f"{v}")

# tangent prefactor J sigma_vo / (K_sigma k Delta s) at the finest step of the full-KG run
put("prefactor", rf"{float(MACROS['rateSvo']) * J_mbe_fus / (0.1 * 0.001):.1e}"
    .replace("e+0", r"\times10^{") + "}")

# ============================================================== parametric sweep
# the 15 configurations shown in the paper (Online Resource 1, maps): nu = 1..6 at s_z = 1.0 mm
# and s_z = 1.0..2.8 mm at nu = 3 (= densify_runner._FOLDERS). The MASTER csv also holds other
# exploratory folders, never shown, whose saccular runs still probe element 90 (81 deg from the
# insult centre): they are excluded, so every sweep-wide number comes from the shown, apex-probed maps.
PAPER = {(float(n), 1.0) for n in range(1, 7)} | {(3.0, round(1.0 + 0.2 * i, 1)) for i in range(10)}
rows = [r for r in csv.DictReader(open(os.path.join(P, "MASTER_PARAMETRIC_STUDY.csv")))
        if r["alph"] == "ON" and "unit-frac" not in r["folder"]
        and (float(r["z3_nu"]), round(float(r["z4_zod"]), 1)) in PAPER]
GRID, JMAX = 400, 4.0          # admissible: converged and J <= 4 (no converged state has J < 1)


def conf(typ, nu, sz):
    return [r for r in rows if r["type"] == typ and float(r["z3_nu"]) == nu
            and abs(float(r["z4_zod"]) - sz) < 1e-9]


def adm(g):
    return [r for r in g if float(r["z2_J"]) <= JMAX]


def axes_limit(g):
    a = adm(g)
    de = max((float(r["x_elastin_deg"]) for r in a if float(r["y_collagen_deg"]) == 0.0), default=0)
    dc = max((float(r["y_collagen_deg"]) for r in a if float(r["x_elastin_deg"]) == 0.0), default=0)
    return de, dc


de1, dc1 = axes_limit(conf("axis", 3.0, 1.0))
de28, dc28 = axes_limit(conf("axis", 3.0, 2.8))
put("parFusDeOne", f"{100 * de1:.0f}")
_, dcs = axes_limit(conf("asym", 3.0, 1.0))
put("parSacDcOne", f"{100 * dcs:.0f}")
put("parFusLamrTrellis", f"{max(float(r['z1_lamb_r']) for z in (1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8) for r in adm(conf('axis', 3.0, z))):.1f}")
put("parFusDcOne", f"{100 * dc1:.0f}")
put("parFusDeLong", f"{100 * de28:.0f}")
put("parFusDcLong", f"{100 * dc28:.0f}")
for typ, name in (("axis", "Fus"), ("asym", "Sac")):
    g = adm(conf(typ, 3.0, 1.0))
    put(f"par{name}Adm", f"{100 * len(g) / GRID:.0f}")
    lr = np.array([float(r["z1_lamb_r"]) for r in g])
    put(f"par{name}LamrMax", f"{lr.max():.1f}")
    put(f"par{name}Clinical", f"{100 * np.mean(lr >= 1.5):.0f}")
    # linear fit of J over the admissible states of the whole sweep of this morphology
    ga = adm([r for r in rows if r["type"] == typ])
    A = np.array([[1.0, float(r["x_elastin_deg"]), float(r["y_collagen_deg"])] for r in ga])
    c = np.linalg.lstsq(A, np.array([float(r["z2_J"]) for r in ga]), rcond=None)[0]
    put(f"par{name}SensDe", f"{c[1]:.1f}")
    put(f"par{name}SensDc", f"{c[2]:.1f}")
    lr_all = [(float(r["z1_lamb_r"]), r) for r in ga]
    best = max(lr_all, key=lambda p: p[0])
    put(f"par{name}LamrSweep", f"{best[0]:.1f}")
    if typ == "asym":
        put("parSacLamrSweepNu", f"{float(best[1]['z3_nu']):.0f}")
        put("parSacLamrSweepSz", f"{float(best[1]['z4_zod']):.1f}")
ratios = []
for typ in ("axis", "asym"):
    ga = adm([r for r in rows if r["type"] == typ])
    A = np.array([[1.0, float(r["x_elastin_deg"]), float(r["y_collagen_deg"])] for r in ga])
    c = np.linalg.lstsq(A, np.array([float(r["z2_J"]) for r in ga]), rcond=None)[0]
    ratios.append(c[2] / c[1])
put("parSensRatio", f"{max(ratios):.1f}")

# beyond the admissible region: the csv holds only the converged grid points, so a missing
# point is a run that did not converge. At nu = 3, s_z = 1.0 mm: share of the grid converged
# with J > 4, share not converged, largest converged J. Over the 15 configurations: no
# non-converged point has an admissible neighbour on the grid (a band of converged J > 4
# states always separates the two), which the text states.
for typ, name in (("axis", "Fus"), ("asym", "Sac")):
    g = conf(typ, 3.0, 1.0)
    put(f"par{name}Beyond", f"{100 * sum(float(r['z2_J']) > JMAX for r in g) / GRID:.0f}")
    put(f"par{name}NoConv", f"{100 * (GRID - len(g)) / GRID:.0f}")
    put(f"par{name}JConvMax", f"{max(float(r['z2_J']) for r in g):.0f}")
_step = 1.0 / 19.0
for (nu, sz) in PAPER:
    for typ in ("axis", "asym"):
        pts = {(round(float(r["x_elastin_deg"]) / _step), round(float(r["y_collagen_deg"]) / _step)):
               float(r["z2_J"]) for r in conf(typ, nu, sz)}
        for i in range(20):
            for j in range(20):
                if (i, j) in pts:
                    continue
                nb = [pts.get(q) for q in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1))]
                if any(v is not None and v <= JMAX for v in nb):
                    raise SystemExit(f"non-converged point next to an admissible one: {typ} nu={nu} s_z={sz} {i},{j}")
# relative volume at the edge of convergence: J at the converged grid points that have a
# non-converged neighbour on the grid, over the 15 configurations (range and median), and
# the largest converged J of the whole sweep
for typ, name in (("axis", "Fus"), ("asym", "Sac")):
    edge, jall = [], []
    for (nu, sz) in PAPER:
        pts = {(round(float(r["x_elastin_deg"]) / _step), round(float(r["y_collagen_deg"]) / _step)):
               float(r["z2_J"]) for r in conf(typ, nu, sz)}
        jall += list(pts.values())
        edge += [v for (i, j), v in pts.items()
                 if any(q not in pts and 0 <= q[0] < 20 and 0 <= q[1] < 20
                        for q in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)))]
    assert min(edge) > JMAX
    put(f"par{name}JEdgeMin", f"{min(edge):.0f}" if min(edge) >= 5 else f"{min(edge):.1f}")
    put(f"par{name}JEdgeMed", f"{np.median(edge):.0f}")
    put(f"par{name}JEdgeMax", f"{max(edge):.0f}")
    put(f"par{name}JConvMaxAll", f"{max(jall):.0f}")

# sensitivity to the transition decay nu at s_z = 1.0 mm
xs = sorted({float(r["x_elastin_deg"]) for r in rows})
ys = sorted({float(r["y_collagen_deg"]) for r in rows})
DE = min(xs, key=lambda v: abs(v - 0.20))
DC = min(ys, key=lambda v: abs(v - 0.10))
put("nuDe", f"{100 * DE:.0f}\\%")   # in percent, as the maps and Section 5.5
put("nuDc", f"{100 * DC:.0f}\\%")
for typ, name in (("axis", "Fus"), ("asym", "Sac")):
    Jn, fr = [], []
    for nu in range(1, 7):
        g = conf(typ, float(nu), 1.0)
        hit = [r for r in g if abs(float(r["x_elastin_deg"]) - DE) < 1e-9
               and abs(float(r["y_collagen_deg"]) - DC) < 1e-9]
        Jn.append(float(hit[0]["z2_J"]))
        fr.append(100.0 * len(adm(g)) / GRID)
    put(f"nu{name}JOne", f"{Jn[0]:.2f}")
    put(f"nu{name}JSix", f"{Jn[-1]:.2f}")
    put(f"nu{name}Spread", f"{100 * (Jn[-1] - Jn[0]) / Jn[0]:.0f}")
    put(f"nu{name}IncTail", f"{max(abs(np.diff(Jn[3:]))):.2f}")
    put(f"nu{name}AdmOne", f"{fr[0]:.0f}")
    put(f"nu{name}AdmSix", f"{fr[-1]:.0f}")
    if typ == "axis":
        put("nuFusJThree", f"{Jn[2]:.2f}")
        put("nuFusJFive", f"{Jn[4]:.2f}")
        put("nuFusThreeFive", pct((Jn[4] - Jn[2]) / Jn[2], 1))

# ======================================================================== output
with _open(os.path.join(HERE, "numbers.tex"), "w") as f:
    f.write("% generated by make_numbers.py; do not edit\n")
    for k in sorted(MACROS):
        f.write(f"\\newcommand{{\\{k}}}{{{MACROS[k]}}}\n")
with _open(os.path.join(HERE, "numbers-inputs.txt"), "w") as f:
    f.write("\n".join(sorted(INPUTS)) + "\n")
print(f"numbers.tex: {len(MACROS)} macros from {len(INPUTS)} files")
for k in sorted(MACROS):
    print(f"  {k:24s} {MACROS[k]}")
