"""Single generator of every result quoted in main.tex. Writes numbers.tex and the bodies
of the two convergence tables (macros convUniRows and convAneuRows).

Sources, all under ../aneurysm-knockout/ (AK). Since 2026-10-02 (night), every run that enters
Stage II without the fibre-angle option is read from AK/alp1/ (A1): the same input files with
<alp_update>1</alp_update>, run with the corrected code (~/FEBio-alpcheck: angle from F Fo^-1, two
tangent terms), in folders that mirror the ones listed here; runs carried to 5 s_load in */to5/, timed
runs ending at the stopping time in A1/speedup-fair/stop/ and A1/uniax-lo1/timing-stop/, the uniaxial
plot series in A1/uniax-lo1/prep/. The lesion fields, the sweep, k_J and mesh of the lesions and the
passive runs (no Stage II) are read from their original folders.
  AK/uniax-lo1/              uniaxial coupon, runs with lo_opt = 1 (the passive runs are read
                             from AK/uniax-fractions/): FEBio logs (elem1_*, node7_*), the plot-ready
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
import json
import os
import re
import math

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
AK = os.path.join(ROOT, "aneurysm-knockout")
# The uniaxial inputs are read from AK/uniax-lo1/ (runs with lo_opt = 1: fibre stress evaluated at
# the original homeostatic state, the rule of the lesion runs; adopted by the author on
# 2026-10-02). UNIAX_LO0=1 reads the earlier runs (lo_opt = 0, AK/uniax-fractions/) instead and
# writes numbers.lo0.tex, leaving numbers.tex alone.
LO1 = os.environ.get("UNIAX_LO0") != "1"
MACROS = {}

_open, INPUTS = builtins.open, set()


# s1/ (2026-10-03): all runs repeated with Stage I as in Section 5.1 (F_g = 1, F_e = F, S = J S_e + penalty;
# ~/FEBio-stage1) in AK/s1/, same relative paths. A read under AK goes to AK/s1/ when that copy exists; only the
# parametric maps and the unused lad*/lat* block keep the earlier runs (checked at the end).
S1 = os.path.join(AK, "s1")
S1_KEEP = ("severity-fig/ladder-numbers.json",
           "tancheck/")   # tangent check against finite differences (its own FEBio copy, no s1 counterpart)
OUTSIDE_S1 = set()


def _in_s1(path):
    q = os.path.abspath(path)
    if q.startswith(AK + os.sep) and not q.startswith(S1 + os.sep):
        r = os.path.join(S1, os.path.relpath(q, AK))
        if os.path.exists(r):
            return r
        OUTSIDE_S1.add(os.path.relpath(q, AK))
    return path


def _recording_open(file, mode="r", *a, **k):
    if isinstance(file, str) and "r" in mode and "+" not in mode:
        file = _in_s1(file)
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


A1 = ak("alp1")
A2 = ak("alp2")      # rate comparison with fus cc_in 0.95, sac cc_in 0.9 (2026-10-03)


# =================================================================== uniaxial coupon
U = os.path.join(A1, "uniax-lo1") if LO1 else ak("uniax-fractions")
U_PREP = os.path.join(U, "prep") if LO1 else U      # plot series {MBE,sp1,sp10,sp100}_{c,a}.txt
U_PASS = ak("uniax-fractions")              # passive runs: no growth, unaffected by lo_opt
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
    src = U_PASS if tag.startswith("PASS") else U
    tb, eb = one(blocks(os.path.join(src, f"elem1_{tag}.txt")))
    tnb, nb = one(blocks(os.path.join(src, f"node7_{tag}.txt")))
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
    path = os.path.join(U_PREP, f"{tag}_{kind}.txt")
    path = _in_s1(path)
    INPUTS.add(os.path.relpath(path, HERE))         # np.loadtxt does not pass through open()
    d = np.loadtxt(path)   # x, sigma, lam/lam_o, lam_g
    return d[d[:, 0] >= 0.0]


fe = {t: series(t, "c") for t in ("MBE", "sp1", "sp10", "sp100")}
sa = {t: series(t, "a") for t in ("MBE", "sp1", "sp10", "sp100")}
SLOAD = {"MBE": 10.0, "sp1": 1.0, "sp10": 10.0, "sp100": 100.0}
# Section 6.1.1 (2026-10-03, asked by the session of the simulated review, sentence approved by the author): largest
# axial Cauchy stress of the finite-element curve during the ramp, 0 <= s/s_load <= 1, over sigma_o, in %
for _t, _n in (("sp1", "One"), ("sp10", "Ten"), ("sp100", "Hundred")):
    _r = fe[_t][fe[_t][:, 0] <= 1.0 + 1e-9]
    put(f"uniRampSig{_n}", f"{100.0 * (_r[:, 1].max() / sig_o - 1.0):.0f}")
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
put("uniFeSemiSig", pct(dsig, 2))
put("uniFeSemiLam", pct(dlam, 2))


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
# the one-day step is the run u_ds1.0t1 (2026-10-02, must point at t = 1): in u_ds1.0 the first step
# was cut to 0.6 day, Stage I ended at 0.6 F_o and J came out about twice too large
DS = ["0.005", "0.01", "0.02", "0.05", "0.1", "0.2", "0.5"] + (["1.0t1"] if LO1 else [])
Jds = {}
for ds in DS:
    tb, eb = one(blocks(os.path.join(U, "ds-study", f"elem1_ds{ds}.txt")))
    Jds[ds] = eb[-1, 3]
    assert not LO1 or abs(tb[np.argmin(abs(tb - 1.0))] - 1.0) < 1e-9, ds        # a step ends at t = 1
put("uniDsChange", pct(abs(Jds["0.5"] - Jds["0.1"]) / Jds["0.1"], 2))    # two decimals: 0.05 %
if LO1:
    put("uniDsOneChange", pct(abs(Jds["1.0t1"] - Jds["0.1"]) / Jds["0.1"], 2))
    put("uniDsOneSteps", re.search(r"Number of time steps completed[\. ]*: (\d+)",
                                   open(os.path.join(U, "ds-study", "u_ds1.0t1.log"), errors="replace").read()).group(1))

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
    kj[("uniax", k)] = kj_case(os.path.join(U if LO1 else AK, "kJ-study", f"uniax_kJ{k}_elem1_MBE.txt"))
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


ZREF = {}      # reference axial coordinate of each node, per deck (filled by ref_radii)


def ref_radii(feb):
    rad, inn = {}, False
    zref = ZREF.setdefault(feb, {})
    for line in open(feb):
        s = line.strip()
        if s.startswith("<Nodes"):
            inn = True
            continue
        if s.startswith("</Nodes>"):
            break
        if inn and s.startswith("<node"):
            nid = int(s.split('id="')[1].split('"')[0])
            x, y, z = [float(v) for v in s.split(">", 1)[1].split("<")[0].split(",")]
            rad[nid] = np.hypot(x, y)
            zref[nid] = round(z, 6)
    return rad


def width(rad, zref, pts, wall):
    """Diameter of the outer or inner wall: the largest distance between two nodes of that wall
    in a cross-section (nodes with the same reference axial coordinate), over all cross-sections.
    Until 2026-10-02 the diameter was twice the largest radial position of an outer node, which
    in the saccular lesion (one side bulges) is larger than the width of the section."""
    ids = sorted(rad)
    Rr = np.array([rad[i] for i in ids])
    sel = [i for i, r in zip(ids, Rr) if (r > Rr.max() - 0.01 if wall == "outer" else r < Rr.min() + 1e-3)]
    best = 0.0
    for z in sorted({zref[i] for i in sel}):
        q = np.array([pts[i][:2] for i in sel if zref[i] == z])
        best = max(best, np.sqrt(((q[:, None, :] - q[None, :, :]) ** 2).sum(-1)).max())
    return best


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
        dmax=width(rad, ZREF[os.path.join(P, f"{tag}.feb")], nb[-1][1], "outer"),
        dhom=width(rad, ZREF[os.path.join(P, f"{tag}.feb")], at(nb, 1.0), "outer"),
        dratio=(width(rad, ZREF[os.path.join(P, f"{tag}.feb")], nb[-1][1], "inner")
                / width(rad, ZREF[os.path.join(P, f"{tag}.feb")], at(nb, 1.0), "inner")))
for key, f in field.items():
    put(f"field{key}Jmax", f"{f['Jmax']:.2f}")
    put(f"field{key}Jmin", f"{f['Jmin']:.2f}")
    put(f"field{key}Lamr", f"{f['lamr']:.2f}")
    put(f"field{key}Diam", f"{f['dmax']:.1f}")
    put(f"field{key}DiamRatio", f"{f['dratio']:.2f}")
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
# the text gives the range of the inner-wall diameter ratio as FusEl ... FusCc
assert min(f["dratio"] for f in field.values()) == field["FusEl"]["dratio"]
assert max(f["dratio"] for f in field.values()) == field["FusCc"]["dratio"]
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
# undeformed geometry of the mesh (the in-vivo geometry of [58]) and the size of F_o in the lesion material
# (Table 1): every nodal radius changes by less than 1e-4 in Stage I, and the strain energy and J at the end of
# Stage I are the same in all elements (Section 5.4; round rABC-6, marks 59-74)
put("aneuAzero", f"{Rr.min():.3f}")
put("aneuWzero", f"{Rr.max() - Rr.min():.2f}")
# s1/ (2026-10-03): with Stage I as Section 5.1 the largest change is 1.07e-4 (it was below 1e-4 before);
# the bound is the measured maximum rounded up to two significant digits
_dr = float(np.max(abs(rc / Rr - 1.0)))
_e = int(np.floor(np.log10(_dr)))
_m = np.ceil(_dr / 10 ** (_e - 1)) / 10
put("aneuFoBound", f"{_m:.1f}\\times10^{{{_e}}}" if _m != 1 else f"10^{{{_e}}}")
_hom = at(blocks(os.path.join(P, "matrix-fus-el-ALL.txt")), 1.0)
for _k in (0, 1):                                   # J and W_R of every element at t = 1
    _v = np.array([v[_k] for v in _hom.values()])
    assert _v.max() / _v.min() - 1.0 < 1e-6, (_k, _v.min(), _v.max())

# ================================================================== rate comparison
F = os.path.join(A2, "full-KG")
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


# the s_load = 1 day runs are those of apex/ds-aneu-t1/ (2026-10-02, must point at t = 1: with a step
# of 0.001 day the earlier runs, in apex/ds-aneu/, stored sigma_vo one step before the end of Stage I)
# rate comparison read at the apex (node 20, element 10) for both lesions: full-KG/apex/,
# copies of the rate decks with only the two primary loggers moved (28 sep); the earlier
# outputs (node 1209 / element 90, 81 deg from the saccular insult centre) stay in full-KG/
RATE = {
    "Fus": dict(mbe=("apex/node-axisym-mbe.txt", "apex/elem-axisym-mbe.txt"),
                one=("apex/ds-aneu-t1/node-ds001.txt", "apex/ds-aneu-t1/elem-ds001.txt"),
                hundred=("apex/to5/node-axisym-sp100.txt", "apex/to5/elem-axisym-sp100.txt")),
    "Sac": dict(mbe=("apex/node-axis-FINALmbe.txt", "apex/elem-axis-FINALmbe.txt"),
                one=("apex/ds-aneu-t1/node-axis-FINALfk-SDS.txt", "apex/ds-aneu-t1/elem-axis-FINALfk-SDS.txt"),
                hundred=("apex/to5/node-axis-sp100.txt", "apex/to5/elem-axis-sp100.txt")),
}
# G&R time, in days, at which three rate-dependent runs end before s = 5 s_load
def _end_days(elem_f):
    return blocks(os.path.join(F, elem_f))[-1][0] - ONSET


put("rateFusTenEnd", f"{_end_days('apex/to5/elem-axisym-sp10.txt'):.0f}")
put("rateFusHundredEnd", f"{_end_days('apex/to5/elem-axisym-sp100.txt'):.0f}")
put("rateSacHundredEnd", f"{_end_days('apex/to5/elem-axis-sp100.txt'):.0f}")
assert _end_days("apex/elem-axis-sp10.txt") >= 50.0          # the saccular 10-day run reaches 5 s_load
assert _end_days("apex/ds-aneu-t1/elem-ds001.txt") >= 5.0 - 1e-9 and _end_days("apex/ds-aneu-t1/elem-axis-FINALfk-SDS.txt") >= 5.0 - 1e-9
for m, cfg in RATE.items():
    mb = rate_series(os.path.join(F, cfg["mbe"][0]), os.path.join(F, cfg["mbe"][1]), 1.0)
    s1 = rate_series(os.path.join(F, cfg["one"][0]), os.path.join(F, cfg["one"][1]), 1.0)
    s100 = rate_series(os.path.join(F, cfg["hundred"][0]),
                       os.path.join(F, cfg["hundred"][1]), 100.0)
    Jm, lrm, svo = mb["J"][-1], mb["lr"][-1], mb["svo"]
    if m == "Fus":                                  # radial stretch of the apex in Stage I (rate material)
        _d = at(blocks(os.path.join(F, cfg["mbe"][0])), 1.0)[20]
        put("rateLamRoPct", pct(np.hypot(_d[0], _d[2]) / np.hypot(_d[0] - _d[1], _d[2] - _d[3]) - 1.0, 1))
    put(f"rate{m}J", f"{Jm:.3f}")
    put(f"rate{m}Lamr", f"{lrm:.3f}")
    put(f"rate{m}OneFrac", f"{100 * (s1['J'][-1] - 1) / (Jm - 1):.0f}")
    put(f"rate{m}OneOvershoot", pct(np.max(s1["sv"]) / svo - 1.0, 1))
    put(f"rate{m}OneSvEnd", pct(s1["sv"][-1] / svo - 1.0, 1))
    _ten = {"Fus": ("apex/to5/node-axisym-sp10.txt", "apex/to5/elem-axisym-sp10.txt"),
            "Sac": ("apex/node-axis-sp10.txt", "apex/elem-axis-sp10.txt")}[m]
    s10 = rate_series(os.path.join(F, _ten[0]), os.path.join(F, _ten[1]), 10.0)
    assert abs(s10["x"][-1] - 5.0) < 1e-6 and abs(s100["x"][-1] - 5.0) < 1e-6
    assert abs(s10["sv"][-1] / svo - 1.0) < 0.05 and abs(s100["sv"][-1] / svo - 1.0) < 0.05
    put(f"rate{m}OnePeakLamr", f"{np.max(s1['lr']):.3f}")
    put(f"rate{m}HundredSvDev", pct(np.max(abs(s100["sv"] / s100["svo"] - 1.0)), 1))
    bad = np.nonzero(abs(s100["J"] - Jm) / (Jm - 1.0) >= 0.05)[0]
    i = bad[-1] + 1 if len(bad) else 0
    put(f"rate{m}HundredSettle", f"{s100['x'][min(i, len(s100['x']) - 1)] * 100:.0f}")
    if m == "Fus":
        put("rateSvo", f"{svo:.0f}")
        J_mbe_fus = Jm
    else:
        J_mbe_sac = Jm
# elastin-only lesions in the rate-comparison material: alp1/ (cc_in 1.0), apex radius increase at the end;
# Figure 4(a) and (c) for comparison (rABC-6 marks 77-78, last paragraph of Section 5.4)
for m, cfg in RATE.items():
    _e = rate_series(os.path.join(A1, "full-KG", cfg["mbe"][0]), os.path.join(A1, "full-KG", cfg["mbe"][1]), 1.0)
    put(f"rateEl{m}LamrPct", pct(_e["lr"][-1] - 1.0, 1))
put("fieldFusElLamrPct", pct(field["FusEl"]["lamr"] - 1.0, 0))
put("fieldSacElLamrPct", pct(field["SacEl"]["lamr"] - 1.0, 0))
# sigma_vo of the calibrated material (fusiform deck at k_J = 1e6), element average. The
# original homeostatic state is uniform along the segment: element 1 (one end) and element 10
# (apex) carry the same value
put("aneuSvo", f"{kj[('aneu', '1e6')][1]:.0f}")
assert abs(kj_case(ak("kJ-study", "aneu_kJ1e6_matrix-fus-el-eB.txt"))[1] - kj[("aneu", "1e6")][1]) < 1e-6

# time-step study of the rate-dependent lesions, s_load = 1 day: the runs of apex/ds-aneu-t1/
# (2026-10-02, must point at t = 1). The earlier study, in ds-aneu/, had no must point; with its
# three steps Stage II began one step early, and the 139, 102 and 75 % it gave were that
# defect and not a time-integration error
DT1 = os.path.join(F, "apex", "ds-aneu-t1")
dsJ = {ds: one(blocks(os.path.join(DT1, f"elem-ds{ds}.txt"))) for ds in ("01", "005", "001")}
frac = {ds: (np.interp(ONSET + 5.0, ts, v[:, 3]) - 1.0) / (J_mbe_fus - 1.0) for ds, (ts, v) in dsJ.items()}
assert len({f"{100 * f:.0f}" for f in frac.values()}) == 1, frac     # the three steps give the same integer
put("dsAneuFrac", f"{100 * frac['001']:.0f}")
assert MACROS["dsAneuFrac"] == MACROS["rateFusOneFrac"]
# largest difference in J over the computation, relative to the final volume increase
t_ref, v_ref = dsJ["001"]
inc = np.interp(ONSET + 5.0, t_ref, v_ref[:, 3]) - 1.0
err = {ds: np.max(abs(dsJ[ds][1][:, 3] - np.interp(dsJ[ds][0], t_ref, v_ref[:, 3]))) / inc for ds in ("01", "005")}
put("dsAneuCoarseErr", f"{100 * err['01']:.2g}")      # two significant figures: 0.15 and 0.068,
put("dsAneuMidErr", f"{100 * err['005']:.2g}")        # so that their ratio reads as dsAneuRatio
put("dsAneuRatio", f"{err['01'] / err['005']:.1f}")
sacJ = {tag: one(blocks(os.path.join(DT1, f"elem-axis-FINALfk-SDS{tag}.txt"))) for tag in ("01", "")}
fsac = {tag: (np.interp(ONSET + 5.0, ts, v[:, 3]) - 1.0) / (J_mbe_sac - 1.0) for tag, (ts, v) in sacJ.items()}
assert len({f"{100 * f:.0f}" for f in fsac.values()}) == 1, fsac     # the two steps give the same integer
put("dsSacFrac", f"{100 * fsac['']:.0f}")
assert MACROS["dsSacFrac"] == MACROS["rateSacOneFrac"]

# terminal state: equilibrated against rate-dependent (the timed run, Delta s = 5 days to
# s = 130 days), and the timed equilibrated solution (one pressurization step and ten
# insult increments) against 522 steps. Both timed runs: speedup-fair/ (2026-10-02), apex element 10
_, v = one(blocks(os.path.join(A2, "speedup-fair", "stop", "T-full-ds5-r1_elem-fus-sl100-2st-ds5.txt")))
J_full = v[-1, 3]
_, v = one(blocks(os.path.join(A2, "speedup-fair", "stop", "T-mbe-1p10-r1_elem-fus-sl100-2st-ds5.txt")))
J_mbe20 = v[-1, 3]
put("termJ", f"{J_mbe_fus:.4f}")
put("termMbeFull", pct(abs(J_full - J_mbe_fus) / J_mbe_fus, 4))   # same lesion, apex element 10
put("termMbeSteps", pct(abs(J_mbe20 - J_mbe_fus) / J_mbe_fus, 3))
steps522 = int(re.search(r"Number of time steps completed[\. ]*: (\d+)",
                         open(os.path.join(F, "apex", "rate-fusiform-MBE.log")).read()).group(1))
put("termMbeManySteps", f"{steps522}")

# ============================================================ timing and tables
# AK/timing/timing.csv only for UNIAX_LO0=1: with LO1 its three uniaxial cases are replaced below by
# timing-stop/ and its aneurysm cases are not used (s1/, 2026-10-03)
tim = {} if LO1 else {r["case"]: r for r in csv.DictReader(open(ak("timing", "timing.csv")))}
if LO1:
    # uniax-lo1/timing-t11/ (2026-10-02): the equilibrated and hyperelastic runs end at t = 11, when
    # the load reaches 3 F_o, and carry a must point at t = 1. In uniax-lo1/timing/ they ran 40 more
    # steps at constant load and ended Stage I at t = 0.8627 (sigma_o 222.9 instead of 260.6 kPa)
    tim.update({r["case"]: r for r in csv.DictReader(open(os.path.join(U, "timing-stop", "timing.csv")))})
for c, name in (("uniax-MBE", "UniMbe"), ("uniax-full", "UniFull"), ("uniax-HYP", "UniHyp")):
    w = float(tim[c]["wall_s"])
    put(f"time{name}", f"{w:.3f}" if w < 1 else f"{w:.0f}")
    put(f"steps{name}", tim[c]["steps"])
    put(f"newton{name}", tim[c]["newton"])
# aneurysm timing, measured again on 2026-10-02 with the same pressurization in both models
# (one step) and the full-KG run at Delta s = 5 days, which reproduces the equilibrated J
# (speedup-fair/timing_fair.csv; second batch timing_fair2.csv: Delta s = 1 day over the same
# 130 days and the hyperelastic pressurization in one step, as in the two growth models)
# A2/speedup-fair/stop/timing_stop.csv (2026-10-03): the full-KG runs end at the stopping time
# s = 145 days (J at the apex changes by less than 1e-4 (J - 1) over the last s_load/20 = 5 days, after
# the end of the load), with Delta s = 5 and 1 day, and mbe-1p10 timed in the same session;
# A2/speedup-fair/timing_fair2.csv gives the hyperelastic case and mbe-1p10-b (speed check)
fair = {r["case"]: r for r in csv.DictReader(open(os.path.join(A2, "speedup-fair", "stop", "timing_stop.csv")))}
fair2 = {r["case"]: r for r in csv.DictReader(open(os.path.join(A2, "speedup-fair", "timing_fair2.csv")))}
# the two batches ran at the same speed: mbe-1p10 was timed in both
assert abs(float(fair2["mbe-1p10-b"]["wall_min_s"]) / float(fair["mbe-1p10"]["wall_min_s"]) - 1.0) < 0.03
for c, name, src in (("mbe-1p10", "AneuMbe", fair), ("full-ds5", "AneuFull", fair),
                     ("hyp-1", "AneuHyp", fair2), ("full-ds1", "AneuFullFine", fair)):
    tim["aneu-" + name] = dict(wall_s=src[c]["wall_min_s"], steps=src[c]["steps"], newton=src[c]["newton"])
    put(f"time{name}", f"{float(src[c]['wall_min_s']):.0f}")
    put(f"steps{name}", src[c]["steps"])
    put(f"newton{name}", src[c]["newton"])
assert fair["full-ds5"]["per_stage_steps_iters"].startswith("1/") and fair["mbe-1p10"]["per_stage_steps_iters"].startswith("1/")
a_m, a_f = tim["aneu-AneuMbe"], tim["aneu-AneuFull"]
put("speedTime", f"{float(a_f['wall_s']) / float(a_m['wall_s']):.1f}")
put("speedSteps", f"{int(a_f['steps']) / int(a_m['steps']):.1f}")
put("speedNewton", f"{int(a_f['newton']) / int(a_m['newton']):.1f}")
put("itPerStepMbe", f"{int(a_m['newton']) / int(a_m['steps']):.1f}")
put("itPerStepFull", f"{int(a_f['newton']) / int(a_f['steps']):.1f}")

# stopping criterion of the timed rate-dependent runs (Section 5.3, 2026-10-02 night): the first G&R time
# s >= s_load at which J has changed by less than STOP_TOL (J - 1) over the last s_load/20 (J of the apex
# element in the aneurysm, of the element in the uniaxial test). The stopping time is found on untimed
# runs carried beyond it (A2/speedup-fair/stop/find/, A1/uniax-lo1/timing-stop/find/); the timed runs
# must end at that time.
STOP_TOL = 1e-4


def stop_time(elem_f, s_load):
    ts, v = one(blocks(elem_f))
    J = {round(t, 6): x[3] for t, x in zip(ts, v)}
    for t in sorted(J):
        prev = round(t - s_load / 20.0, 6)
        if t - ONSET >= s_load - 1e-9 and prev in J and abs(J[t] - J[prev]) < STOP_TOL * (J[t] - 1.0):
            return t - ONSET
    raise SystemExit(f"{elem_f}: stopping criterion not met")


stop_aneu = stop_time(os.path.join(A2, "speedup-fair", "stop", "find", "elem-fus-sl100-2st-ds5.txt"), 100.0)
stop_uni = stop_time(os.path.join(U, "timing-stop", "find", "elem1_ds0.5.txt"), 10.0)
for f in ("T-full-ds5-r1_elem-fus-sl100-2st-ds5.txt", "T-full-ds1-r1_elem-fus-sl100-2st-ds1.txt"):
    assert abs(blocks(os.path.join(A2, "speedup-fair", "stop", f))[-1][0] - ONSET - stop_aneu) < 1e-9, f
assert abs(blocks(os.path.join(U, "timing-stop", "elem1_ds0.5.txt"))[-1][0] - T_O - stop_uni) < 1e-9
put("stopTol", "10^{-4}")
put("stopAneu", f"{stop_aneu:.0f}")
put("stopUni", f"{stop_uni:.1f}")
# no wall-time ratio for the uniaxial test: with one element the fixed cost of a run dominates the
# wall time (author, 2026-10-02); the text gives steps, Newton iterations and the two wall times

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
    folder = os.path.join(U, "timing-stop") if LO1 and name == "convUniRows" else os.path.join(A2, "timing")
    cols = [representative(os.path.join(folder, log), t0, t1) for log, t0, t1 in runs]
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
              ("../speedup-fair/T-hyp-1-r1.log", 0.0, 1.0)], "convAneuRows")
for name, v in zip(("UniMbe", "UniFull", "UniHyp", "AneuMbe", "AneuFull", "AneuHyp"),
                   it_u + it_a):
    put(f"repIt{name}", f"{v}")

# saccular rate comparison (rABC-7, 2026-10-04): A2/speedup-fair-sac/, decks of setup_sac.py, the lesion of
# Figure 8 (65 % loss of elastic-fiber integrity, 10 % loss of collagen cross-linking). Steps, Newton
# iterations, apex J and residuals do not depend on the load of the machine; the wall times of that batch
# are not used (the parametric maps were running: 1.8 s per iteration against 0.45 s in the fusiform runs)
SAC = os.path.join(A2, "speedup-fair-sac")
sac_t = {r["case"]: r for r in csv.DictReader(open(os.path.join(SAC, "loaded", "timing_sac.csv")))}
for c, name in (("mbe-1p10-sac", "SacMbe"), ("full-ds5-sac", "SacFull")):
    put(f"steps{name}", sac_t[c]["steps"])
    put(f"newton{name}", sac_t[c]["newton"])
# speed-up in steps and Newton iterations, rounded as speedSteps / speedNewton (rABC-8, paper-review session)
put("speedStepsSac", f"{int(sac_t['full-ds5-sac']['steps']) / int(sac_t['mbe-1p10-sac']['steps']):.1f}")
put("speedNewtonSac", f"{int(sac_t['full-ds5-sac']['newton']) / int(sac_t['mbe-1p10-sac']['newton']):.1f}")
# wall times of the saccular rate comparison from the batch run on the idle machine (clean/, 3 repetitions, one
# thread; rABC-8, paper-review session), rounded as timeAneuMbe/Full and speedTime; steps and Newton iterations
# must equal those of loaded/
sac_c = {r["case"]: r for r in csv.DictReader(open(os.path.join(SAC, "clean", "timing_sac.csv")))}
for c in ("mbe-1p10-sac", "full-ds5-sac"):
    assert (sac_c[c]["steps"], sac_c[c]["newton"]) == (sac_t[c]["steps"], sac_t[c]["newton"]), c
put("timeSacMbe", f"{float(sac_c['mbe-1p10-sac']['wall_min_s']):.0f}")
put("timeSacFull", f"{float(sac_c['full-ds5-sac']['wall_min_s']):.0f}")
put("speedTimeSac", f"{float(sac_c['full-ds5-sac']['wall_min_s']) / float(sac_c['mbe-1p10-sac']['wall_min_s']):.1f}")
# mean speed-up of the two lesions (abstract; rABC-8, paper-review session): mean of the unrounded ratios
put("speedStepsMean", f"{(int(a_f['steps']) / int(a_m['steps']) + int(sac_t['full-ds5-sac']['steps']) / int(sac_t['mbe-1p10-sac']['steps'])) / 2:.1f}")
put("speedTimeMean", f"{(float(a_f['wall_s']) / float(a_m['wall_s']) + float(sac_c['full-ds5-sac']['wall_min_s']) / float(sac_c['mbe-1p10-sac']['wall_min_s'])) / 2:.1f}")

# tangent check (AK/tancheck/, rABC-8): relative difference between the code tangent and central finite
# differences of the stress, lines "FDCHK group total ... rel_to_total=<x>"; largest over every line of the
# uniaxial logs and of the lesion logs, rounded UP to one significant digit (the text says "below")
def fd_max(logs):
    vals = [float(re.search(r"rel_to_total=(\S+)", l).group(1))
            for f in logs for l in open(os.path.join(AK, "tancheck", f)) if l.startswith("FDCHK group total")]
    assert vals, logs
    return max(vals)


def up1(x):
    e = math.floor(math.log10(x))
    m = math.ceil(round(x / 10 ** e, 9))
    if m == 10:
        m, e = 1, e + 1
    return rf"{m}\times10^{{{e}}}"


put("fdUniErr", up1(fd_max(["a-full-t5/u_ds0.5.log", "a-full-t20/u_ds0.5.log", "a-full-t5-alp0/u_ds0.5.log",
                            "b-mbe-t5/conv-uniax-MBE.log"])))
put("fdLesionErr", up1(fd_max(["c-lesion-t10/a_ds1.0.log", "d-sac-full-t11/sac-full.log",
                               "e-sac-mbe-t1.5/sac-mbe.log", "f-fus-mbe-t1.5/fus-mbe.log"])))
stop_sac = stop_time(os.path.join(SAC, "find", "elem-sac-sl100-2st-ds5.txt"), 100.0)
_, v = one(blocks(os.path.join(SAC, "loaded", "T-full-ds5-sac-r1_elem-sac-sl100-2st-ds5.txt")))
J_full_sac = v[-1, 3]
assert abs(blocks(os.path.join(SAC, "loaded", "T-full-ds5-sac-r1_elem-sac-sl100-2st-ds5.txt"))[-1][0]
           - ONSET - stop_sac) < 1e-9
_, v = one(blocks(os.path.join(SAC, "loaded", "T-mbe-1p10-sac-r1_elem-sac-sl100-2st-ds5.txt")))
J_mbe_sac_t = v[-1, 3]
put("stopSac", f"{stop_sac:.0f}")
put("termSacJ", f"{J_mbe_sac_t:.4f}")
put("termSacMbeFull", pct(abs(J_full_sac - J_mbe_sac_t) / J_mbe_sac_t, 4))   # apex element 10
assert f"{J_mbe_sac_t:.3f}" == MACROS["rateSacJ"], (J_mbe_sac_t, MACROS["rateSacJ"])   # the Figure 8 lesion


def tex_sci(x):
    m, e = f"{x:.1E}".split("E")
    return m if int(e) == 0 else f"{m}\\times10^{{{int(e)}}}"


fr_sac, _ = representative(os.path.join(SAC, "loaded", "T-mbe-1p10-sac-r1.log"), 1.0, 2.0)
put("convSacRes", ", ".join(tex_sci(x) for x in fr_sac))
# saccular convergence table (rABC-8, paper-review session): the median-iteration step of the MBE-KG run over
# (1, 2] and of the full-KG run over (1, 101], same cells as convAneuRows, two columns
sac_cols = [representative(os.path.join(SAC, "loaded", "T-mbe-1p10-sac-r1.log"), 1.0, 2.0),
            representative(os.path.join(SAC, "loaded", "T-full-ds5-sac-r1.log"), 1.0, 101.0)]
_rows = []
for i in range(max(len(c[0]) for c in sac_cols)):
    cells = []
    for fr, en in sac_cols:
        cells += [sci(fr[i]), sci(en[i])] if i < len(fr) else ["--", "--"]
    _rows.append(f"Iteration {i} & {cells[0]} & {cells[1]} &  & {cells[2]} & {cells[3]} \\\\")
put("convSacRows", "\n".join(_rows))
put("repItSacMbe", f"{len(sac_cols[0][0]) - 1}")
put("repItSacFull", f"{len(sac_cols[1][0]) - 1}")

# tangent prefactor J sigma_vo / (K_sigma k Delta s) at the finest step of the full-KG run
put("prefactor", rf"{float(MACROS['rateSvo']) * J_mbe_fus / (0.1 * 0.001):.1e}"
    .replace("e+0", r"\times10^{") + "}")

# axial extent of the defect relative to the segment length (rABC-8, paper-review session): l_0 = axial length
# of the mesh (max z - min z); s_z of the Figure 4 lesion decks, and the s_z range of the paper profiles of the
# maps (PAPER_PROFILES of AK/s1/param-x/run_maps.py, the nu = 5 sweep), two significant digits
def sig2(x):
    return f"{x:.{max(0, 1 - math.floor(math.log10(abs(x))))}f}"


_sz, _l0 = set(), set()
for tag in CASES.values():
    _deck = open(os.path.join(P, f"{tag}.feb")).read()
    _sz.add(float(re.search(r"<z_spread>([^<]+)</z_spread>", _deck).group(1)))
    _z = [float(m.group(1)) for m in re.finditer(r'<node id="[\s\d]+">[^,]+,[^,]+,([^<]+)<', _deck)]
    _l0.add(round(max(_z) - min(_z), 9))
assert len(_sz) == 1 and len(_l0) == 1, (_sz, _l0)
SZ_LES, L0 = _sz.pop(), _l0.pop()
put("szRatioLesion", sig2(SZ_LES / L0))

# mesh and segment-length check of the Figure 4 lesions (AK/s1/mesh-study2/summary.csv, written by summary.py;
# point 16 of the simulated review, paper-review session, 4 Oct 2026). Meshes: gate (the original mesh rebuilt),
# z2 (axial elements halved), L30 (segment length 2 l_0). Changes in % relative to gate, per lesion.
_ms = {(r["case"], r["mesh"]): r for r in csv.DictReader(open(os.path.join(S1, "mesh-study2", "summary.csv")))}
_cases = sorted({c for c, m in _ms})
assert _cases == ["fus-el", "fus-elcc", "sac-el", "sac-elcc"], _cases
def _rel(col, mesh):
    return [100 * (float(_ms[(c, mesh)][col]) / float(_ms[(c, "gate")][col]) - 1) for c in _cases]
put("meshAxJmax", f"{max(abs(v) for v in _rel('Jmax', 'z2')):.1f}")
_x = max(abs(v) for v in _rel("apex_ratio", "z2"))          # bound, rounded UP to one significant digit ("less than")
_e = math.floor(math.log10(_x)); _m = math.ceil(round(_x / 10 ** _e, 9))
if _m == 10:
    _m, _e = 1, _e + 1
put("meshAxLamr", f"{_m * 10 ** _e:.{max(0, -_e)}f}")
for col, name in (("Jmax", "segJmax"), ("apex_ratio", "segLamr"), ("Jmin", "segJmin")):
    v = _rel(col, "L30")
    put(f"{name}Lo", f"{min(v):.0f}")
    put(f"{name}Hi", f"{max(v):.0f}")
put("szRatioLong", sig2(SZ_LES / (2 * L0)))
import importlib.util
_spec = importlib.util.spec_from_file_location("run_maps", os.path.join(S1, "param-x", "run_maps.py"))
_rm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_rm)
_szs = [z for nu, z in _rm.PAPER_PROFILES if (nu, z) in [(5.0, q) for q in _rm.SZ]]     # the s_z sweep, 1.0 to 3.0 mm
assert len(_szs) == len(_rm.SZ)
put("szRatioMin", sig2(min(_szs) / L0))
put("szRatioMax", sig2(max(_szs) / L0))

# ============================================================== parametric sweep
# the 11 profiles shown in the paper (Online Resource 1, maps; AK/s1/param-x, Stage I of Section 5.1, author
# 4 Oct 2026): s_z = 1.0, 1.4, 1.8, 2.2, 2.6, 3.0 mm at nu = 5 and nu = 1..6 at s_z = 3 mm, (5, 3) shared
# (PAPER_PROFILES of run_maps.py; the computed s_z = 2.8 map is not shown).
# The s1 MASTER csv also holds the same s_z sweep at nu = 3 (kept as data, not shown): it is excluded.
# BASE = the profile of latorre2020numerical and of the lesions of Figures 4-8 (nu = 5, s_z = 3 mm).
PAPER = set(_rm.PAPER_PROFILES)
BASE = (5.0, 3.0)
assert BASE in PAPER and len(PAPER) == 11
rows = [r for r in csv.DictReader(open(os.path.join(P, "MASTER_PARAMETRIC_STUDY.csv")))
        if r["alph"] == "ON" and "unit-frac" not in r["folder"]
        and (float(r["z3_nu"]), round(float(r["z4_zod"]), 1)) in PAPER]
GRID = 400
LAM_MAX = 2.5                  # admissible: converged and lambda_r/lambda_ro <= 2.5 (apex)


put("parLamMax", f"{LAM_MAX:.1f}")


def conf(typ, nu, sz):
    return [r for r in rows if r["type"] == typ and float(r["z3_nu"]) == nu
            and abs(float(r["z4_zod"]) - sz) < 1e-9]


def admissible(r):
    return float(r["z1_lamb_r"]) <= LAM_MAX


def adm(g):
    return [r for r in g if admissible(r)]


def axes_limit(g):
    a = adm(g)
    de = max((float(r["x_elastin_deg"]) for r in a if float(r["y_collagen_deg"]) == 0.0), default=0)
    dc = max((float(r["y_collagen_deg"]) for r in a if float(r["x_elastin_deg"]) == 0.0), default=0)
    return de, dc


de1, dc1 = axes_limit(conf("axis", 5.0, 1.0))           # s_z sweep at nu = 5: shortest and longest extent
de28, dc28 = axes_limit(conf("axis", 5.0, max(_rm.SZ)))   # longest s_z shown (3.0 mm since 4 Oct 2026)
put("parFusDeOne", f"{100 * de1:.0f}")
_, dcs = axes_limit(conf("asym", 5.0, 1.0))
put("parSacDcOne", f"{100 * dcs:.0f}")
deb, dcb = axes_limit(conf("axis", *BASE))               # base profile
put("parFusDeBase", f"{100 * deb:.0f}")
put("parFusDcBase", f"{100 * dcb:.0f}")
_, dcsb = axes_limit(conf("asym", *BASE))
put("parSacDcBase", f"{100 * dcsb:.0f}")
put("parFusLamrTrellis", f"{max(float(r['z1_lamb_r']) for z in _rm.SZ for r in adm(conf('axis', 5.0, z))):.1f}")
put("parFusDcOne", f"{100 * dc1:.0f}")
put("parFusDeLong", f"{100 * de28:.0f}")
put("parFusDcLong", f"{100 * dc28:.0f}")
for typ, name in (("axis", "Fus"), ("asym", "Sac")):
    g = adm(conf(typ, *BASE))
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
    put(f"par{name}LamrSweepNu", f"{float(best[1]['z3_nu']):.0f}")
    put(f"par{name}LamrSweepSz", f"{float(best[1]['z4_zod']):.1f}")
ratios = []
for typ in ("axis", "asym"):
    ga = adm([r for r in rows if r["type"] == typ])
    A = np.array([[1.0, float(r["x_elastin_deg"]), float(r["y_collagen_deg"])] for r in ga])
    c = np.linalg.lstsq(A, np.array([float(r["z2_J"]) for r in ga]), rcond=None)[0]
    ratios.append(c[2] / c[1])
put("parSensRatio", f"{max(ratios):.1f}")
# the same fit for each paper profile separately: range over the 12 profiles of (dJ/dd_c)/(dJ/dd_e)
# (point 9 of the simulated review, paper-review session, 4 Oct 2026)
for typ, name in (("axis", "Fus"), ("asym", "Sac")):
    rp = []
    for prof in sorted(PAPER):
        ga = adm(conf(typ, *prof))
        assert len(ga) >= 10, (typ, prof, len(ga))
        A = np.array([[1.0, float(r["x_elastin_deg"]), float(r["y_collagen_deg"])] for r in ga])
        c = np.linalg.lstsq(A, np.array([float(r["z2_J"]) for r in ga]), rcond=None)[0]
        rp.append(c[2] / c[1])
    put(f"parSensRatio{name}Min", f"{min(rp):.2f}")
    put(f"parSensRatio{name}Max", f"{max(rp):.2f}")

# beyond the admissible region: the csv holds only the converged grid points, so a missing
# point is a run that did not converge. At the base profile (nu = 5, s_z = 3 mm): share of the grid converged
# with lambda_r/lambda_ro > LAM_MAX, share not converged, largest converged J. Over the 11 paper
# profiles, the non-converged points that have an admissible neighbour on the grid were
# recomputed with a tenth of the step (param-x/edge-retry/results.json): each one must either
# converge with lambda_r/lambda_ro > LAM_MAX or pass LAM_MAX before Newton fails, which the text states.
for typ, name in (("axis", "Fus"), ("asym", "Sac")):
    g = conf(typ, *BASE)
    put(f"par{name}Beyond", f"{100 * sum(not admissible(r) for r in g) / GRID:.0f}")
    put(f"par{name}NoConv", f"{100 * (GRID - len(g)) / GRID:.0f}")
    put(f"par{name}JConvMax", f"{max(float(r['z2_J']) for r in g):.0f}")
_step = 1.0 / 19.0
RETRY = [q for q in json.load(open(os.path.join(P, "edge-retry", "results.json")))
         if (q["nu"], q["sz"]) in PAPER]      # the s_z = 2.8 retries are outside the paper profiles since 4 Oct 2026
edge_retry = []
_noconv = 0                                         # grid points without a converged solution, all profiles
for (nu, sz) in PAPER:
    for typ in ("axis", "asym"):
        pts = {(round(float(r["x_elastin_deg"]) / _step), round(float(r["y_collagen_deg"]) / _step)):
               float(r["z1_lamb_r"]) for r in conf(typ, nu, sz)}
        _noconv += 400 - len(pts)
        for i in range(20):
            for j in range(20):
                if (i, j) in pts:
                    continue
                nb = [pts.get(q) for q in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1))]
                if any(v is not None and v <= LAM_MAX for v in nb):
                    hit = [q for q in RETRY if (q["type"], q["nu"], q["sz"], q["i"], q["j"]) == (typ, nu, sz, i, j)]
                    if not hit or hit[0]["lam_last"] <= LAM_MAX:
                        raise SystemExit(f"non-converged point next to an admissible one, not recomputed: {typ} nu={nu} s_z={sz} {i},{j}")
                    edge_retry.append(hit[0])
assert len(edge_retry) == len(RETRY), (len(edge_retry), len(RETRY))
WORDS = "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen seventeen eighteen nineteen twenty".split()
put("parEdgeN", WORDS[len(edge_retry)])
put("parNoConvN", f"{_noconv}")
put("parProfilesN", f"{len(PAPER)}")
_conv = [q for q in edge_retry if q["converged"]]
_fail = [q for q in edge_retry if not q["converged"]]
put("parEdgeConvN", WORDS[len(_conv)])
put("parEdgeFailN", WORDS[len(_fail)])
assert len(_conv) == 1                     # the text gives the ratio of the one converged point
put("parEdgeConvLam", f"{_conv[0]['lam_last']:.1f}")
# smallest ratio reached before Newton fails, rounded DOWN to one decimal (the text says "exceeds")
put("parEdgeFailLam", f"{math.floor(10 * min(q['lam_last'] for q in _fail)) / 10:.1f}")
# relative volume of the end ring of the segment (element 1, z = 0 to 0.9 mm, axial displacement fixed) at the
# last converged increment of the failed retries; the ring is uniform around the circumference and the two ends
# behave alike (xplt files read on 5 Oct 2026). Every failed retry reports negative Jacobians in its log.
_jend = []
for q in _fail:
    _d = os.path.join(P, "edge-retry", f"{q['type']}_nu{q['nu']:g}_sz{q['sz']:g}_{q['i']}_{q['j']}_dt0.1_b")
    _last = open(os.path.join(_d, f"elem-{q['type']}-B_MBE.txt")).read().split("*Step")[-1]
    _row = [l for l in _last.splitlines() if l and l[0].isdigit()][0].split(",")
    assert int(_row[0]) == 1
    _jend.append(float(_row[4]))
    assert "negative jacobians" in open(os.path.join(_d, "_run.log")).read().lower(), _d
assert max(_jend) < 0.01, _jend
put("parEdgeJEnd", f"{max(_jend):.2g}")
# relative volume at the edge of convergence: J at the converged grid points that have a
# non-converged neighbour on the grid, over the 11 paper profiles (range and median), and
# the largest converged J of the whole sweep
for typ, name in (("axis", "Fus"), ("asym", "Sac")):
    edge, edge_lr, jall = [], [], []
    for (nu, sz) in PAPER:
        pts = {(round(float(r["x_elastin_deg"]) / _step), round(float(r["y_collagen_deg"]) / _step)):
               (float(r["z2_J"]), float(r["z1_lamb_r"])) for r in conf(typ, nu, sz)}
        jall += [v[0] for v in pts.values()]
        on_edge = [v for (i, j), v in pts.items()
                   if any(q not in pts and 0 <= q[0] < 20 and 0 <= q[1] < 20
                          for q in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)))]
        edge += [v[0] for v in on_edge]
        edge_lr += [v[1] for v in on_edge]
    put(f"par{name}JEdgeMin", f"{min(edge):.0f}" if min(edge) >= 5 else f"{min(edge):.1f}")
    put(f"par{name}JEdgeMed", f"{np.median(edge):.0f}")
    put(f"par{name}JEdgeMax", f"{max(edge):.0f}")
    put(f"par{name}JConvMaxAll", f"{max(jall):.0f}")

# sensitivity to the transition decay nu at s_z = 3 mm
xs = sorted({float(r["x_elastin_deg"]) for r in rows})
ys = sorted({float(r["y_collagen_deg"]) for r in rows})
DE = min(xs, key=lambda v: abs(v - 0.20))
DC = min(ys, key=lambda v: abs(v - 0.10))
put("nuDe", f"{100 * DE:.0f}\\%")   # in percent, as the maps and Section 5.5
put("nuDc", f"{100 * DC:.0f}\\%")
for typ, name in (("axis", "Fus"), ("asym", "Sac")):
    Jn, fr = [], []
    for nu in range(1, 7):
        g = conf(typ, float(nu), 3.0)
        hit = [r for r in g if abs(float(r["x_elastin_deg"]) - DE) < 1e-9
               and abs(float(r["y_collagen_deg"]) - DC) < 1e-9]
        Jn.append(float(hit[0]["z2_J"]))
        fr.append(100.0 * len(adm(g)) / GRID)
    put(f"nu{name}JOne", f"{Jn[0]:.2f}")
    put(f"nu{name}JSix", f"{Jn[-1]:.2f}")
    put(f"nu{name}Spread", f"{100 * (Jn[-1] - Jn[0]) / Jn[0]:.0f}")
    put(f"nu{name}IncTail", f"{max(abs(np.diff(Jn[3:]))):.2f}")
    # largest J over the admissible states of each nu map; the text gives the range over nu
    jadm = [max(float(r["z2_J"]) for r in adm(conf(typ, float(nu), 3.0))) for nu in range(1, 7)]
    put(f"nu{name}JAdmLo", f"{min(jadm):.1f}")
    put(f"nu{name}JAdmHi", f"{max(jadm):.1f}")
    put(f"nu{name}AdmOne", f"{fr[0]:.0f}")
    put(f"nu{name}AdmSix", f"{fr[-1]:.0f}")
    if typ == "axis":
        put("nuFusJThree", f"{Jn[2]:.3f}")
        put("nuFusJFive", f"{Jn[4]:.3f}")
        put("nuFusThreeFive", pct((Jn[4] - Jn[2]) / Jn[2], 1))

# ======================================================================== severity ladder
# AK/severity-fig/ladder-numbers.json (make_ladder_numbers.py, from ladder.csv): lesions with
# d_e = d_c up to the largest admissible value, and the severities of latorre2020numerical
import json

with open(os.path.join(AK, "severity-fig", "ladder-numbers.json")) as f:
    LAD = json.load(f)
for key, name in (("fusiform", "Fus"), ("saccular", "Sac")):
    L = LAD[key]
    put(f"lad{name}D", f"{L['d_max']:.0f}")
    put(f"lad{name}J", f"{L['J_apex_dmax']:.2f}")
    put(f"lad{name}Lamr", f"{L['lamr_dmax']:.2f}")
    put(f"lad{name}DBelow", f"{L['d_below']:.0f}")
    put(f"lad{name}JBelow", f"{L['J_apex_below']:.2f}")
    put(f"lad{name}DAbove", f"{L['d_above']:.0f}")
    put(f"lad{name}JAbove", f"{L['J_apex_above']:.2f}")
    assert L["converged_above"]
    lat = {(c["d_e"], c["d_c"]): c for c in L["latorre"]}
    for tag, k in (("MildE", (32, 0)), ("Mild", (32, 9)), ("SevE", (62, 0)), ("Sev", (62, 23))):
        put(f"lat{name}{tag}J", f"{lat[k]['J_apex']:.2f}" if lat[k]["converged"] else "none")
# the statements of Section 6.2.1 on the severities of latorre2020numerical
latF = {(c["d_e"], c["d_c"]): c for c in LAD["fusiform"]["latorre"]}
latS = {(c["d_e"], c["d_c"]): c for c in LAD["saccular"]["latorre"]}
assert all(c["converged"] for c in latS.values())
assert latF[(32, 0)]["converged"] and latF[(32, 9)]["converged"]
assert not latF[(62, 0)]["converged"] and not latF[(62, 23)]["converged"]

# ======================================================================== Figure 6 (severity ladder, d_e = d_c)
# rABC-6 (2026-10-03, proposal of review/rabc6_proposals_fable.md): the eight decks drawn by
# AK/severity-fig/render_ladder_matrix.py, read as the Figure 4 fields: J over the elements (*-ALL.txt),
# apex radius ratio at node 20 (*-NODES.txt); checked against ladder.csv. The lad* macros above (criterion
# J <= 4, ladder-numbers.json) are not used in the text.
SF = ak("severity-fig")
LADDER = {"Fus": (("A", "lad-fus-d13p000-lev2", 13), ("B", "lad-fus-d17p000-lev", 17),
                  ("C", "lad-fus-d19p000-v5", 19), ("D", "lad-fus-d20p000-v5", 20)),
          "Sac": (("E", "lad-sac-d22p000-lev2", 22), ("F", "lad-sac-d27p000-lev2", 27),
                  ("G", "lad-sac-d32p000-lev2", 32), ("H", "lad-sac-d37p000-lev2", 37))}
with open(os.path.join(SF, "ladder.csv")) as f:
    _lcsv = {r["deck"]: r for r in csv.DictReader(f)}
for name, cases in LADDER.items():
    for tag, deck, d in cases:
        r = _lcsv[deck + ".feb"]
        assert float(r["d_e_percent"]) == d == float(r["d_c_percent"]) and r["converged"] == "yes"
        J = [v[0] for v in blocks(os.path.join(SF, deck + "-ALL.txt"))[-1][1].values()]
        lamr = node_ratio(blocks(os.path.join(SF, deck + "-NODES.txt")), 20)
        assert abs(max(J) - float(r["J_max_elements"])) < 5e-5
        assert abs(lamr - float(r["lambda_r_over_lambda_ro_node20"])) < 5e-5
        put(f"lad{tag}Deg", f"{d:.0f}")            # d_e = d_c of the panel, %
        put(f"lad{tag}Jmax", f"{max(J):.2f}")      # largest J over the elements (panel label)
        put(f"lad{tag}Jmin", f"{min(J):.2f}")      # smallest J over the elements (contracted ends)
        put(f"lad{tag}Lamr", f"{lamr:.2f}")        # apex radius ratio, node 20 (as fieldFusElLamr)
# the next computed equal losses: lambda_r/lambda_ro above LAM_MAX (ladder.csv)
for name, deck in (("Fus", "lad-fus-d21p000-bis.feb"), ("Sac", "lad-sac-d38p000-v5.feb")):
    r = _lcsv[deck]
    assert r["converged"] == "yes" and float(r["lambda_r_over_lambda_ro_node20"]) > LAM_MAX
    put(f"lad{name}DNext", f"{float(r['d_e_percent']):.0f}")
assert latF[(32, 9)]["J_apex"] > 4 > latF[(32, 0)]["J_apex"]
assert max(c["J_apex"] for c in latS.values()) < 4

# ======================================================================== Figures 7 and 8 (profiles at fixed loss)
# 2026-10-04 (author): the losses of Figure 6(d) and (h), d_e = d_c = 20 % (fusiform) and 37 % (saccular), with
# s_z = 1.2, 1.8, 2.4, 3.0 mm at nu = 5 (Figure 7) and nu = 1, 3, 5, 7 at s_z = 3 mm (Figure 8); decks in
# AK/profile-fig (make_decks.py), the nu = 5, s_z = 3 mm panels are the Figure 6 runs. Apex radius ratio at node 20,
# rounded from the full value (2.3053 -> 2.31).
PF = ak("profile-fig")
_PREF = {"Fus": os.path.join(SF, "lad-fus-d20p000-v5"), "Sac": os.path.join(SF, "lad-sac-d37p000-lev2")}
_PSER = {"Sz": [(5, z) for z in (1.2, 1.8, 2.4, 3.0)], "Nu": [(n, 3.0) for n in (1, 3, 5, 7)]}
for name in ("Fus", "Sac"):
    for ser, prof in _PSER.items():
        lam = []
        for nu, sz in prof:
            st = _PREF[name] if (nu, sz) == (5, 3.0) else os.path.join(
                PF, f"prof-{name.lower()}-nu{nu}-sz{sz:.1f}".replace(".", "p"))
            feb = open(st + ".feb").read()
            assert f"<z_spread>{sz:g}</z_spread>" in feb or f"<z_spread>{sz:.1f}</z_spread>" in feb, st
            assert f"<z_nu>{nu}</z_nu>" in feb and f"<t_nu>{nu}</t_nu>" in feb, st
            assert "N O R M A L   T E R M I N A T I O N" in open(st + ".log").read(), st
            nb = blocks(st + "-NODES.txt")
            assert abs(nb[-1][0] - 11.0) < 1e-9, st
            lam.append(node_ratio(nb, 20))
        assert all(b > a for a, b in zip(lam, lam[1:])), (name, ser, lam)   # "rises" in the text
        put(f"prof{ser}{name}Lo", f"{lam[0]:.2f}")    # apex radius ratio, first column
        put(f"prof{ser}{name}Hi", f"{lam[-1]:.2f}")   # apex radius ratio, last column
        # text: only the saccular lesion with nu = 7 exceeds the reporting limit
        assert all(v <= LAM_MAX for v in lam) == (not (name == "Sac" and ser == "Nu")), (name, ser, lam)
        if name == "Sac" and ser == "Nu":
            assert lam[-1] > LAM_MAX and all(v <= LAM_MAX for v in lam[:-1])

# ============================================================ rate-dependent lesion (Section 5.3)
# AK/s1/fullkg-sign (author 5 Oct 2026): full-KG of the fusiform lesion of Figure 4(a), s_load = 100 days, from the MBE deck
# param-x/matrix-fus-el.feb with MBE = 0 (setup_run.py). a: deposition stretches of Table 2; b: unit elastin deposition
# stretches. Element 10 is at the apex (z = 7.2 mm), node 20 likewise; -ALL.txt holds J of every element.
def _fkg_blocks(path):
    out = []
    for blk in open(path).read().split("*Step")[1:]:
        t = float(re.search(r"\*Time\s*=\s*(\S+)", blk).group(1))
        out.append((t, np.array([[float(v) for v in l.split(",")] for l in blk.splitlines() if l and l[0].isdigit()])))
    return out
_FK = os.path.join(S1, "fullkg-sign")
_la = open(os.path.join(_FK, "a", "fus-el-full.log")).read()
_lb = open(os.path.join(_FK, "b", "fus-el-full.log")).read()
assert "E R R O R   T E R M I N A T I O N" in _la and "N O R M A L   T E R M I N A T I O N" in _lb
_ea = _fkg_blocks(os.path.join(_FK, "a", "matrix-fus-el-eA.txt"))
_eb = _fkg_blocks(os.path.join(_FK, "b", "matrix-fus-el-eA.txt"))
_ja = [b[0, 4] for t, b in _ea if t > 1 + 1e-9]
assert all(y < x for x, y in zip(_ja, _ja[1:])) and _ja[0] < 1.0      # the apex volume decreases from the first increment
_feb = open(os.path.join(_FK, "a", "fus-el-full.feb")).read()
assert "<gr_end>101.0</gr_end>" in _feb and "<MBE>0</MBE>" in _feb
_zc = float(re.search(r"<z_center>([^<]+)</z_center>", _feb).group(1))
_nod = {int(m.group(1)): tuple(map(float, m.group(2).split(","))) for m in re.finditer(r'<node id="\s*(\d+)">([^<]+)</node>', _feb)}
_els = {int(m.group(1)): [int(v) for v in m.group(2).split(",")] for m in re.finditer(r'<elem id="\s*(\d+)">([^<]+)</elem>', _feb)}
_all = _fkg_blocks(os.path.join(_FK, "a", "matrix-fus-el-ALL.txt"))
assert abs(_all[-1][0] - _ea[-1][0]) < 1e-9
_row = _all[-1][1][_all[-1][1][:, 1].argmax()]
_zmax = np.mean([_nod[n][2] for n in _els[int(_row[0])][:8]])
put("fkgFailS", f"{_ea[-1][0]:.0f}")
put("fkgFailPct", f"{_ea[-1][0] - 1:.0f}")                            # degradation period: s = 1 to 101 days
put("fkgApexJ", f"{_ja[-1]:.2f}")
put("fkgOffApexJ", f"{_row[1]:.1f}")
put("fkgOffApexDz", f"{abs(_zmax - _zc):.1f}")
put("fkgUnitJ", f"{_eb[-1][1][0, 4]:.3f}")

# ======================================================================== output
_bad = sorted(f for f in OUTSIDE_S1 if not f.startswith(S1_KEEP))
assert not _bad, f"read outside s1/: {_bad}"
with _open(os.path.join(HERE, "numbers.tex" if LO1 else "numbers.lo0.tex"), "w") as f:
    f.write("% generated by make_numbers.py; do not edit\n")
    for k in sorted(MACROS):
        f.write(f"\\newcommand{{\\{k}}}{{{MACROS[k]}}}\n")
with _open(os.path.join(HERE, "numbers-inputs.txt" if LO1 else "numbers-inputs.lo0.txt"), "w") as f:
    f.write("\n".join(sorted(INPUTS)) + "\n")
print(f"numbers.tex: {len(MACROS)} macros from {len(INPUTS)} files")
for k in sorted(MACROS):
    print(f"  {k:24s} {MACROS[k]}")
