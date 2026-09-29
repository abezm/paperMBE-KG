#!/usr/bin/env python3
"""k_J (volumetric penalty) sensitivity study -- referee item M5.

Stage I (t <= gr_start) uses a penalized volumetric response
    S = Se + C^-1 * k_J * ln(Jdep*J),   Jdep = 0.9999   (hard-coded in the .cpp)
and sigma_vo -- the reference scalar of the whole growth law -- is captured
from THAT stress.  The theory, however, assumes det(Fe) = 1.  So sigma_vo,
and hence every downstream result, formally depends on k_J.  All the paper's
decks use k_J = 1e6 but never report it.

This sweeps k_J over four orders of magnitude on
  (a) the fusiform aneurysm (matrix-fus-el.feb, 11 steps, gr_start = 1)
  (b) the uniaxial coupon    (uniax-MBE.feb,     50 steps, gr_start = 1)
and reports, for each:
    J at end of Stage I      -> how well incompressibility is actually enforced
    sigma_vo = (sx+sy+sz)/3  -> the reference scalar
    final J                  -> the actual result quantity (growth/dilatation)
"""
import re, shutil, subprocess, sys, time
from pathlib import Path

FEBIO = __import__("os").environ.get("FEBIO", "febio4")
HERE = Path(__file__).parent
SRC = {
    "aneu": HERE.parent / "param-x" / "matrix-fus-el.feb",
    "uniax": HERE.parent / "uniax-fractions" / "uniax-MBE.feb",
}
KJ = ["1e3", "1e4", "1e5", "1e6", "1e7", "1e8"]   # 1e6 = the paper's value


def build(case, kj):
    src = SRC[case]
    txt = src.read_text()
    tag = f"{case}_kJ{kj}"
    # 1) the penalty itself
    txt, n = re.subn(r"<k_J>[^<]*</k_J>", f"<k_J>{kj}</k_J>", txt)
    assert n == 1, f"expected exactly one <k_J> tag, found {n}"
    # 2) redirect every logfile output so runs do not overwrite each other
    txt = re.sub(r'file="([^"]+)\.txt"', rf'file="{tag}_\1.txt"', txt)
    feb = HERE / f"{tag}.feb"
    feb.write_text(txt)
    return feb, tag


def gr_start(feb):
    m = re.search(r"<gr_start>([^<]+)</gr_start>", feb.read_text())
    return float(m.group(1))


def read_elem(path):
    """FEBio logfile format: '*Time = <t>' blocks, then 'id,v1,v2,...' rows."""
    out, t = {}, None
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith("*Time"):
            t = float(line.split("=")[1])
        elif line and line[0].isdigit() and t is not None:
            f = line.split(",")
            out.setdefault(t, {})[int(f[0])] = [float(x) for x in f[1:]]
    return out


def main():
    rows = []
    for case in ("aneu", "uniax"):
        gs = gr_start(SRC[case])
        for kj in KJ:
            feb, tag = build(case, kj)
            t0 = time.time()
            r = subprocess.run([FEBIO, "-i", feb.name, "-o", f"{tag}.log"],
                               cwd=HERE, capture_output=True, text=True)
            wall = time.time() - t0
            log = (HERE / f"{tag}.log").read_text()
            ok = "N O R M A L" in log
            if not ok:
                rows.append((case, kj, None, None, None, None, "DIVERGED"))
                print(f"{case:5s} k_J={kj:4s}  *** DID NOT CONVERGE ***", flush=True)
                continue
            # first element_data file for this case
            cand = sorted(HERE.glob(f"{tag}_elem1_*.txt")) + \
                   sorted(HERE.glob(f"{tag}_*-eB.txt"))
            d = read_elem(cand[0])
            eid = sorted(d[min(d)].keys())[0]
            ts = sorted(d)
            # Stage-I end == the step at gr_start
            t_s1 = min(ts, key=lambda t: abs(t - gs))
            sx, sy, sz, J1 = d[t_s1][eid][:4]
            svo = (sx + sy + sz) / 3.0
            Jf = d[ts[-1]][eid][3]
            rows.append((case, kj, J1, svo, Jf, wall, "ok"))
            print(f"{case:5s} k_J={kj:4s}  J(stageI)={J1:.8f}  "
                  f"sigma_vo={svo:9.4f} kPa  J(final)={Jf:.6f}  [{wall:.0f}s]",
                  flush=True)

    print("\n=== summary (relative to the paper's k_J = 1e6) ===")
    for case in ("aneu", "uniax"):
        base = [r for r in rows if r[0] == case and r[1] == "1e6" and r[6] == "ok"]
        if not base:
            continue
        _, _, J1b, svob, Jfb, _, _ = base[0]
        print(f"\n{case}:  baseline sigma_vo = {svob:.4f} kPa, J_final = {Jfb:.6f}")
        print(f"{'k_J':>6} {'J(stageI)':>12} {'sigma_vo':>11} {'d_svo %':>9} "
              f"{'J(final)':>10} {'d_Jf %':>8}")
        for r in rows:
            if r[0] != case:
                continue
            if r[6] != "ok":
                print(f"{r[1]:>6} {'DIVERGED':>12}")
                continue
            _, kj, J1, svo, Jf, _, _ = r
            print(f"{kj:>6} {J1:12.8f} {svo:11.4f} "
                  f"{100*(svo-svob)/abs(svob):9.3f} {Jf:10.6f} "
                  f"{100*(Jf-Jfb)/Jfb:8.3f}")


if __name__ == "__main__":
    main()
