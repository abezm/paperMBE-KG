#!/usr/bin/env python3
"""Timing of the uniaxial test with the equilibrated and hyperelastic runs ending at t = 11.

The load reaches 3 F_o at t = 11 and stays constant afterwards. The equilibrated (MBE) and the
hyperelastic solutions do not depend on time, so the 40 steps that ../timing/ ran after t = 11
repeated the same state. Here the two runs end at t = 11 and carry a must point at t = 1, so that
Stage I ends at F_o (in ../timing/ their first step was cut and Stage I ended at t = 0.8627).
The rate-dependent run (u_ds0.5.feb) is copied unchanged and timed again in the same session.

Protocol of ../run_all.py: one thread, median of 15 runs. Executable unmodified.
usage (cwd = timing-t11):  run_t11.py

alp1/ (2 Oct 2026): copy of ../../uniax-lo1/timing-t11/run_t11.py; executable is the copy ~/FEBio-stage1;
the input files are already in this folder (t = 11, must point, <alp_update>1</alp_update>), read unchanged.
timing-stop/ (2 Oct 2026): the full-KG run (u_ds0.5) ends at the stopping time s = 81.5 days (find/; J changes
by less than 1e-4 (J - 1) over the last s_load/20 = 0.5 day, after the end of the load), 165 steps.
"""
import csv, os, re, statistics, subprocess, time
from pathlib import Path

FEBIO = __import__("os").environ.get("FEBIO", "febio4")
HERE = Path(__file__).resolve().parent
SRC = HERE
ENV1 = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
REPS = 15
MUST = ('\t\t<load_controller id="2" name="MustPoints" type="loadcurve">\n\t\t\t<interpolate>LINEAR</interpolate>\n'
        '\t\t\t<extend>CONSTANT</extend>\n\t\t\t<points>\n\t\t\t\t<pt>0,1.0</pt>\n\t\t\t\t<pt>1.0,1.0</pt>\n'
        '\t\t\t\t<pt>10000,1.0</pt>\n\t\t\t</points>\n\t\t</load_controller>\n\t</LoadData>')


def rep(s, a, b, name):
    assert s.count(a) == 1, (name, a, s.count(a))
    return s.replace(a, b)


def t11(text, name):
    s = rep(text, "<time_steps>50</time_steps>", "<time_steps>11</time_steps>", name)
    s = rep(s, "<dtmax>1.0</dtmax>", '<dtmax lc="2">1.0</dtmax>', name)
    if 'load_controller id="2"' in s:                      # unused curve of the MBE deck: its values become 1.0
        for t in ("0", "1.0", "10000"):
            s = rep(s, f"<pt>{t},0.25</pt>", f"<pt>{t},1.0</pt>", name)
    else:
        s = rep(s, "\t</LoadData>", MUST, name)
    return s


def summary(log):
    txt = Path(log).read_text(errors="replace")
    ok = "N O R M A L   T E R M I N A T I O N" in txt
    g = lambda pat: (re.search(pat, txt) or [None, None])[1]
    return ok, g(r"Number of time steps completed[\. ]*: (\d+)"), g(r"Total number of equilibrium iterations[\. ]*: (\d+)")


rows = []
for case, deck, change in (("uniax-MBE", "conv-uniax-MBE.feb", False), ("uniax-full", "u_ds0.5.feb", False),
                           ("uniax-HYP", "conv-uniax-HYP.feb", False)):
    s0 = (SRC / deck).read_text(encoding="ISO-8859-1")
    (HERE / deck).write_text(t11(s0, deck) if change else s0, encoding="ISO-8859-1")
    log = HERE / (Path(deck).stem + ".log")
    walls = []
    for _ in range(REPS):
        t0 = time.perf_counter()
        subprocess.run([FEBIO, "-i", deck, "-o", log.name], cwd=HERE, env=ENV1, capture_output=True, text=True)
        walls.append(time.perf_counter() - t0)
    ok, steps, newton = summary(log)
    if not ok:
        raise SystemExit(f"{deck}: did not terminate normally")
    rows.append(dict(case=case, reps=REPS, stat="median", wall_s=f"{statistics.median(walls):.4f}",
                     wall_all=" ".join(f"{w:.4f}" for w in walls), steps=steps, newton=newton))
    print(rows[-1], flush=True)
with open(HERE / "timing.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
