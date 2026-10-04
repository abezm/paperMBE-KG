#!/usr/bin/env python3
"""One-thread wall time of the fusiform s_load = 100 full-KG run at ds = 1, 2, 5 days, and of
the MBE deck with the same one-step pressurization as full-KG (mbe-1p10) and the paper's MBE deck (mbe-20), under the same conditions (same conventions as ../timing/timing.py:
OMP_NUM_THREADS = MKL_NUM_THREADS = 1, sequential, minimum of 3 repetitions).

Every repetition runs its own copy of the deck (T-<case>-r<k>.feb) with unique logger names,
and refuses to start if any of its outputs exists. Writes timing_fair2.csv (second batch: full-KG ds = 1 over the full 130 days of growth, the hyperelastic pressurization in 10 and in 1 increments, and mbe-1p10 again as the session reference).

alp1/ (2 Oct 2026): copy of the script of the same name in ../../speedup-fair/; executable is the copy
~/FEBio-stage1, decks carry <alp_update>1</alp_update>; only the cases the paper reads are kept.
alp2/ (3 Oct 2026): severities fus cc_in 0.95; the full-KG Delta s = 1 day run is timed in stop/.
"""
import csv
import os
import re
import subprocess
import time
from pathlib import Path

FEBIO = __import__("os").environ.get("FEBIO", "febio4")   # alp1/: copy with F Fo^-1 and tangent terms A, B
HERE = Path(__file__).parent
ENV = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
SRC = [("hyp-1", HERE / "hyp-1.feb"),
       ("mbe-1p10-b", HERE / "mbe-1p10.feb")]
REPS = 3

rows = []
for case, src in SRC:
    txt = src.read_text(encoding="ISO-8859-1")
    walls, stats = [], None
    for k in range(1, REPS + 1):
        name = f"T-{case}-r{k}"
        t = re.sub(r'file="([^"]+)"', lambda m: f'file="{name}_{m.group(1)}"', txt)
        deck = HERE / f"{name}.feb"
        outs = [deck.with_suffix(".log"), deck.with_suffix(".xplt")] + \
               [HERE / f for f in re.findall(r'file="([^"]+)"', t)]
        clash = [str(o) for o in outs if o.exists()] + ([str(deck)] if deck.exists() else [])
        if clash:
            raise SystemExit(f"refusing {name}: exists {clash}")
        deck.write_text(t, encoding="ISO-8859-1")
        t0 = time.perf_counter()
        subprocess.run([FEBIO, "-i", deck.name, "-silent"], cwd=HERE, env=ENV,
                       capture_output=True, text=True)
        walls.append(time.perf_counter() - t0)
        log = deck.with_suffix(".log").read_text()
        if "N O R M A L" not in log:
            raise SystemExit(f"{name}: did not terminate normally")
        # FEBio prints one block per stage, then (multi-stage only) a SUMMARY block with totals:
        # [-1] is therefore the total over all stages; the per-stage blocks precede it.
        ns = [int(x) for x in re.findall(r"Number of time steps completed[\. ]*: (\d+)", log)]
        ni = [int(x) for x in re.findall(r"Total number of equilibrium iterations[\. ]*: (\d+)", log)]
        per = list(zip(ns[:-1], ni[:-1])) if len(ns) > 1 else [(ns[0], ni[0])]
        assert sum(a for a, _ in per) == ns[-1] and sum(b for _, b in per) == ni[-1], (ns, ni)
        stats = (ns[-1], ni[-1], " | ".join(f"{a}/{b}" for a, b in per))
        print(case, k, f"{walls[-1]:.2f} s", stats, flush=True)
    rows.append(dict(case=case, deck=src.name, reps=REPS, wall_min_s=f"{min(walls):.2f}",
                     wall_all=" ".join(f"{w:.2f}" for w in walls),
                     steps=stats[0], newton=stats[1], per_stage_steps_iters=stats[2]))
with open(HERE / "timing_fair2.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
for r in rows:
    print(r)
