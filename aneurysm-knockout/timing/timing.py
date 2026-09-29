#!/usr/bin/env python3
"""Wall time, step count and Newton-iteration count of the six convergence-table runs.

The decks are copies of those that produced the convergence tables
(conv-tables/conv-*.feb, uniax-fractions/ds-study/u_ds0.5.feb, full-KG/ds-study/a_ds1.0.feb),
run here so that no earlier log is overwritten. FEBio runs on one thread
(OMP_NUM_THREADS = MKL_NUM_THREADS = 1). Each uniaxial run is repeated 15 times and the
median is kept; each aneurysm run 3 times and the minimum is kept.

Writes timing.csv: case, reps, wall_s (median or min), wall_all, steps, newton, linsolve_s
"""
import csv
import os
import re
import statistics
import subprocess
import time
from pathlib import Path

FEBIO = __import__("os").environ.get("FEBIO", "febio4")
HERE = Path(__file__).parent
ENV = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
RUNS = [  # case, deck, reps, statistic
    ("uniax-MBE", "conv-uniax-MBE.feb", 15, "median"),
    ("uniax-full", "u_ds0.5.feb", 15, "median"),
    ("uniax-HYP", "conv-uniax-HYP.feb", 15, "median"),
    ("aneu-MBE", "conv-aneu-MBE.feb", 3, "min"),
    ("aneu-full", "a_ds1.0.feb", 3, "min"),
    ("aneu-HYP", "conv-aneu-HYP.feb", 3, "min"),
]


def hms(s):
    h, m, x = s.split(":")
    return 3600 * int(h) + 60 * int(m) + float(x)


rows = []
for case, deck, reps, stat in RUNS:
    log = HERE / f"{Path(deck).stem}.log"
    walls = []
    for _ in range(reps):
        t0 = time.perf_counter()
        subprocess.run([FEBIO, "-i", deck, "-o", log.name], cwd=HERE, env=ENV,
                       capture_output=True, text=True)
        walls.append(time.perf_counter() - t0)
    txt = log.read_text()
    if "N O R M A L" not in txt:
        raise SystemExit(f"{deck}: did not terminate normally")
    steps = int(re.search(r"Number of time steps completed[\. ]*: (\d+)", txt).group(1))
    newton = int(re.search(r"Total number of equilibrium iterations[\. ]*: (\d+)", txt).group(1))
    lin = re.search(r"Time in linear solver: ([\d:.]+)", txt)
    wall = statistics.median(walls) if stat == "median" else min(walls)
    rows.append(dict(case=case, reps=reps, stat=stat, wall_s=f"{wall:.4f}",
                     wall_all=" ".join(f"{w:.4f}" for w in walls), steps=steps,
                     newton=newton, linsolve_s=hms(lin.group(1)) if lin else ""))
    print(rows[-1], flush=True)

with open(HERE / "timing.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print("wrote", HERE / "timing.csv")
