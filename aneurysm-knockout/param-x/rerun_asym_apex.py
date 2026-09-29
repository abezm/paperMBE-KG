"""Re-run the saccular (asym) configurations of the paper with the probe at the apex.

TAA-asym-MBE-BASE.feb recorded element 90 and node 1209, 72 deg from the centre of the
insult; it now records element 10 and node 20, the apex, as the fusiform deck does. Every
point of the 20x20 (mu_in, cc_in) grid is run again, and the four summary_asym_*.csv of
each folder are rewritten with the converged points only (the format of run_febio_many.py
and densify_runner.py). The previous files are in ARCHIVE/2026-09-27_pre-sac-apex/.

    python3 rerun_asym_apex.py pilot      # 4 points of nu=3, s_z=1.0, no file written
    python3 rerun_asym_apex.py all        # the 15 configurations, 4 in parallel
"""
import multiprocessing as mp
import os
import shutil
import sys
import time

import numpy as np

import densify_runner as dr

GRID = np.linspace(0.0, 1.0, 20)
CONFIGS = [("asym", nu, zod) for (nu, zod) in dr._FOLDERS]
PROG = os.path.join(dr.PX, "rerun_asym_apex_progress.txt")


def run_config(type_, nu, zod, points=None, write=True):
    tag = f"{type_},{nu:g},{zod:g}"
    tree, root, mat, outs = dr.configure(type_, nu, zod)
    fdir = os.path.join(dr.PX, dr.folder_of(nu, zod))
    keys = {k: [f for f in outs if part in f and side in f][0]
            for k, part, side in (("elem_A", "elem", "-A"), ("elem_B", "elem", "-B"),
                                  ("node_A", "node", "-A"), ("node_B", "node", "-B"))}
    workdir = os.path.join(dr.PX, "_work", f"apex_{type_}_{nu:g}_{zod:g}")
    os.makedirs(workdir, exist_ok=True)
    pts = points or [(mu, cc) for mu in GRID for cc in GRID]
    rows = {k: [] for k in keys}
    t0, conv = time.time(), 0
    for k, (mu, cc) in enumerate(pts, 1):
        r = dr.run_point(tree, mat, outs, mu, cc, workdir)
        if r is not None:
            conv += 1
            for key, fn in keys.items():
                rows[key].append([mu, cc, *r[fn]])
        if k % 50 == 0 or k == len(pts):
            print(f"  {tag}: {k}/{len(pts)} run, {conv} conv, {(time.time() - t0) / 60:.1f} min",
                  flush=True)
    shutil.rmtree(workdir, ignore_errors=True)
    if write:
        for key, data in rows.items():
            with open(os.path.join(fdir, f"summary_{type_}_{key}.csv"), "w") as f:
                for row in data:
                    f.write(",".join("%.18e" % v for v in row) + "\n")
        with open(PROG, "a") as f:
            f.write(f"{tag},{conv}\n")
    print(f"DONE {tag}: {conv}/{len(pts)} converged, {(time.time() - t0) / 60:.1f} min", flush=True)
    return rows


if __name__ == "__main__":
    if sys.argv[1] == "pilot":
        rows = run_config("asym", 3.0, 1.0, points=[(0.8421052631578947, 0.0), (1.0, 1.0),
                                                   (0.0, 1.0), (0.5263157894736842, 0.8947368421052632)],
                          write=False)
        for a, n in zip(rows["elem_A"], rows["node_A"]):
            print("mu=%.3f cc=%.3f elem=%d J=%.4f node=%d" % (a[0], a[1], a[2], a[6], n[2]))
    elif sys.argv[1] == "all":
        done = set()
        if os.path.exists(PROG):
            done = {",".join(l.split(",")[:3]) for l in open(PROG) if l.strip()}
        todo = [c for c in CONFIGS if f"{c[0]},{c[1]:g},{c[2]:g}" not in done]
        print(f"RERUN ASYM APEX: {len(todo)}/{len(CONFIGS)} configs, {dr.NWORKERS} parallel x "
              f"{dr.THREADS} threads", flush=True)
        with mp.Pool(dr.NWORKERS) as pool:
            pool.starmap(run_config, todo, chunksize=1)
        print("ALL DONE", flush=True)
