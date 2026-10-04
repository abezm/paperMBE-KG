#!/usr/bin/env python3
"""s1/ (3 Oct 2026): the parametric maps with Stage I as in Section 5.1 (executable ~/FEBio-stage1).

Copy of ../../param-x/rerun_asym_apex.py for both morphologies, on the 20 x 20 grid. Profiles:
s_z = 1.0, 1.4, 1.8, 2.2, 2.6, 2.8 mm at nu = 3 (author, 3 Oct 2026; kept as data, not in the paper since 4 Oct);
the same s_z at nu = 5, the value of latorre2020numerical (author, 4 Oct 2026); nu = 1..6 at s_z = 3 mm.
PAPER_PROFILES = the 12 shown in the paper (nu = 5 s_z sweep + nu sweep);
of (mu_in, cc_in) = linspace(0, 1, 20). Decks from TAA-{axis,asym}-MBE-BASE.feb of this folder (copied
unchanged from ../../param-x/) through densify_runner.configure. Output, per configuration, in
<folder>/: summary_{type}_{elem,node}_{A,B}.csv (converged points only, the format of the earlier runs)
and points_{type}.csv (mu_in, cc_in, converged, wall seconds for every point). Progress in progress.txt.

    python3 run_maps.py pilot     # axis, nu = 3, s_z = 1.0 (one configuration)
    python3 run_maps.py all       # every configuration not yet in progress.txt, 4 in parallel
"""
import multiprocessing as mp
import os
import shutil
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "param-x"))
import densify_runner as dr  # noqa: E402

dr.FEBIO = __import__("os").environ.get("FEBIO", "febio4")
dr.PX = HERE                         # decks (TAA-*-MBE-BASE.feb) and outputs in s1/param-x/
GRID = np.linspace(0.0, 1.0, 20)
SZ_RUN = (1.0, 1.4, 1.8, 2.2, 2.6, 2.8)          # s_z sweeps computed (nu = 3 and nu = 5)
SZ = (1.0, 1.4, 1.8, 2.2, 2.6, 3.0)              # s_z shown in the paper at nu = 5 (author, 4 Oct 2026: up to 3 mm;
                                                 # the 3.0 map is the nu sweep one, the 2.8 map stays as data)
PAPER_PROFILES = list(dict.fromkeys([(5.0, z) for z in SZ] + [(float(n), 3.0) for n in range(1, 7)]))
PROFILES = list(dict.fromkeys([(3.0, z) for z in SZ_RUN] + [(5.0, z) for z in SZ_RUN] + [(float(n), 3.0) for n in range(1, 7)]))
# saccular maps (~4 h) before fusiform ones (~2.5 h), so the last ones to finish are short
CONFIGS = [(t, nu, zod) for t in ("asym", "axis") for (nu, zod) in PROFILES]
assert len(CONFIGS) == 36 and len(set(CONFIGS)) == 36 and len(PAPER_PROFILES) == 11
PROG = os.path.join(HERE, "progress.txt")


def run_config(type_, nu, zod):
    tag = f"{type_},{nu:g},{zod:g}"
    tree, root, mat, outs = dr.configure(type_, nu, zod)
    fdir = os.path.join(HERE, dr.folder_of(nu, zod))
    os.makedirs(fdir, exist_ok=True)
    keys = {k: [f for f in outs if part in f and side in f][0]
            for k, part, side in (("elem_A", "elem", "-A"), ("elem_B", "elem", "-B"),
                                  ("node_A", "node", "-A"), ("node_B", "node", "-B"))}
    workdir = os.path.join(HERE, "_work", f"{type_}_{nu:g}_{zod:g}")
    os.makedirs(workdir, exist_ok=True)
    rows = {k: [] for k in keys}
    pts = []
    t0, conv = time.time(), 0
    for k, (mu, cc) in enumerate([(mu, cc) for mu in GRID for cc in GRID], 1):
        t1 = time.time()
        r = dr.run_point(tree, mat, outs, mu, cc, workdir)
        pts.append((mu, cc, int(r is not None), time.time() - t1))
        if r is not None:
            conv += 1
            for key, fn in keys.items():
                rows[key].append([mu, cc, *r[fn]])
        if k % 50 == 0:
            print(f"  {tag}: {k}/400 run, {conv} conv, {(time.time() - t0) / 60:.1f} min", flush=True)
    shutil.rmtree(workdir, ignore_errors=True)
    for key, data in rows.items():
        with open(os.path.join(fdir, f"summary_{type_}_{key}.csv"), "w") as f:
            for row in data:
                f.write(",".join("%.18e" % v for v in row) + "\n")
    with open(os.path.join(fdir, f"points_{type_}.csv"), "w") as f:
        f.write("mu_in,cc_in,converged,wall_s\n")
        for p in pts:
            f.write("%.18e,%.18e,%d,%.2f\n" % p)
    with open(PROG, "a") as f:
        f.write(f"{tag},{conv},{(time.time() - t0) / 60:.1f} min\n")
    print(f"DONE {tag}: {conv}/400 converged, {(time.time() - t0) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    if sys.argv[1] == "pilot":
        run_config("axis", 3.0, 1.0)
    elif sys.argv[1] == "all":
        done = set()
        if os.path.exists(PROG):
            done = {",".join(l.split(",")[:3]) for l in open(PROG) if l.strip()}
        todo = [c for c in CONFIGS if f"{c[0]},{c[1]:g},{c[2]:g}" not in done]
        print(f"MAPS s1: {len(todo)}/{len(CONFIGS)} configs, {dr.NWORKERS} parallel x {dr.THREADS} threads", flush=True)
        with mp.Pool(dr.NWORKERS) as pool:
            pool.starmap(run_config, todo, chunksize=1)
        print("ALL DONE", flush=True)
