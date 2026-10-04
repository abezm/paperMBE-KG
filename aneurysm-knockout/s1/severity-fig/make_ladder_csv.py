#!/usr/bin/env python3
"""s1/severity-fig/ladder.csv from the ten decks of this folder (the eight panels of Figure 6 and the next equal
losses, 21 % fusiform and 38 % saccular), with the columns that make_numbers.py reads: deck, d_e_percent,
d_c_percent, converged (normal termination and last time = 11), J_max_elements (largest J over the elements at the
last time, *-ALL.txt) and lambda_r_over_lambda_ro_node20 (radial position of node 20 at the last time over t = 1,
*-NODES.txt)."""
import csv, os, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))

def blocks(p):
    out, t, d = [], None, None
    for ln in open(p):
        ln = ln.strip()
        if ln.startswith("*Time"):
            if d is not None: out.append((t, d))
            t, d = float(ln.split("=")[1]), {}
        elif ln and ln[0].isdigit():
            q = ln.split(","); d[int(q[0])] = [float(x) for x in q[1:]]
    if d is not None: out.append((t, d))
    return out

rows = []
for f in sorted(os.listdir(HERE)):
    if not f.endswith(".feb"): continue
    deck = f[:-4]
    m = re.match(r"lad-(fus|sac)-d(\d+)p000", deck)
    d = int(m.group(2))
    log = open(os.path.join(HERE, deck + ".log")).read()
    nb, eb = blocks(os.path.join(HERE, deck + "-NODES.txt")), blocks(os.path.join(HERE, deck + "-ALL.txt"))
    conv = "N O R M A L   T E R M I N A T I O N" in log and abs(eb[-1][0] - 11.0) < 1e-9
    h = min(nb, key=lambda b: abs(b[0] - 1.0))[1]; fin = nb[-1][1]
    lam = np.hypot(*fin[20][:2]) / np.hypot(*h[20][:2])
    rows.append(dict(morphology={"fus": "fusiform", "sac": "saccular"}[m.group(1)], d_e_percent=d, d_c_percent=d,
                     deck=f, converged="yes" if conv else "no", t_last_converged=eb[-1][0],
                     J_max_elements=f"{max(v[0] for v in eb[-1][1].values()):.6f}",
                     lambda_r_over_lambda_ro_node20=f"{lam:.6f}"))
with open(os.path.join(HERE, "ladder.csv"), "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
for r in rows: print(r["deck"], r["converged"], r["J_max_elements"], r["lambda_r_over_lambda_ro_node20"])
