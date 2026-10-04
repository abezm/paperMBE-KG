#!/usr/bin/env python3
"""Point 16 summary: J max, J min and apex radius ratio at t = 11 for gate/, z2/, L30/ (four lesions of Figure 4).
Checks that every gate run reproduces param-x (same mesh, rebuilt by gen.py). Writes summary.csv.
Apex: inner-wall node nearest to (z_center, insult direction +y); ratio = r(t = 11) / r(t = 1), as lambda_r/lambda_ro."""
import csv
import numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
PX = HERE.parent / "param-x"


def blocks(f):
    out = []
    for b in Path(f).read_text().split("*Step")[1:]:
        L = b.strip().splitlines()
        tm = float([l for l in L if "Time" in l][0].split("=")[1])
        out.append((tm, {int(l.split(",")[0]): [float(x) for x in l.split(",")[1:]] for l in L if l and l[0].isdigit()}))
    return out


rows = []
for case in ("fus-el", "fus-elcc", "sac-el", "sac-elcc"):
    ref = np.array([v[0] for v in blocks(PX / f"matrix-{case}-ALL.txt")[-1][1].values()])
    for d, zc in (("gate", 7.5), ("z2", 7.5), ("L30", 15.0)):
        A = blocks(HERE / d / f"{case}-ALL.txt")
        assert abs(A[-1][0] - 11.0) < 1e-9, (case, d, A[-1][0])
        J = np.array([v[0] for v in A[-1][1].values()])
        if d == "gate":
            assert abs(J.max() - ref.max()) < 1e-6 and abs(J.min() - ref.min()) < 1e-6, (case, J.max(), ref.max())
        I = blocks(HERE / d / f"{case}-INNER.txt")
        t1 = [dd for tm, dd in I if abs(tm - 1.0) < 1e-9][0]
        tf = I[-1][1]
        apex = min(t1, key=lambda n: (abs(t1[n][2] - zc), -t1[n][1] / np.hypot(t1[n][0], t1[n][1])))
        rows.append(dict(case=case, mesh=d, Jmax=J.max(), Jmin=J.min(),
                         apex_ratio=np.hypot(*tf[apex][:2]) / np.hypot(*t1[apex][:2])))
with open(HERE / "summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    for r in rows:
        w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})
print("summary.csv written; the four gates reproduce param-x")
