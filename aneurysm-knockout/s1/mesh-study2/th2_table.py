"""Table of th2 (circumferential x2) against gate, same quantities as RESULT.md: J max/min over elements at t = 11,
apex radius ratio = max over inner-wall nodes at z = z_center of r(t = 11)/r(t = 1)."""
import re, numpy as np
def rd(f):
    out = {}; t = None
    for l in open(f):
        m = re.match(r"\*Time\s*=\s*([\d.eE+-]+)", l)
        if m: t = float(m.group(1)); out[t] = []; continue
        if l[0].isdigit(): out[t].append([float(v) for v in l.split(",")])
    return {k: np.array(v) for k, v in out.items()}
for c in ("fus-el", "fus-elcc", "sac-el", "sac-elcc"):
    row = []
    for d in ("gate", "th2"):
        J = rd(f"{d}/{c}-ALL.txt"); X = rd(f"{d}/{c}-INNER.txt")
        j = J[11.0][:, 1]; a, b = X[1.0], X[11.0]
        k = np.abs(a[:, 3] - 7.5) < 1e-6
        ratio = (np.hypot(b[k, 1], b[k, 2]) / np.hypot(a[k, 1], a[k, 2])).max()
        row.append((j.max(), j.min(), ratio))
    (g, h) = row
    print(f"{c:9s} Jmax {g[0]:.4f}/{h[0]:.4f} ({100*(h[0]/g[0]-1):+.2f}%)  Jmin {g[1]:.4f}/{h[1]:.4f} ({100*(h[1]/g[1]-1):+.2f}%)  ratio {g[2]:.4f}/{h[2]:.4f} ({100*(h[2]/g[2]-1):+.2f}%)")
