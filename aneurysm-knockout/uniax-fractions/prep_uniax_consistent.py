#!/usr/bin/env python3
"""Prepare plot-ready data for the uniaxial figure from the CONSISTENT knockout
runs (MBE + full-KG s_load=1,10,100). Emits <case>_c.txt with columns:
    s/s_load , sigma_axial[kPa] , lambda/lambda_o , lambda_g   (AXIAL coupon)
Reads raw FEBio logfiles directly (numpy only, no matplotlib)."""
import numpy as np


def read_log(fn, ncol):
    t, out = None, []
    for line in open(fn):
        s = line.strip()
        if s.startswith("*Time"):
            t = float(s.split("=")[1])
        elif s and not s.startswith("*"):
            p = s.split(",")
            if len(p) >= ncol and t is not None:
                out.append([t] + [float(x) for x in p[1:ncol]])
    return np.array(out)


def case(elem_f, node_f):
    e = read_log(elem_f, 6)   # elemID, sx, sy, sz, J, Fgxx
    n = read_log(node_f, 4)   # nodeID, ux, uy, uz
    m = min(len(e), len(n))
    # AXIAL coupon (cube off the z axis, pulled along z): the axial component is z,
    # so the stress is sz (col 3) and the displacement uz (col 3), not sx / ux.
    return dict(t=e[:m, 0], sxx=e[:m, 3], J=e[:m, 4], lg=e[:m, 5], ux=n[:m, 3])


files = {
    "MBE":   ("elem1_MBE.txt",     "node7_MBE.txt",     10.0),
    "sp1":   ("elem1_sp1ko.txt",   "node7_sp1ko.txt",    1.0),
    "sp10":  ("elem1_sp10ko.txt",  "node7_sp10ko.txt",  10.0),
    "sp100": ("elem1_sp100ko.txt", "node7_sp100ko.txt", 100.0),
}
init_time = 1.0

mbe = case(*files["MBE"][:2])
i1 = int(np.argmin(np.abs(mbe["t"] - init_time)))
lam_o = 1.0 + mbe["ux"][i1]
print("lambda_o =", round(lam_o, 4))

for tag, (ef, nf, sl) in files.items():
    c = case(ef, nf)
    x = (c["t"] - init_time) / sl
    lam_ratio = (1.0 + c["ux"]) / lam_o
    data = np.column_stack([x, c["sxx"], lam_ratio, c["lg"]])
    np.savetxt(f"{tag}_c.txt", data, fmt="%.6g",
               header="s/s_load  sigma_xx  lam/lam_o  lambda_g")
    print(f"  {tag:5s} -> {tag}_c.txt  end: s/sload={x[-1]:.2f} "
          f"sxx={c['sxx'][-1]:.1f} lam/lamo={lam_ratio[-1]:.3f} lg={c['lg'][-1]:.3f}")
