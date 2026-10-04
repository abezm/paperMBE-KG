"""Collect the retries of run_edge.py (folders *_b) into results.json: converged or not, apex
lambda_r/lambda_ro and J at the last converged step, and the fraction of the insult ramp
(t = 1..11) at which lambda_r/lambda_ro first exceeds 2.5."""
import glob, json, numpy as np
R0 = 0.6468   # reference radius of the apex node, = y - uy in every log

def rd(f):
    out, t = [], None
    for l in open(f):
        if l.startswith("*Time"): t = float(l.split("=")[1])
        elif not l.startswith("*") and l.strip(): out.append([t] + [float(x) for x in l.split(",")])
    return np.array(out)

res = []
for d in sorted(glob.glob("*_b")):
    typ, nu, sz, i, j = d.split("_")[:5]
    n = rd(glob.glob(d + "/node-*-A_MBE.txt")[0]); e = rd(glob.glob(d + "/elem-*-A_MBE.txt")[0])
    lam = np.hypot(n[:, 2], n[:, 4]) / np.hypot(n[:, 2] - n[:, 3], n[:, 4] - n[:, 5])
    ok = "N O R M A L   T E R M I N A T I O N" in open(d + "/_run.log").read()
    k = int(np.argmax(lam > 2.5))
    res.append(dict(type=typ, nu=float(nu[2:]), sz=float(sz[2:]), i=int(i), j=int(j), converged=ok,
                    lam_last=float(lam[-1]), J_last=float(e[-1, 5]), ramp_last=float((n[-1, 0] - 1) / 10),
                    ramp_cross_2p5=float((n[k, 0] - 1) / 10)))
    print(res[-1])
json.dump(res, open("results.json", "w"), indent=1)
