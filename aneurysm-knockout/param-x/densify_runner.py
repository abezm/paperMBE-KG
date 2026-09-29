"""
Densify the (d_e, d_c) parametric grid from 20x20 -> 39x39 by inserting midpoints,
reusing the existing 400 converged points. Mapping (confirmed 2026-07-15):
  z_nu = nu, t_nu = nu, z_spread = zod (s_z), z_center = 7.5, alp_update = 1 (alph ON),
  asymm_aneurysm = 0 (axis) / 1 (asym). mu_in = 1 - d_elastin, cc_in = 1 - d_collagen.
Modes:
  validate <type> <nu> <zod> [n]  -> reproduce n existing points, compare J to the folder CSV
  densify  <type> <nu> <zod>      -> run the 1121 midpoints, append to the folder's summary CSVs
"""
import xml.etree.ElementTree as ET
import subprocess, numpy as np, os, sys, time

PX = os.path.dirname(os.path.abspath(__file__)) + "/"
FEBIO = __import__("os").environ.get("FEBIO", "febio4")
Z_CENTER = 7.5
DENS_OLD = 20          # existing grid: linspace(0,1,20)
DENS_NEW = 39          # refined grid: linspace(0,1,39); [::2] == the 20-grid
THREADS  = int(os.environ.get("DENSIFY_THREADS", "4"))   # OMP/MKL threads per FEBio solve
NWORKERS = int(os.environ.get("DENSIFY_WORKERS", "4"))   # configs run in parallel
# 4 workers x 4 threads = 16 cores; leaves ~4 for the desktop.

def folder_of(nu, zod):
    zs = ('%g' % zod)
    return f"nu={('%g'%nu)}_alph=ON_zod={zs}"

def configure(type_, nu, zod):
    base = f"{PX}/TAA-{type_}-MBE-BASE.feb"
    tree = ET.parse(base); root = tree.getroot()
    mat = root.find('.//material[@name="Material1"]')
    def setv(tag, val):
        n = mat.find(tag)
        if n is None: raise RuntimeError(f"missing tag <{tag}> in {base}")
        n.text = str(val)
    setv('z_nu', ('%g'%nu)); setv('t_nu', ('%g'%nu)); setv('z_spread', ('%g'%zod))
    setv('z_center', ('%g'%Z_CENTER)); setv('alp_update', 1)
    setv('asymm_aneurysm', 0 if 'axis' in type_ else 1)
    outs = [d.attrib['file'] for d in root.find('.//Output/logfile') if 'file' in d.attrib]
    return tree, root, mat, outs

def run_point(tree, mat, outs, mu, cc, workdir):
    """Run one point inside an ISOLATED workdir (FEBio writes its fixed-name output
    files there), so parallel same-type configs never clobber each other."""
    mat.find('mu_in').text = ('%.18e'%mu); mat.find('cc_in').text = ('%.18e'%cc)
    tree.write(os.path.join(workdir, "_run.feb"), xml_declaration=True)
    env = {**os.environ, "OMP_NUM_THREADS": str(THREADS), "MKL_NUM_THREADS": str(THREADS)}
    subprocess.run([FEBIO, "-i", "_run.feb"], stdout=subprocess.DEVNULL, cwd=workdir, env=env)
    log = os.path.join(workdir, "_run.log")
    ok = False
    try:
        with open(log) as f:
            if "N O R M A L   T E R M I N A T I O N" in f.read(): ok = True
    except FileNotFoundError:
        return None
    if not ok: return None
    res = {}
    for fn in outs:
        d = np.loadtxt(os.path.join(workdir, fn), delimiter=',', comments='*')
        res[fn] = d[-1, :] if d.ndim > 1 else d
    return res

def A_elem(outs):  return [f for f in outs if 'elem' in f and '-A' in f][0]

def validate(type_, nu, zod, n=6):
    tree, root, mat, outs = configure(type_, nu, zod)
    ref = np.loadtxt(f"{PX}/{folder_of(nu,zod)}/summary_{type_}_elem_A.csv",
                     delimiter=',', comments='#')
    idx = np.linspace(0, len(ref)-1, n).astype(int)
    workdir = f"{PX}/_work/_validate_{type_}"; os.makedirs(workdir, exist_ok=True)
    key = A_elem(outs)
    print(f"VALIDATE {type_} nu={nu} zod={zod}: reproducing {n} existing points")
    print(f"{'mu':>8}{'cc':>8}{'J_ref':>14}{'J_new':>14}{'rel_err':>12}")
    worst = 0.0
    for i in idx:
        mu, cc, Jref = ref[i,0], ref[i,1], ref[i,6]
        r = run_point(tree, mat, outs, mu, cc, workdir)
        if r is None:
            print(f"{mu:8.4f}{cc:8.4f}{Jref:14.6f}{'DID-NOT-CONVERGE':>26}"); worst = 9e9; continue
        Jnew = r[key][4]
        rel = abs(Jnew-Jref)/max(abs(Jref),1e-9)
        worst = max(worst, rel)
        print(f"{mu:8.4f}{cc:8.4f}{Jref:14.6f}{Jnew:14.6f}{rel:12.2e}")
    print(f"WORST relative error: {worst:.2e}  -> {'PASS' if worst < 1e-4 else 'FAIL'}")
    return worst

def midpoints():
    """The 1121 refined-grid points that are NOT on the existing 20x20 even sub-grid."""
    g = np.linspace(0.0, 1.0, DENS_NEW)          # 39 values; g[::2] == linspace(0,1,20)
    pts = []
    for i in range(DENS_NEW):
        for j in range(DENS_NEW):
            if not (i % 2 == 0 and j % 2 == 0):  # skip even-even (already run)
                pts.append((g[i], g[j]))          # (mu, cc)
    return pts

def converged_nodes(fEA):
    """Set of 20-grid (imu,icc) indices that converged (present in the existing summary).
    A point absent from the summary did not terminate normally = failed."""
    if not os.path.exists(fEA):
        return set()
    d = np.loadtxt(fEA, delimiter=',', comments='#')
    if d.ndim == 1: d = d[None, :]
    s = set()
    for mu, cc in d[:, 0:2]:
        fi, fj = mu * (DENS_OLD - 1), cc * (DENS_OLD - 1)   # to 20-grid index space
        ii, jj = round(fi), round(fj)
        if abs(fi - ii) < 1e-3 and abs(fj - jj) < 1e-3:     # genuine even-node only
            s.add((int(ii), int(jj)))
    return s

def midpoints_frontier(cset):
    """Refined-grid midpoints that TOUCH at least one converged existing node
    (>=1 even-node neighbour within +/-1 index). Midpoints whose every neighbour
    failed are deep in the white region and are skipped."""
    if not cset:               # nothing converged for this config: nothing to refine
        return []
    g = np.linspace(0.0, 1.0, DENS_NEW)
    pts = []
    for i in range(DENS_NEW):
        for j in range(DENS_NEW):
            if i % 2 == 0 and j % 2 == 0:     # existing node, reuse
                continue
            di = (-1, 1) if i % 2 else (0,)
            dj = (-1, 1) if j % 2 else (0,)
            touch = False
            for a in di:
                for b in dj:
                    I, J = i + a, j + b
                    if 0 <= I < DENS_NEW and 0 <= J < DENS_NEW:
                        if (I // 2, J // 2) in cset:
                            touch = True
            if touch:
                pts.append((g[i], g[j]))
    return pts

PROG = f"{PX}/densify_progress.txt"
def _done_set():
    if not os.path.exists(PROG): return set()
    return set(l.strip() for l in open(PROG) if l.strip())

def densify(type_, nu, zod, limit=None, pilot=False):
    tag = f"{type_},{nu:g},{zod:g}"
    if not pilot and tag in _done_set():
        print(f"SKIP {tag} (already done)"); return
    tree, root, mat, outs = configure(type_, nu, zod)
    fdir = f"{PX}/{folder_of(nu,zod)}"
    if not os.path.isdir(fdir):
        print(f"!! missing folder {fdir}; skipping {tag}"); return
    eA=[f for f in outs if 'elem' in f and '-A' in f][0]; eB=[f for f in outs if 'elem' in f and '-B' in f][0]
    nA=[f for f in outs if 'node' in f and '-A' in f][0]; nB=[f for f in outs if 'node' in f and '-B' in f][0]
    workdir = f"{PX}/_work/{type_}_{nu:g}_{zod:g}"; os.makedirs(workdir, exist_ok=True)
    if pilot:
        pts = midpoints()
    else:
        cset = converged_nodes(f"{fdir}/summary_{type_}_elem_A.csv")
        pts = midpoints_frontier(cset)
        skipped = (DENS_NEW*DENS_NEW - DENS_OLD*DENS_OLD) - len(pts)
        print(f"  {tag}: frontier-skip -> running {len(pts)}/1121 midpoints "
              f"({skipped} deep-white skipped)", flush=True)
    if limit: pts = pts[:limit]
    rows = {'elem_A':[], 'elem_B':[], 'node_A':[], 'node_B':[]}
    t0=time.time(); conv=0
    for k,(mu,cc) in enumerate(pts, 1):
        r = run_point(tree, mat, outs, mu, cc, workdir)
        if r is not None:
            conv += 1
            rows['elem_A'].append([mu,cc,*r[eA]]); rows['elem_B'].append([mu,cc,*r[eB]])
            rows['node_A'].append([mu,cc,*r[nA]]); rows['node_B'].append([mu,cc,*r[nB]])
        if k % 25 == 0 or k == len(pts):
            el=time.time()-t0; eta=el/k*(len(pts)-k)
            print(f"  {tag}: {k}/{len(pts)} run, {conv} conv, {el/60:.1f} min, ETA {eta/60:.1f} min", flush=True)
    import shutil; shutil.rmtree(workdir, ignore_errors=True)
    suffix = "_MID_pilot" if pilot else ""
    # append converged midpoints to the folder's summary CSVs (or *_pilot for a dry run)
    fmap = {'elem_A':f"summary_{type_}_elem_A", 'elem_B':f"summary_{type_}_elem_B",
            'node_A':f"summary_{type_}_node_A", 'node_B':f"summary_{type_}_node_B"}
    for key,base in fmap.items():
        out = f"{fdir}/{base}{suffix}.csv"
        with open(out, 'a') as f:
            for row in rows[key]:
                f.write(','.join('%.18e'%v for v in row) + '\n')
    if not pilot:
        with open(PROG,'a') as f: f.write(tag+'\n')
    print(f"DONE {tag}: {conv}/{len(pts)} converged, appended{' (PILOT)' if pilot else ''}, {(time.time()-t0)/60:.1f} min", flush=True)

# the 15 unique (nu, zod) folders the paper's figures use, x 2 types
_FOLDERS = list(dict.fromkeys(
    [(nu,1.0) for nu in (1,2,3,4,5,6)] +                       # Figs 10/11 nu grid
    [(3.0,z) for z in (1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8)]))  # Figs 8/9 s_z grid
ALL_CONFIGS = [(t,nu,zod) for (nu,zod) in _FOLDERS for t in ('axis','asym')]

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "validate":
        type_, nu, zod = sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
        n = int(sys.argv[5]) if len(sys.argv) > 5 else 6
        validate(type_, nu, zod, n)
    elif mode == "pilot":
        densify(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]),
                limit=int(sys.argv[5]) if len(sys.argv)>5 else 12, pilot=True)
    elif mode == "densify":
        densify(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]))
    elif mode == "hero":
        import multiprocessing as mp
        hero = [('axis',3.0,1.0), ('asym',3.0,1.0)]
        todo = [c for c in hero if f"{c[0]},{c[1]:g},{c[2]:g}" not in _done_set()]
        print(f"DENSIFY HERO: {len(todo)} configs, {NWORKERS} parallel x {THREADS} threads", flush=True)
        with mp.Pool(NWORKERS) as pool:
            pool.starmap(densify, todo, chunksize=1)
        print("HERO CONFIGS DONE", flush=True)
    elif mode == "all":
        import multiprocessing as mp
        todo = [c for c in ALL_CONFIGS if f"{c[0]},{c[1]:g},{c[2]:g}" not in _done_set()]
        print(f"DENSIFY ALL: {len(todo)}/{len(ALL_CONFIGS)} configs remaining, "
              f"{NWORKERS} parallel x {THREADS} threads", flush=True)
        with mp.Pool(NWORKERS) as pool:
            pool.starmap(densify, todo, chunksize=1)
        print("ALL CONFIGS DONE", flush=True)
