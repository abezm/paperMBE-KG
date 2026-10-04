"""s1 copy (3 Oct 2026) of ../../../param-x/edge-retry/run_edge.py, executable ~/FEBio-stage1.

The cases are no longer listed by hand: they are the grid points of the 12 paper profiles of ../run_maps.py (PAPER_PROFILES) that did
not converge in the sweep (step 1, one retry) and have a 4-neighbour with lambda_r/lambda_ro <= 2.5, the rule of
make_numbers.py. Each is run again from t = 1 with step DT and max_retries 10, in its own folder; collect.py
(copied unchanged) writes results.json.

    python3 run_edge_s1.py list        # print the cases only
    python3 run_edge_s1.py [DT]        # run them (DT = 0.1 by default)
"""
import csv, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "param-x"))
sys.path.insert(0, os.path.join(HERE, ".."))
import densify_runner as D          # noqa: E402
from run_maps import PAPER_PROFILES as PROFILES      # noqa: E402

D.FEBIO = __import__("os").environ.get("FEBIO", "febio4")
D.PX = os.path.join(HERE, "..")
LAM_MAX = 2.5
STEP = 1.0 / 19.0

rows = list(csv.DictReader(open(os.path.join(HERE, "..", "MASTER_PARAMETRIC_STUDY.csv"))))
CASES = []
for (nu, sz) in PROFILES:
    for typ in ("axis", "asym"):
        pts = {(round(float(r["x_elastin_deg"]) / STEP), round(float(r["y_collagen_deg"]) / STEP)): float(r["z1_lamb_r"])
               for r in rows if r["type"] == typ and float(r["z3_nu"]) == nu and abs(float(r["z4_zod"]) - sz) < 1e-9}
        assert pts, (typ, nu, sz)
        for i in range(20):
            for j in range(20):
                if (i, j) in pts:
                    continue
                nb = [pts.get(q) for q in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1))]
                if any(v is not None and v <= LAM_MAX for v in nb):
                    CASES.append((typ, nu, sz, i, j))


def run(case, DT):
    typ, nu, sz, i, j = case
    tree, root, mat, outs = D.configure(typ, nu, sz)
    c = root.find("Control")
    # stage I (t = 0..1) keeps the sweep's single step; the smaller step applies from t = 1
    ts = c.find("time_stepper")
    ts.find("max_retries").text = "10"; ts.find("dtmin").text = str(DT / 64)
    for p in root.findall('.//load_controller[@id="2"]//point'):
        t0 = float(p.text.split(",")[0])
        p.text = f"{t0:g},{1.0 if t0 == 0 else DT}"
    wd = os.path.join(HERE, f"{typ}_nu{nu:g}_sz{sz:g}_{i}_{j}_dt{DT:g}_b")
    if os.path.exists(wd):
        raise SystemExit(f"refusing to overwrite {wd}")
    os.makedirs(wd)
    mat.find('mu_in').text = '%.18e' % (1 - i / 19); mat.find('cc_in').text = '%.18e' % (1 - j / 19)
    tree.write(f"{wd}/_run.feb", xml_declaration=True)
    env = {**os.environ, "OMP_NUM_THREADS": "4", "MKL_NUM_THREADS": "4"}
    subprocess.run([D.FEBIO, "-i", "_run.feb"], stdout=subprocess.DEVNULL, cwd=wd, env=env)
    ok = "N O R M A L   T E R M I N A T I O N" in open(f"{wd}/_run.log").read()
    print(os.path.basename(wd), "CONVERGED" if ok else "FAILED", flush=True)


if __name__ == "__main__":
    print(f"{len(CASES)} cases:", *CASES, sep="\n  ")
    if sys.argv[1:] != ["list"]:
        DT = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
        with ThreadPoolExecutor(4) as ex:
            list(ex.map(lambda k: run(k, DT), CASES))
