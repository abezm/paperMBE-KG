"""Short full-KG runs of deck a with the MBE degradation ramp (gr_end = 11, as ../param-x/matrix-fus-el.feb), so that
t = 2, 3 have the same degradation as the MBE run; k_parameter = 0.1, 1, 10, 100 (k * ds = k, ds = 1). Author 5 Oct 2026.
If the MBE root attracts the rate law, J at the apex tends to the MBE value (1.042 at t = 2) as k grows."""
import os, re, subprocess
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
FEBIO = os.environ.get("FEBIO", __import__("os").environ.get("FEBIO", "febio4"))
KS = ["0.1", "1", "10", "100"]

def sub(s, old, new):
    assert s.count(old) == 1, old
    return s.replace(old, new)

def run(k):
    s = open(os.path.join(HERE, "a", "fus-el-full.feb")).read()
    s = sub(s, "<time_steps>501</time_steps>", "<time_steps>3</time_steps>")
    s = sub(s, "<gr_end>101.0</gr_end>", "<gr_end>11.0</gr_end>")
    s = sub(s, "<k_parameter>0.1</k_parameter>", f"<k_parameter>{k}</k_parameter>")
    d = os.path.join(HERE, f"k{k}"); os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "fus-el-full.feb"), "w").write(s)
    env = {**os.environ, "OMP_NUM_THREADS": "4", "MKL_NUM_THREADS": "4"}
    subprocess.run([FEBIO, "-i", "fus-el-full.feb"], cwd=d, env=env, stdout=open(os.path.join(d, "stdout.txt"), "w"))

with ThreadPoolExecutor(4) as ex:
    list(ex.map(run, KS))
