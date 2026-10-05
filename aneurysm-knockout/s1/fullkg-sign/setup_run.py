"""Full-KG of the Figure 4 fusiform lesion (elastin loss only) with the s1 code, s_load = 100 (author 5 Oct 2026).
Deck = ../param-x/matrix-fus-el.feb (MBE) with MBE = 0, gr_end = 101 and the control of
../../alp1/fig4-full/sl100-b (max_retries 10, opt_iter 10, dtmax 1, max_refs 100; dtmin 0.01 here).
a: Fig. 4 parameters; b: unit elastin deposition stretches; c: mass fractions 0.45/0.45/0.10."""
import os, re, subprocess
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "param-x", "matrix-fus-el.feb")
FEBIO = os.environ.get("FEBIO", __import__("os").environ.get("FEBIO", "febio4"))
import sys
DRY = "dry" in sys.argv
CASES = {"a": {}, "b": {"Ge_t": "1.0", "Ge_z": "1.0"}, "c": {"f_e": ".45", "f_c": ".45"}}

def sub(s, old, new):
    assert s.count(old) == 1, old
    return s.replace(old, new)

def make(tag):
    s = open(SRC).read()
    s = sub(s, "<time_steps>11</time_steps>", "<time_steps>501</time_steps>")
    s = sub(s, "<max_retries>1</max_retries>", "<max_retries>10</max_retries>")
    s = sub(s, "<opt_iter>5</opt_iter>", "<opt_iter>10</opt_iter>")
    s = sub(s, "<dtmin>1</dtmin>", "<dtmin>0.01</dtmin>")
    s = sub(s, '<dtmax lc="2"></dtmax>', "<dtmax>1</dtmax>")
    s = sub(s, "<max_refs>15</max_refs>", "<max_refs>100</max_refs>")
    s = sub(s, "<gr_end>11.0</gr_end>", "<gr_end>101.0</gr_end>")
    s = sub(s, "<MBE>1</MBE>", "<MBE>0</MBE>")
    for k, v in CASES[tag].items():
        s = re.sub(rf"<{k}>[^<]*</{k}>", f"<{k}>{v}</{k}>", s, count=1)
    d = os.path.join(HERE, tag); os.makedirs(d, exist_ok=True)
    open(os.path.join(d, "fus-el-full.feb"), "w").write(s)
    if DRY:
        return
    env = {**os.environ, "OMP_NUM_THREADS": "4", "MKL_NUM_THREADS": "4"}
    subprocess.run([FEBIO, "-i", "fus-el-full.feb"], cwd=d, env=env, stdout=open(os.path.join(d, "stdout.txt"), "w"))

with ThreadPoolExecutor(3) as ex:
    list(ex.map(make, CASES))
