#!/usr/bin/env python3
"""Round rABC-8 (2026-10-04): finite-difference tangent check on the lesions not covered by setup.py.
Saccular full-KG, saccular MBE-KG, fusiform MBE-KG; apex element 10; all integration points.
Sources are read only; uses make() of setup.py logic (copied, no s1/ writes)."""
from pathlib import Path
S1 = Path(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "../s1/alp2"))
HERE = Path(__file__).resolve().parent

def make(folder, src, out, t, elem=10, gp=-1, extra=()):
    s = src.read_text()
    a = "<alp_update>1</alp_update>"
    assert s.count(a) == 1, src
    add = (f"{a}\n\t\t\t<fd_check>1</fd_check>\n\t\t\t<fd_time>{t}</fd_time>\n"
           f"\t\t\t<fd_elem>{elem}</fd_elem>\n\t\t\t<fd_gp>{gp}</fd_gp>\n\t\t\t<fd_eps>1e-6</fd_eps>")
    s = s.replace(a, add)
    for old, new in extra:
        assert s.count(old) == 1, (src, old)
        s = s.replace(old, new)
    d = HERE / folder
    d.mkdir(exist_ok=False)
    (d / out).write_text(s)

# saccular full-KG, s_load = 100 d, ds = 5 d: Stage II times 6, 11, ...; stop after t = 11
make("d-sac-full-t11", S1 / "speedup-fair-sac/sac-sl100-2st-ds5.feb", "sac-full.feb", 11.0,
     extra=[("<time_steps>26</time_steps>", "<time_steps>2</time_steps>")])
# saccular MBE-KG: Stage II pseudo-time 1.1 ... 2.0; check at 1.5, stop there
make("e-sac-mbe-t1.5", S1 / "speedup-fair-sac/mbe-1p10-sac.feb", "sac-mbe.feb", 1.5,
     extra=[("<time_steps>10</time_steps>", "<time_steps>5</time_steps>")])
# fusiform MBE-KG: 20 steps of 0.1 from t = 0; check at 1.5, stop there
make("f-fus-mbe-t1.5", S1 / "timing/conv-aneu-MBE.feb", "fus-mbe.feb", 1.5,
     extra=[("<time_steps>20</time_steps>", "<time_steps>15</time_steps>")])
