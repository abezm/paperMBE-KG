#!/usr/bin/env python3
"""Copies the decks into run folders of tancheck/ and adds the finite-difference check parameters.
Sources are read only (s1/ is never written)."""
from pathlib import Path
S1 = Path(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "../s1"))
HERE = Path(__file__).resolve().parent
UNI = S1 / "alp1/uniax-lo1/timing-stop"
LES = S1 / "alp2/timing"

def make(folder, src, out, t, elem=1, gp=-1, extra=()):
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
    d.mkdir(exist_ok=True)
    (d / out).write_text(s)

make("a-full-t5",       UNI / "u_ds0.5.feb", "u_ds0.5.feb", 5.0)
make("a-full-t20",      UNI / "u_ds0.5.feb", "u_ds0.5.feb", 20.0)
make("a-full-t5-alp0",  UNI / "u_ds0.5.feb", "u_ds0.5.feb", 5.0,
     extra=[("<alp_update>1</alp_update>", "<alp_update>0</alp_update>")])
make("b-mbe-t5",        UNI / "conv-uniax-MBE.feb", "conv-uniax-MBE.feb", 5.0)
make("c-lesion-t10",    LES / "a_ds1.0.feb", "a_ds1.0.feb", 10.0, elem=10, gp=-1)
