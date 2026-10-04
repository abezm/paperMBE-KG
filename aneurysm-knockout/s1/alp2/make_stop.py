#!/usr/bin/env python3
"""alp2/ (3 Oct 2026): stopping time of the timed aneurysm runs from the untimed run speedup-fair/stop/find/
(J at the apex element 10 changes by less than 1e-4 (J - 1) over the last s_load/20 = 5 days, s >= s_load = 100 days;
the criterion of make_numbers.py), then the timed decks ending at that time and the timing script.
usage (cwd = alp2): python3 make_stop.py
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
SF = HERE / "speedup-fair"


def blocks(path):
    out, t = [], None
    for ln in open(path):
        s = ln.strip()
        if s.startswith("*Time"):
            t = float(s.split("=")[1]); out.append((t, {}))
        elif s and s[0].isdigit() and t is not None:
            p = s.split(","); out[-1][1][int(p[0])] = [float(v) for v in p[1:]]
    return out


J = {round(t, 6): d[10][3] for t, d in blocks(SF / "stop" / "find" / "elem-fus-sl100-2st-ds5.txt")}
stop = next(t - 1.0 for t in sorted(J) if t - 1.0 >= 100.0 - 1e-9 and round(t - 5.0, 6) in J
            and abs(J[t] - J[round(t - 5.0, 6)]) < 1e-4 * (J[t] - 1.0))
n5 = round(stop / 5.0)
assert abs(n5 * 5.0 - stop) < 1e-9
print(f"stopping time s = {stop:g} days ({n5} steps of 5 days, {round(stop)} of 1 day)")
s5 = (SF / "fus-sl100-2st-ds5.feb").read_text(encoding="ISO-8859-1")
s1 = (SF / "fus-sl100-2st-ds1-130.feb").read_text(encoding="ISO-8859-1")
assert s5.count("<time_steps>26</time_steps>") == 1 and s1.count("<time_steps>130</time_steps>") == 1
for name, txt in (("fus-sl100-2st-ds5-stop.feb", s5.replace("<time_steps>26</time_steps>", f"<time_steps>{n5}</time_steps>")),
                  ("fus-sl100-2st-ds1-stop.feb", s1.replace("<time_steps>130</time_steps>", f"<time_steps>{round(stop)}</time_steps>")),
                  ("mbe-1p10.feb", (SF / "mbe-1p10.feb").read_text(encoding="ISO-8859-1"))):
    out = SF / "stop" / name
    assert not out.exists(), out
    out.write_text(txt, encoding="ISO-8859-1")
src = (HERE.parent / "alp1" / "speedup-fair" / "stop" / "timing_stop.py").read_text()
for a, b in (('HERE / "fus-sl100-2st-ds5-155.feb"', 'HERE / "fus-sl100-2st-ds5-stop.feb"'),
             ('HERE / "fus-sl100-2st-ds1-155.feb"', 'HERE / "fus-sl100-2st-ds1-stop.feb"'),
             ("stop/ (2 Oct 2026): the full-KG runs end at the stopping time s = 155 days",
              f"alp2/ (3 Oct 2026): severities fus cc_in 0.95, sac 0.9; the full-KG runs end at the stopping time s = {stop:g} days")):
    assert src.count(a) == 1, a
    src = src.replace(a, b)
(SF / "stop" / "timing_stop.py").write_text(src)
