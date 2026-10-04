#!/usr/bin/env python3
"""Saccular timing (2026-10-03, simulated-review point 15): the fusiform timing decks of ../speedup-fair/ with the
saccular lesion of setup.py (mu_in .35, cc_in 0.9, asymm_aneurysm 1); output files renamed fus -> sac.
Sources are only read; an existing destination stops the script.  usage (cwd = this folder): python3 setup_sac.py
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SF = HERE.parent / "speedup-fair"
PAIRS = [(SF / "stop" / "find" / "fus-sl100-2st-ds5-500.feb", HERE / "find" / "sac-sl100-2st-ds5-500.feb"),
         (SF / "mbe-1p10.feb", HERE / "mbe-1p10-sac.feb"),
         (SF / "fus-sl100-2st-ds5.feb", HERE / "sac-sl100-2st-ds5.feb")]
SUBS = [("<mu_in>.8</mu_in>", "<mu_in>.35</mu_in>"), ("<cc_in>0.95</cc_in>", "<cc_in>0.9</cc_in>"),
        ("<asymm_aneurysm>0</asymm_aneurysm>", "<asymm_aneurysm>1</asymm_aneurysm>")]
for src, out in PAIRS:
    assert not out.exists(), out
    s = src.read_text(encoding="ISO-8859-1")
    for a, b in SUBS:
        assert s.count(a) == 1, (src, a)
        s = s.replace(a, b)
    s = re.sub(r'file="([^"]*)fus-', r'file="\1sac-', s)
    assert "fus-" not in "".join(re.findall(r'file="([^"]+)"', s)), out
    out.parent.mkdir(exist_ok=True)
    out.write_text(s, encoding="ISO-8859-1")
    print(out.relative_to(HERE))
