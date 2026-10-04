#!/usr/bin/env python3
"""alp2/ (3 Oct 2026): the rate-comparison runs of alp1/ with the lesion severities adopted by the author:
fusiform 20 % elastin loss + 5 % loss of collagen cross-linking (cc_in 0.95), saccular 65 % elastin loss + 10 %
(cc_in 0.9). Each input file is the alp1 one with only <cc_in> changed; the folder structure mirrors alp1/.
Sources are only read; a destination that exists stops the script.  usage (cwd = alp2): python3 setup.py
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
A1 = HERE.parent / "alp1"
DECKS = {
    "full-KG/apex": ["rate-fusiform-MBE.feb", "rate-saccular-MBE.feb", "rate-saccular-full-sload10.feb"],
    "full-KG/apex/to5": ["rate-fusiform-full-sload10.feb", "rate-fusiform-full-sload100.feb",
                         "rate-saccular-full-sload100.feb"],
    "full-KG/apex/ds-aneu-t1": ["fus_ds001.feb", "fus_ds005.feb", "fus_ds01.feb", "sac_ds001.feb", "sac_ds01.feb"],
    "speedup-fair": ["mbe-1p10.feb", "fus-sl100-2st-ds5.feb", "fus-sl100-2st-ds1-130.feb", "hyp-1.feb"],
    "timing": ["a_ds1.0.feb", "conv-aneu-MBE.feb"],
}
CC = {("0", ".8"): "0.95", ("1", ".35"): "0.9"}          # (asymm_aneurysm, mu_in) -> cc_in
n = 0
for sub, names in DECKS.items():
    (HERE / sub).mkdir(parents=True, exist_ok=True)
    for name in names:
        src, out = A1 / sub / name, HERE / sub / name
        assert not out.exists(), out
        s = src.read_text(encoding="ISO-8859-1")
        asym = s.split("<asymm_aneurysm>")[1].split("<")[0].strip()
        mu = s.split("<mu_in>")[1].split("<")[0].strip()
        assert s.count("<cc_in>1.</cc_in>") == 1, src
        out.write_text(s.replace("<cc_in>1.</cc_in>", f"<cc_in>{CC[(asym, mu)]}</cc_in>"), encoding="ISO-8859-1")
        n += 1
print(f"{n} input files written")
