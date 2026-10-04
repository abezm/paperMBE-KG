#!/usr/bin/env python3
"""s1/profile-fig (2026-10-04): decks for the profile figure (Fig. 6 at fixed d_e = d_c, profile varied).
Copies of ../severity-fig/lad-fus-d20p000-v5.feb (fusiform, 20 %) and lad-sac-d37p000-lev2.feb (saccular, 37 %);
only z_spread, z_nu, t_nu and the output file names change. nu sets z_nu and t_nu together, as in the maps
(param-x/densify_runner.configure). The nu = 5, s_z = 3 mm panels reuse the two source runs."""
import os, re
HERE = os.path.dirname(os.path.abspath(__file__))
SRC = {"fus": "lad-fus-d20p000-v5", "sac": "lad-sac-d37p000-lev2"}
CASES = [(m, 5, sz) for m in SRC for sz in (1.0, 1.2, 1.8, 2.4, 2.6)] + [(m, nu, 3.0) for m in SRC for nu in (1, 3, 6, 7)]

def tag(m, nu, sz):
    return f"prof-{m}-nu{nu}-sz{sz:.1f}".replace(".", "p")

for m, nu, sz in CASES:
    src = SRC[m]
    txt = open(os.path.join(HERE, "..", "severity-fig", src + ".feb")).read()
    assert txt.count(src) == 6
    for key, old, new in (("z_spread", "3.0", f"{sz:g}"), ("z_nu", "5", f"{nu:g}"), ("t_nu", "5", f"{nu:g}")):
        pat = f"<{key}>{old}</{key}>"
        assert txt.count(pat) == 1, (key, src)
        txt = txt.replace(pat, f"<{key}>{new}</{key}>")
    out = tag(m, nu, sz)
    open(os.path.join(HERE, out + ".feb"), "w").write(txt.replace(src, out))
    print(out)
