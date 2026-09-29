#!/usr/bin/env python3
"""Close referee item M8 ON THE DECKS THAT ACTUALLY PRODUCE THE PAPER'S CLINICAL NUMBERS.

The first mesh study refined `rate-fusiform-MBE` -- the RATE-COMPARISON config: unit G^e,
J ~ 1.02, a tiny deformation.  The paper's field maps and parametric study do NOT come from
there: they run the CALIBRATED material (G^e = 1.9/1.62) and reach J = 1.65 (fusiform) /
1.91 (saccular).  Mesh convergence at J=1.02 says little about J=1.9, where the elements are
actually strained.  So refine THOSE.

The clinical decks were verified to carry the IDENTICAL mesh (same 4920 nodes, same
coordinates, same NodeSet/Surface names) as the deck gen_mesh.py was built against, so the
generator applies unchanged with the source swapped.

VALIDATION GATE (same discipline as before): the regenerated nr=1 must reproduce the
published apex values, J_max = 1.654 (fusiform) and 1.910 (saccular).  If it does not,
STOP -- the refined meshes mean nothing.
"""
import re, subprocess, sys, time
from pathlib import Path

import gen_mesh

FEBIO = __import__("os").environ.get("FEBIO", "febio4")
HERE = Path(__file__).parent
PX = HERE.parent / "param-x"

CASES = [("fus", PX / "matrix-fus-el.feb", 1.654),
         ("sac", PX / "matrix-sac-el.feb", 1.910)]
TOL = 0.005          # the published values are quoted to 3 decimals


def jmax(tag):
    f = HERE / f"{tag}_{tag}_ALLJ.txt"
    block = f.read_text().split("*Step")[-1]
    return max(float(l.split(",")[1]) for l in block.splitlines() if l and l[0].isdigit())


def run(feb):
    log = HERE / f"{feb.stem}.log"
    t0 = time.perf_counter()
    subprocess.run([FEBIO, "-i", str(feb), "-o", str(log)], cwd=HERE, capture_output=True, text=True)
    ok = "N O R M A L" in log.read_text()
    return ok, time.perf_counter() - t0


assert gen_mesh.verify_template(), "TEMPLATE no valida: ABORTAR"
print()

for name, src, published in CASES:
    gen_mesh.SRC = src                       # same mesh, different material/insult
    print(f"=== {name}  ({src.name}, G^e=1.9/1.62 calibrado)   publicado: J_max = {published}")
    results = []
    for nr in (1, 2, 3):
        tag = f"{name}_nr{nr}"
        gen_mesh.emit(nr, HERE / f"{tag}.feb")
        ok, w = run(HERE / f"{tag}.feb")
        if not ok:
            print(f"  nr={nr}: NO TERMINO NORMAL -- abortando {name}")
            break
        j = jmax(tag)
        results.append((nr, j))
        if nr == 1:
            d = abs(j - published)
            print(f"  COMPUERTA nr=1: J_max = {j:.5f}  vs publicado {published}  "
                  f"({'PASA' if d <= TOL else 'FALLA'})")
            if d > TOL:
                print("  La malla regenerada NO reproduce el resultado publicado. ABORTAR: "
                      "los refinamientos no significarian nada.")
                sys.exit(1)
        else:
            print(f"  nr={nr}: J_max = {j:.5f}   dif vs nr=1 = {(j/results[0][1]-1)*100:+.3f}%   ({w:.0f}s)")
    print()

print("Si las diferencias son <1%, los numeros clinicos del paper son independientes de malla y M8 cierra.")
