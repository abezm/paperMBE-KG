#!/usr/bin/env python3
"""Referee run (author 9 Oct 2026): circumferential refinement x2 of the four Figure 4 lesions
(s1/param-x/matrix-{fus,sac}-{el,elcc}.feb, read only). Every circumferential element split in two:
the 40 angular node layers (9 deg) get their mid-layers (4.5 deg), 80 layers, 40 hex27 around; radial and axial
layers unchanged (1 x 40 x 20). Uses build/emit of gen.py unchanged; only the angular layers are refined."""
import numpy as np
import gen

_layers = gen.layers

def layers_th2(t):
    rv, tv, zv = _layers(t)
    assert len(tv) == 40 and np.allclose(np.diff(tv), 9.0)
    mid = 0.5 * (tv + np.append(tv[1:], tv[0] + 360.0))
    return rv, np.sort(np.concatenate([tv, mid % 360.0])), zv

gen.layers = layers_th2

if __name__ == "__main__":
    for case in ("fus-el", "fus-elcc", "sac-el", "sac-elcc"):
        src = gen.PX / f"matrix-{case}.feb"
        _, _, z0 = _layers(src.read_text())
        assert len(z0) == 41
        gen.emit(src, gen.HERE / "th2" / f"{case}.feb", z0)
