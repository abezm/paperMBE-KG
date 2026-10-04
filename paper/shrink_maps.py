"""Write reduced copies of the parametric maps (fig-param/*.png) to fig-param-s/, for
supplementary.tex. Each panel prints at 0.45 textwidth (~7.6 cm); at 1/3 of the original
width it keeps ~650 dpi at that size, and a 256-colour palette (no dithering) holds the
contour maps without visible change. The originals are not modified.

Run: python3 shrink_maps.py      (after the maps are regenerated)
     python3 shrink_maps.py ../aneurysm-knockout/s1/param-x/maps-s1b fig-param-s1b-s
                                  (maps of the article; the prefix "param_" of plot_maps_s1.py is dropped)
"""
import glob
import os
import sys

from PIL import Image

SRC, DST = (sys.argv[1], sys.argv[2]) if len(sys.argv) == 3 else ("fig-param", "fig-param-s")
SCALE = 1.0 / 3.0
os.makedirs(DST, exist_ok=True)
tot_in = tot_out = 0
for p in sorted(glob.glob(os.path.join(SRC, "*.png"))):
    name = os.path.basename(p)
    q = os.path.join(DST, name[len("param_"):] if name.startswith("param_") else name)
    im = Image.open(p)
    w, h = im.size
    im = im.convert("RGB").resize((round(w * SCALE), round(h * SCALE)), Image.LANCZOS)
    im.quantize(colors=256, method=Image.MEDIANCUT, dither=Image.NONE).save(q, optimize=True)
    tot_in += os.path.getsize(p)
    tot_out += os.path.getsize(q)
print(f"{DST}/: {tot_in / 1e6:.0f} MB -> {tot_out / 1e6:.1f} MB")
