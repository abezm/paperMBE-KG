#!/usr/bin/env python3
"""Simulated review, point 16 (2026-10-04, approved by the author): axial mesh refinement, segment length and
end condition of the four lesions of Figure 4 (s1/param-x/matrix-{fus,sac}-{el,elcc}.feb, read only).

Cases written here (each in its own folder, nothing in param-x is touched):
  z2   : the same tube with every axial element split in two (40 axial elements; the hex27 node layers of the
         original are kept and the mid-layers are added, so the cosine grading of the original is kept)
  L30  : the z2 mesh stretched to l_0 = 30 mm (same element lengths as the original), insult centred at 15 mm
  gate : the original mesh rebuilt by this script (20 axial elements); it must reproduce the param-x results
The end-condition case is written by gen_ends.py (it needs the axial reaction of a fixed-end run).

hex27 node ordering: the TEMPLATE of ../../mesh-study/gen_mesh.py, verified there against all 400 elements of
the original mesh. Radial layers (nr = 1) and circumferential layers are taken verbatim from the source deck.
"""
import re
import sys
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
PX = HERE.parent / "param-x"
TEMPLATE = [(0,2,0),(0,0,0),(2,0,0),(2,2,0),(0,2,2),(0,0,2),(2,0,2),(2,2,2),
            (0,1,0),(1,0,0),(2,1,0),(1,2,0),(0,1,2),(1,0,2),(2,1,2),(1,2,2),
            (0,2,1),(0,0,1),(2,0,1),(2,2,1),(0,1,1),(1,0,1),(2,1,1),(1,2,1),
            (1,1,0),(1,1,2),(1,1,1)]
FACE = [(0,2,0),(0,0,0),(0,0,2),(0,2,2),(0,1,0),(0,0,1),(0,1,2),(0,2,1),(0,1,1)]


def layers(t):
    nodes = np.array([[float(v) for v in m.group(1).split(",")]
                      for m in re.finditer(r'<node id="\s*\d+"\s*>\s*([^<]+)</node>', t)])
    r = np.hypot(nodes[:, 0], nodes[:, 1])
    th = np.degrees(np.arctan2(nodes[:, 1], nodes[:, 0])) % 360
    return np.unique(r.round(6)), np.unique(th.round(4)), np.unique(nodes[:, 2].round(9))


def build(rv, tv, zv):
    NR_L, NTH_L, NZ_L = len(rv), len(tv), len(zv)
    nr, nth, nz = (NR_L - 1) // 2, NTH_L // 2, (NZ_L - 1) // 2
    nid, coords = {}, []
    for i in range(NR_L):
        for j in range(NTH_L):
            for k in range(NZ_L):
                nid[(i, j, k)] = len(coords) + 1
                a = np.radians(tv[j])
                coords.append((rv[i] * np.cos(a), rv[i] * np.sin(a), zv[k]))
    elems = [[nid[(2*ei + d[0], (2*ej + d[1]) % NTH_L, 2*ek + d[2])] for d in TEMPLATE]
             for ei in range(nr) for ej in range(nth) for ek in range(nz)]
    ends = [0, NZ_L - 1]
    fixz = [nid[(i, j, k)] for i in range(NR_L) for j in range(NTH_L) for k in ends]
    th_at = lambda deg: int(np.argmin(np.abs(((tv - deg + 180) % 360) - 180)))
    fixx = [nid[(0, th_at(d), k)] for d in (90.0, 270.0) for k in ends]
    fixy = [nid[(0, th_at(d), k)] for d in (0.0, 180.0) for k in ends]
    faces = [[nid[(d[0], (2*ej + d[1]) % NTH_L, 2*ek + d[2])] for d in FACE]
             for ej in range(nth) for ek in range(nz)]
    # end faces (outward normal +z at z = L): quad9 on the k = NZ_L-1 layer
    endL = [[nid[(2*ei + a, (2*ej + b) % NTH_L, NZ_L - 1)] for a, b in
             ((0,0),(2,0),(2,2),(0,2),(1,0),(2,1),(1,2),(0,1),(1,1))]
            for ei in range(nr) for ej in range(nth)]
    inner = [nid[(0, j, k)] for j in range(NTH_L) for k in range(NZ_L)]
    return coords, elems, sorted(fixx), sorted(fixy), sorted(fixz), faces, endL, inner


def emit(src, dst, zv, z_center=None):
    t = src.read_text()
    rv, tv, _ = layers(t)
    coords, elems, fx, fy, fz, faces, endL, inner = build(rv, tv, zv)
    nx = "\n".join(f'\t\t\t<node id="{i+1:6d}"> {c[0]:.15f}, {c[1]:.15f}, {c[2]:.15f} </node>' for i, c in enumerate(coords))
    ex = "\n".join(f'\t\t\t<elem id="{i+1:6d}">' + ",".join(f"{n:6d}" for n in e) + "</elem>" for i, e in enumerate(elems))
    t = re.sub(r'(<Nodes name="[^"]*">).*?(</Nodes>)', lambda m: m.group(1) + "\n" + nx + "\n\t\t" + m.group(2), t, flags=re.S)
    t = re.sub(r'(<Elements type="hex27"[^>]*>).*?(</Elements>)', lambda m: m.group(1) + "\n" + ex + "\n\t\t" + m.group(2), t, flags=re.S)
    ns = lambda name, ids: f'<NodeSet name="{name}">\n' + "\n".join(f"\t\t\t{i}," for i in ids).rstrip(",") + "\n\t\t</NodeSet>"
    for nm, ids in (("FixedDisplacement1", fx), ("FixedDisplacement2", fy), ("FixedDisplacement3", fz)):
        assert len(re.findall(rf'<NodeSet name="{nm}"', t)) == 1, nm
        t = re.sub(rf'<NodeSet name="{nm}"\s*>.*?</NodeSet>', ns(nm, ids), t, flags=re.S)
    fx9 = "\n".join(f'\t\t\t<quad9 id="{i+1:6d}">' + ",".join(f"{n:6d}" for n in f) + "</quad9>" for i, f in enumerate(faces))
    t = re.sub(r'(<Surface name="PressureLoad1"\s*>).*?(</Surface>)', lambda m: m.group(1) + "\n" + fx9 + "\n\t\t" + m.group(2), t, flags=re.S)
    if z_center is not None:
        assert t.count("<z_center>7.5</z_center>") == 1
        t = t.replace("<z_center>7.5</z_center>", f"<z_center>{z_center}</z_center>")
    # outputs: J of every element and the position of every inner-wall node (ids differ between meshes)
    tag = dst.stem
    out = (f'\t\t<element_data file="{tag}-ALL.txt" data="J" name="all_dat" delim="," >'
           + ",".join(str(i) for i in range(1, len(elems) + 1)) + " </element_data>\n"
           f'\t\t<node_data file="{tag}-INNER.txt" data="x;y;z" name="inner" delim="," >'
           + ",".join(str(i) for i in inner) + " </node_data>\n")
    t = re.sub(r"\t*<(element|node)_data[^\n]*\n", "", t)
    t = t.replace("</logfile>", out + "\t\t</logfile>", 1)
    assert t.count(f'file="{tag}-ALL.txt"') == 1
    dst.parent.mkdir(exist_ok=True)
    assert not dst.exists(), dst
    dst.write_text(t)
    # end face of z = L for gen_ends.py
    (dst.parent / f"{tag}.endface").write_text("\n".join(",".join(map(str, f)) for f in endL))
    print(f"{dst.relative_to(HERE)}: {len(coords)} nodes, {len(elems)} elements, z in [{zv[0]}, {zv[-1]}]")


if __name__ == "__main__":
    for case in ("fus-el", "fus-elcc", "sac-el", "sac-elcc"):
        src = PX / f"matrix-{case}.feb"
        _, _, z0 = layers(src.read_text())
        assert len(z0) == 41 and abs(z0[0]) < 1e-12 and abs(z0[-1] - 15.0) < 1e-9
        z2 = np.sort(np.concatenate([z0, 0.5 * (z0[1:] + z0[:-1])]))
        emit(src, HERE / "gate" / f"{case}.feb", z0)
        emit(src, HERE / "z2" / f"{case}.feb", z2)
        emit(src, HERE / "L30" / f"{case}.feb", 2.0 * z2, z_center=15.0)
