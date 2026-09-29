#!/usr/bin/env python3
"""Structured hex27 tube mesh generator for the mesh-convergence study (referee item M8).

The paper has ONE mesh (1 radial x 20 circumferential x 20 axial, hex27) and inherits its
mesh-convergence claim from latorre2020numerical, which is a DIFFERENT (constrained-mixture)
model.  Nobody has ever refined THIS model's mesh.  With a single element through the
thickness a through-thickness ("accordion") mode cannot even be represented, so its absence
is not evidence of stability.

Nothing here is guessed:
  * the hex27 local node ordering is DERIVED from the original mesh (element 1) and then
    VERIFIED against all 400 elements;
  * the geometry, boundary conditions and pressure surface are read off the original deck:
      r_in=0.6468, r_out=0.6870, L=15, full 360 deg tube
      z_dof=0 on both end faces (240 nodes)
      x_dof=0 at theta=90,270 (r_in, both ends);  y_dof=0 at theta=0,180 (r_in, both ends)
      pressure on the INNER face, quad9 (400 faces)

VALIDATION GATE: regenerate the CURRENT mesh (nr=1) and require the MBE run to reproduce the
known apex answer (J=1.02340 at the apex).  Only if that passes are nr=2,3 trusted.
"""
import re
import numpy as np
from pathlib import Path

HERE = Path(__file__).parent
SRC = HERE.parent / "full-KG" / "rate-fusiform-MBE.feb"

R_IN, R_OUT, LZ = 0.6468, 0.6870, 15.0
NTH, NZ = 20, 20          # circumferential / axial elements (unchanged)

# The original axial spacing is NOT uniform: it is cosine-graded, refined at the centre
# (spacing 0.2968 at the defect) and coarse at the ends (0.4532).  A naive linspace does
# NOT reproduce it -- the validation gate caught exactly that.  Since we only refine
# RADIALLY, take the theta and z node layers VERBATIM from the original deck and subdivide
# r alone, so nothing but the radial discretisation changes.

# hex27 local ordering, as offsets (d_ir, d_ith, d_iz) from the element's base corner.
TEMPLATE = [(0,2,0),(0,0,0),(2,0,0),(2,2,0),(0,2,2),(0,0,2),(2,0,2),(2,2,2),
            (0,1,0),(1,0,0),(2,1,0),(1,2,0),(0,1,2),(1,0,2),(2,1,2),(1,2,2),
            (0,2,1),(0,0,1),(2,0,1),(2,2,1),(0,1,1),(1,0,1),(2,1,1),(1,2,1),
            (1,1,0),(1,1,2),(1,1,1)]


def parse_original():
    t = SRC.read_text()
    nodes = {int(m.group(1)): [float(v) for v in m.group(2).split(',')]
             for m in re.finditer(r'<node id="\s*(\d+)"\s*>\s*([^<]+)</node>', t)}
    elems = {int(m.group(1)): [int(v) for v in m.group(2).split(',')]
             for m in re.finditer(r'<elem id="\s*(\d+)"\s*>\s*([^<]+)</elem>', t)}
    return t, nodes, elems


def verify_template():
    """The ordering was read off element 1; require it to hold for ALL 400 elements."""
    _, nodes, elems = parse_original()
    ids = sorted(nodes)
    a = np.array([nodes[i] for i in ids])
    r = np.hypot(a[:, 0], a[:, 1]); th = np.degrees(np.arctan2(a[:, 1], a[:, 0])) % 360; z = a[:, 2]
    ur, uth, uz = np.unique(r.round(6)), np.unique(th.round(3)), np.unique(z.round(6))
    idx = {n: (int(np.searchsorted(ur, round(rr, 6))),
               int(np.searchsorted(uth, round(tt, 3))),
               int(np.searchsorted(uz, round(zz, 6))))
           for n, rr, tt, zz in zip(ids, r, th, z)}
    nth_layers = len(uth)
    bad = 0
    for eid, conn in elems.items():
        base_r = min(idx[n][0] for n in conn)
        base_z = min(idx[n][2] for n in conn)
        ths = sorted({idx[n][1] for n in conn})
        # wrap-around column: theta indices are {38,39,0}, so the base is 38, not 39
        base_t = ths[0] if (max(ths) - min(ths)) <= 2 else min(t for t in ths if t >= nth_layers - 2)
        got = []
        for n in conn:
            dr = idx[n][0] - base_r
            dt = (idx[n][1] - base_t) % nth_layers
            dz = idx[n][2] - base_z
            got.append((dr, dt, dz))
        if got != TEMPLATE:
            bad += 1
    print(f"template check: {len(elems)-bad}/{len(elems)} elementos coinciden con el TEMPLATE"
          f"{'  -> OK' if bad == 0 else f'  -> {bad} FALLAN'}")
    return bad == 0


def build(nr):
    """Structured hex27 tube with nr elements through the thickness."""
    _, onodes, _ = parse_original()
    oa = np.array([onodes[i] for i in sorted(onodes)])
    otv = np.unique((np.degrees(np.arctan2(oa[:, 1], oa[:, 0])) % 360).round(3))
    ozv = np.unique(oa[:, 2].round(6))
    NR_L, NTH_L, NZ_L = 2*nr + 1, len(otv), len(ozv)      # node layers
    rv = np.linspace(R_IN, R_OUT, NR_L)     # r IS uniform in the original (0.0201, 0.0201)
    tv, zv = otv, ozv                        # theta and z taken verbatim

    nid = {}
    coords = []
    for i in range(NR_L):
        for j in range(NTH_L):
            for k in range(NZ_L):
                nid[(i, j, k)] = len(coords) + 1
                a = np.radians(tv[j])
                coords.append((rv[i]*np.cos(a), rv[i]*np.sin(a), zv[k]))

    elems = []
    for ei in range(nr):
        for ej in range(NTH):
            for ek in range(NZ):
                b = (2*ei, 2*ej, 2*ek)
                elems.append([nid[(b[0]+d[0], (b[1]+d[1]) % NTH_L, b[2]+d[2])] for d in TEMPLATE])

    def near(v, x): return abs(v - x) < 1e-9
    ends = [k for k in range(NZ_L) if near(zv[k], 0.0) or near(zv[k], LZ)]
    fixz = [nid[(i, j, k)] for i in range(NR_L) for j in range(NTH_L) for k in ends]
    def th_at(deg):
        return int(np.argmin(np.abs(((tv - deg + 180) % 360) - 180)))
    fixx = [nid[(0, th_at(d), k)] for d in (90.0, 270.0) for k in ends]
    fixy = [nid[(0, th_at(d), k)] for d in (0.0, 180.0) for k in ends]

    # inner-surface quad9 faces (local face = the ir=0 face of each inner element)
    FACE = [(0,2,0),(0,0,0),(0,0,2),(0,2,2),(0,1,0),(0,0,1),(0,1,2),(0,2,1),(0,1,1)]
    faces = []
    for ej in range(NTH):
        for ek in range(NZ):
            b = (0, 2*ej, 2*ek)
            faces.append([nid[(b[0]+d[0], (b[1]+d[1]) % NTH_L, b[2]+d[2])] for d in FACE])

    return coords, elems, sorted(fixx), sorted(fixy), sorted(fixz), faces


def emit(nr, dst):
    t, _, _ = parse_original()
    coords, elems, fx, fy, fz, faces = build(nr)

    nodes_xml = "\n".join(f'\t\t\t<node id="{i+1:6d}"> {c[0]:.15f}, {c[1]:.15f}, {c[2]:.15f} </node>'
                          for i, c in enumerate(coords))
    elems_xml = "\n".join(f'\t\t\t<elem id="{i+1:6d}">' + ",".join(f"{n:6d}" for n in e) + '</elem>'
                          for i, e in enumerate(elems))
    t = re.sub(r'(<Nodes name="[^"]*">).*?(</Nodes>)', lambda m: m.group(1)+"\n"+nodes_xml+"\n\t\t"+m.group(2), t, flags=re.S)
    t = re.sub(r'(<Elements type="hex27"[^>]*>).*?(</Elements>)', lambda m: m.group(1)+"\n"+elems_xml+"\n\t\t"+m.group(2), t, flags=re.S)

    def nodeset(name, ids):
        return f'<NodeSet name="{name}">\n' + "\n".join(f"\t\t\t{i}," for i in ids).rstrip(",") + f"\n\t\t</NodeSet>"
    for nm, ids in (("FixedDisplacement1", fx), ("FixedDisplacement2", fy), ("FixedDisplacement3", fz)):
        t = re.sub(rf'<NodeSet name="{nm}"\s*>.*?</NodeSet>', nodeset(nm, ids), t, flags=re.S)

    faces_xml = "\n".join(f'\t\t\t<quad9 id="{i+1:6d}">' + ",".join(f"{n:6d}" for n in f) + '</quad9>'
                          for i, f in enumerate(faces))
    t = re.sub(r'(<Surface name="PressureLoad1"\s*>).*?(</Surface>)',
               lambda m: m.group(1)+"\n"+faces_xml+"\n\t\t"+m.group(2), t, flags=re.S)

    # Log J for EVERY element, so the max dilatation is comparable ACROSS meshes
    # (element ids differ between meshes, so a fixed id like "element 90" is meaningless).
    # This MUST live inside emit(): patching it onto the file afterwards gets silently wiped
    # the next time the mesh is regenerated -- which already caused one stale-file misread.
    allj = ('\t\t<element_data file="{tag}_ALLJ.txt" data="J" name="allJ" delim="," >'
            .format(tag=Path(dst).stem)
            + ",".join(str(i) for i in range(1, len(elems)+1)) + ' </element_data>\n')
    t = t.replace('\t\t<element_data file=', allj + '\t\t<element_data file=', 1)

    tag = Path(dst).stem
    t = re.sub(r'file="([^"]+)\.txt"', rf'file="{tag}_\1.txt"', t)
    Path(dst).write_text(t)
    print(f"  nr={nr}: {len(coords)} nodos, {len(elems)} elems, fixz={len(fz)}, presion={len(faces)} caras -> {Path(dst).name}")


if __name__ == "__main__":
    assert verify_template(), "TEMPLATE no valida: ABORTAR (no generar nada)"
    print()
    for nr in (1, 2, 3):
        emit(nr, HERE / f"mesh_nr{nr}.feb")
