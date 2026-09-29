#!/usr/bin/env python3
"""
Automated matrix figure of the MBE-KG aneurysm field maps.

Renders the DEFORMED wall of each case, laid out as a 2x2 matrix:

        elastin only      elastin + 10% collagen
 fus    matrix-fus-el     matrix-fus-elcc
 sac    matrix-sac-el     matrix-sac-elcc

coloured by relative volume J (matrix_J.png) and, separately, by strain-energy
density per reference volume W_R (matrix_energy.png).

QUADRATIC ELEMENTS: the mesh is hex27 (tri-quadratic). Instead of collapsing
each element to its 8 corner nodes (which discards the curved edges/faces and
under-represents the bulge), every hex27 is tessellated into 8 linear sub-hexes
using ALL 27 nodes arranged on their natural 3x3x3 lattice. This reproduces the
quadratic geometry exactly at the node positions and renders smoothly, without
relying on VTK's quadratic-cell node-ordering conventions.

FIELD INTERPOLATION: FEBio's element_data logs ONE value of J (and sed) per
element, so the underlying field is piecewise-constant per element -- a true
quadratic (per-integration-point) field is not available from these logs. For
display we assign that value to the sub-hexes and then nodal-average
(cell_data_to_point_data) so the colour varies smoothly; set SMOOTH=False to
see the raw per-element (flat) values instead.

Inputs (already produced by FEBio, no FEBio Studio needed):
  * mesh connectivity  -> <case>.feb              (hex27, 27 nodes/elem)
  * deformed node xyz  -> <case>-NODES.txt        (logfile node_data x;y;z)
  * per-element J,sed  -> <case>-ALL.txt          (logfile element_data J;sed)

Requires:  pyvista (pip install pyvista);  headless:  export PYVISTA_OFF_SCREEN=true
Run:  python3 render_matrix.py      Out:  matrix_J.png , matrix_energy.png
"""
import numpy as np
import pyvista as pv

# ----------------------------------------------------------------------------
# (grid-row, grid-col, tag, panel-letter)
CASES = [
    (0, 0, "matrix-fus-el",   "(a)"),
    (0, 2, "matrix-fus-elcc", "(b)"),
    (1, 0, "matrix-sac-el",   "(c)"),
    (1, 2, "matrix-sac-elcc", "(d)"),
]
ROW_LABELS = ["Fusiform\n(20% elastin degradation)", "Saccular\n(65% elastin degradation)"]
COL_HEADERS = ["No loss of\ncollagen cross-linking", "10% loss of\ncollagen cross-linking"]
NROW = 2
WINDOW = (1650, 1250)      # ~square data panels (a long tube fills the diagonal)
SMOOTH = True              # nodal-average the per-element field for smooth contours
NCOLORS = 10               # discrete colour bands (house style: stepped colour bar)
CMAP_J = "coolwarm"        # relative volume: blue -> white -> red (diverging)
CMAP_W = "jet"             # strain energy: blue -> red (rainbow)
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
# oblique 3/4 view; the axis is then rolled onto the panel diagonal to fill it
FOCAL = (0.0, 0.0, 7.5)                    # centre of the l=15 mm segment
VIEWDIR = (0.82, 0.50, 0.28)               # camera offset dir from focal
ROLL = -42                                 # deg: axis on diagonal (BL->TR), frees top-left for text
ZOOM = 1.5                                 # a bit out, so the whole segment is visible
# ----------------------------------------------------------------------------

# FEBio hex27 local node index (0-based) -> its (i,j,k) position on the 3x3x3
# lattice, with i,j,k in {0,1,2}. Corners, then 12 edge mid-nodes, then 6 face
# centres, then the body centre.
LATTICE = {
    0: (0, 0, 0), 1: (2, 0, 0), 2: (2, 2, 0), 3: (0, 2, 0),
    4: (0, 0, 2), 5: (2, 0, 2), 6: (2, 2, 2), 7: (0, 2, 2),
    8: (1, 0, 0), 9: (2, 1, 0), 10: (1, 2, 0), 11: (0, 1, 0),
    12: (1, 0, 2), 13: (2, 1, 2), 14: (1, 2, 2), 15: (0, 1, 2),
    16: (0, 0, 1), 17: (2, 0, 1), 18: (2, 2, 1), 19: (0, 2, 1),
    20: (1, 0, 1), 21: (2, 1, 1), 22: (1, 2, 1), 23: (0, 1, 1),
    24: (1, 1, 0), 25: (1, 1, 2), 26: (1, 1, 1),
}
IJK2LOCAL = {ijk: k for k, ijk in LATTICE.items()}      # inverse map
# the 8 sub-cells of the 2x2x2 lattice, each as VTK-hex-ordered lattice corners
_D = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
      (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
SUBHEX = []
for bk in (0, 1):
    for bj in (0, 1):
        for bi in (0, 1):
            SUBHEX.append([(bi + dx, bj + dy, bk + dz) for dx, dy, dz in _D])


def parse_connectivity(feb):
    """Return (n_elem, 27) int array of 0-based node indices."""
    conn, in_el = [], False
    for line in open(feb):
        s = line.strip()
        if s.startswith("<Elements"):
            in_el = True; continue
        if s.startswith("</Elements>"):
            break
        if in_el and s.startswith("<elem"):
            ints = s.split(">", 1)[1].split("<")[0]
            conn.append([int(v) - 1 for v in ints.split(",")])
    return np.asarray(conn, dtype=np.int64)


def last_block(fn):
    blk, cur = [], []
    for line in open(fn):
        s = line.strip()
        if s.startswith("*Time"):
            if cur: blk.append(cur)
            cur = []
        elif s and not s.startswith("*"):
            cur.append(s)
    if cur: blk.append(cur)
    return blk[-1]


def deformed_points(nodes_txt):
    rows = [r.split(",") for r in last_block(nodes_txt)]
    n = max(int(p[0]) for p in rows)
    pts = np.zeros((n, 3))
    for p in rows:
        pts[int(p[0]) - 1] = (float(p[1]), float(p[2]), float(p[3]))
    return pts


def elem_fields(all_txt, n_elem):
    J = np.zeros(n_elem); sed = np.zeros(n_elem)
    for row in last_block(all_txt):
        p = row.split(",")
        eid = int(p[0]) - 1
        J[eid] = float(p[1]); sed[eid] = float(p[2])
    return J, sed


def build_grid(conn, pts, J, sed):
    """Tessellate each hex27 into 8 linear sub-hexes over the 3x3x3 lattice."""
    n_elem = conn.shape[0]
    cells, Jc, Wr = [], [], []
    for e in range(n_elem):
        loc = conn[e]                                   # 27 global node ids
        for sub in SUBHEX:
            ids = [loc[IJK2LOCAL[ijk]] for ijk in sub]  # 8 global node ids
            cells.append([8] + ids)
            Jc.append(J[e]); Wr.append(sed[e])
    cells = np.asarray(cells, np.int64).ravel()
    ctypes = np.full(len(Jc), pv.CellType.HEXAHEDRON, np.uint8)
    grid = pv.UnstructuredGrid(cells, ctypes, pts)
    grid.cell_data["J"] = np.asarray(Jc)
    grid.cell_data["W_R"] = np.asarray(Wr)             # per reference volume
    if SMOOTH:
        grid = grid.cell_data_to_point_data()          # nodal-average -> smooth
    return grid


def ref_node_radii(feb):
    """Reference (undeformed) cylindrical radius of every node, by 1-based id."""
    rad, inn = {}, False
    for line in open(feb):
        s = line.strip()
        if s.startswith("<Nodes"):
            inn = True; continue
        if s.startswith("</Nodes>"):
            break
        if inn and s.startswith("<node"):
            nid = int(s.split('id="')[1].split('"')[0])
            x, y, _ = [float(v) for v in s.split(">", 1)[1].split("<")[0].split(",")]
            rad[nid] = np.hypot(x, y)
    return rad


def outer_diameter_range(feb, pts):
    """(min, max) deformed OUTER-wall diameter in mm (physical dimension)."""
    rad = ref_node_radii(feb)
    ids = np.array(sorted(rad))
    Rr = np.array([rad[i] for i in ids])
    outer = Rr > Rr.max() - 0.01                     # outer-wall nodes
    Rc = np.hypot(pts[ids - 1, 0], pts[ids - 1, 1])
    do = 2.0 * Rc[outer]
    return float(do.min()), float(do.max())


def add_margin_labels(png):
    """Add white margins with rotated row labels (left) and column headers (top)."""
    from PIL import Image, ImageDraw, ImageFont
    img = Image.open(png).convert("RGB")
    W, H = img.size
    LM, TM = 120, 74
    canvas = Image.new("RGB", (W + LM, H + TM), "white")
    canvas.paste(img, (LM, TM))
    draw = ImageDraw.Draw(canvas)
    hfont = ImageFont.truetype(FONT, 30)
    rfont = ImageFont.truetype(FONT, 28)

    def centered(cx, cy, text, font):
        b = draw.multiline_textbbox((0, 0), text, font=font, align="center")
        draw.multiline_text((cx - (b[2] - b[0]) / 2 - b[0], cy - (b[3] - b[1]) / 2 - b[1]),
                            text, fill="black", font=font, align="center")

    # column headers (data columns are centred at ~0.207 W and ~0.793 W)
    centered(LM + 0.2066 * W, TM / 2, COL_HEADERS[0], hfont)
    centered(LM + 0.7934 * W, TM / 2, COL_HEADERS[1], hfont)
    # rotated row labels in the left margin
    for cyf, text in zip((0.25, 0.75), ROW_LABELS):
        b = draw.multiline_textbbox((0, 0), text, font=rfont, align="center")
        w, h = int(b[2] - b[0]), int(b[3] - b[1])
        tile = Image.new("RGBA", (w + 8, h + 8), (255, 255, 255, 0))
        ImageDraw.Draw(tile).multiline_text((4 - b[0], 4 - b[1]), text,
                                            fill="black", font=rfont, align="center")
        rot = tile.rotate(90, expand=True)
        canvas.paste(rot, (int(LM / 2 - rot.width / 2), int(TM + cyf * H - rot.height / 2)), rot)
    canvas.save(png)


APEX_NODE = 102        # outer-wall node at the apex (above inner node 20)
POINTS = {}


def draw_apex_arrows(png, targets):
    """A thin arrow from the empty lower-right of each panel to the apex, labelled 'apex'.
    The tube runs along the panel diagonal, so the lower-right is free and the arrow does
    not cross the wall."""
    from PIL import Image, ImageDraw, ImageFont
    im = Image.open(png).convert("RGB")
    d = ImageDraw.Draw(im)
    font = ImageFont.truetype(FONT, 24)
    for tx, ty in targets:
        sx, sy = tx + 170, ty + 120                  # tail: down and to the right
        ang = np.arctan2(ty - sy, tx - sx)
        L, W = 20, 7                                 # filled head: length, half-width
        bx, by = tx - L * np.cos(ang), ty - L * np.sin(ang)
        d.line([(sx, sy), (bx, by)], fill="black", width=2)
        nx, ny = -np.sin(ang), np.cos(ang)
        d.polygon([(tx, ty), (bx + W * nx, by + W * ny), (bx - W * nx, by - W * ny)], fill="black")
        d.text((sx + 6, sy - 4), "apex", fill="black", font=font)
    im.save(png)



def main():
    grids, stats, diam = {}, {}, {}
    for _, _, tag, *_ in CASES:
        conn = parse_connectivity(f"{tag}.feb")
        pts = deformed_points(f"{tag}-NODES.txt")
        POINTS[tag] = pts
        J, sed = elem_fields(f"{tag}-ALL.txt", conn.shape[0])
        grids[tag] = build_grid(conn, pts, J, sed)
        stats[tag] = {"J": (J.min(), J.max()), "W_R": (sed.min(), sed.max())}
        diam[tag] = outer_diameter_range(f"{tag}.feb", pts)

    campos = [tuple(np.add(FOCAL, VIEWDIR)), FOCAL, (0.0, 0.0, 1.0)]

    def render(field, cmap, clim, fmt, nlab, fname, label, apex_tags=()):
        name, unit = label.split()[0], label.split()[1]     # e.g. "J", "[-]"
        # data panels in columns 0 and 2; discrete colour bar in the merged col 1
        pl = pv.Plotter(shape=(NROW, 3), window_size=WINDOW,
                        col_weights=[1.0, 0.42, 1.0], groups=[(np.s_[:], 1)],
                        off_screen=True, border=False)
        for r, c, tag, letter in CASES:
            pl.subplot(r, c)
            lo, hi = stats[tag][field]
            dlo, dhi = diam[tag]
            pl.add_mesh(grids[tag], scalars=field, cmap=cmap, clim=clim,
                        n_colors=NCOLORS, show_edges=True, edge_color="black",
                        line_width=1, lighting=False, show_scalar_bar=False)
            rng = (f"{name} = {lo:.2f} – {hi:.2f} {unit}\n"
                   f"Ø = {dlo:.2f} – {dhi:.2f} mm")
            pl.add_text(rng, font_size=11, color="black", position="upper_left")
            pl.add_text(letter, font_size=15, color="black", position="lower_right")
            pl.camera_position = campos
            pl.reset_camera()
            pl.camera.Roll(ROLL)            # roll the axis onto the panel diagonal
            pl.camera.zoom(ZOOM)
        # discrete vertical colour bar, centred in the merged middle column
        pl.subplot(0, 1)
        pl.add_mesh(list(grids.values())[0], scalars=field, cmap=cmap, clim=clim,
                    n_colors=NCOLORS, opacity=0.0, lighting=False,
                    show_scalar_bar=True,
                    scalar_bar_args=dict(title=label, vertical=True,
                                         position_x=0.18, position_y=0.14,
                                         width=0.5, height=0.72, n_labels=nlab,
                                         fmt=fmt, color="black",
                                         title_font_size=24, label_font_size=22))
        # where the numbers are read: the apex (insult centre on the global +Y axis,
        # mid-length), on the outer surface above inner node 20 (outer node 102). Its screen
        # position is taken from each panel's camera and an arrow is drawn to it afterwards.
        pl.render()
        targets = []
        for r, c, tag, _ in CASES:
            if tag not in apex_tags:
                continue
            pl.subplot(r, c)
            ren = pl.renderer
            p = POINTS[tag][APEX_NODE - 1]
            ren.SetWorldPoint(float(p[0]), float(p[1]), float(p[2]), 1.0)
            ren.WorldToDisplay()
            x, y, _ = ren.GetDisplayPoint()
            targets.append((x, WINDOW[1] - y))          # display origin is bottom-left
        pl.screenshot(fname)
        draw_apex_arrows(fname, targets)
        add_margin_labels(fname)                         # row/column titles (PIL)
        from PIL import Image                             # raster PDF, like the other figs
        pdf = fname.replace(".png", ".pdf")
        Image.open(fname).convert("RGB").save(pdf, "PDF", resolution=200.0)
        print("wrote", fname, "and", pdf)

    # the apex is marked once per morphology, in Fig. 4 (a) and (c) only
    render("J",   CMAP_J, (0.6, 2.6), "%.1f", 6, "matrix_J.png",      "J [-]",
           apex_tags=("matrix-fus-el", "matrix-sac-el"))
    render("W_R", CMAP_W, (30, 150),  "%.0f", 7, "matrix_energy.png", "W_R [kPa]")


if __name__ == "__main__":
    main()
