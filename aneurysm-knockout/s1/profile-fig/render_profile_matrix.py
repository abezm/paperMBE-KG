#!/usr/bin/env python3
"""
Profile matrices (2026-10-04, author: "matriz de verdad"): fixed d_e = d_c (20 % fusiform, 37 % saccular),
rows = morphology, columns = one profile value each. Copy of render_profile_lam.py (4 x 4 draft).

  sz:  columns s_z = 1.2, 1.8, 2.4, 3.0 mm at nu = 5
  nu:  columns nu = 1, 3, 5, 7 at s_z = 3 mm
  lad:  Fig. 6 joined, two rows (kept by the author); fus, sac:  Fig. 6 split in two (author, 2026-10-04): one row, columns d_e = d_c of Fig. 6 (nu = 5, s_z = 3 mm)

All four figures share ONE colour bar (author, 2026-10-04): limits = smallest and largest NODAL value over
every panel of the four (all_range()), bands 0.1 wide from 0.6 to 2.6.

The colour is the nodal lambda_r/lambda_ro, not an element average (author, 2026-10-05: the average made the
quadratic hex27 look like first-order elements, and the bar stopped at 2.49 under a 2.57 label). Each hex27 is
drawn as 8 linear sub-hexes over its 27 nodes, so the colour varies inside the element.

Column headers on top, row labels (morphology, d_e = d_c) on the left, one colour bar centred on the
two rows (60 % of their height, as in Fig. 4).

  PYTHONDONTWRITEBYTECODE=1 PYVISTA_OFF_SCREEN=true python3 render_profile_matrix.py sz|nu|fus|sac
Out: profile-{sz,nu}-lam-s1.{png,pdf}, ladder-{fus,sac}-lam-s1.{png,pdf}.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "severity-fig"))
import render_matrix_iso as m          # noqa: E402  (imported once, functions only)
import pyvista as pv                   # noqa: E402
import matplotlib                      # noqa: E402
from PIL import Image, ImageDraw, ImageFont   # noqa: E402

# (row, col, panel key, column header, letter); STEM: panel key -> deck stem, relative to this folder
REF = {"fus": "../severity-fig/lad-fus-d20p000-v5", "sac": "../severity-fig/lad-sac-d37p000-lev2"}
LADDER = {"fus": [(13, "lad-fus-d13p000-lev2"), (17, "lad-fus-d17p000-lev"), (19, "lad-fus-d19p000-v5"),
                  (20, "lad-fus-d20p000-v5")],
          "sac": [(22, "lad-sac-d22p000-lev2"), (27, "lad-sac-d27p000-lev2"), (32, "lad-sac-d32p000-lev2"),
                  (37, "lad-sac-d37p000-lev2")]}


def cases_of(var):
    """(CASES, STEM, column headers, row labels) of one figure."""
    cases, stem, L = [], {}, iter("abcdefgh")
    if var == "lad":                                      # Fig. 6 joined (author: "dejar la fig 6 pegada")
        for r, mo in enumerate(("fus", "sac")):
            for c, (d, deck) in enumerate(LADDER[mo]):
                key = f"lad-{mo}-d{d}"
                stem[key] = "../severity-fig/" + deck
                cases.append((r, c, key, rf"$d_e = d_c$ = {d}%", f"({next(L)})"))
        return cases, stem, [], [("Fusiform",), ("Saccular",)]
    if var in LADDER:
        for c, (d, deck) in enumerate(LADDER[var]):
            key = f"{var}-d{d}"
            stem[key] = "../severity-fig/" + deck
            cases.append((0, c, key, rf"$d_e = d_c$ = {d}%", f"({next(L)})"))
        return cases, stem, [h for *_, h, _l in cases], [None]
    values = {"sz": (1.2, 1.8, 2.4, 3.0), "nu": (1, 3, 5, 7)}[var]   # equal steps (author, 2026-10-04)
    for r, mo in enumerate(("fus", "sac")):
        for c, v in enumerate(values):
            nu, sz = (5, v) if var == "sz" else (v, 3.0)
            key = f"{mo}-{var}-{v:g}"
            stem[key] = REF[mo] if (nu, sz) == (5, 3.0) else f"prof-{mo}-nu{nu}-sz{sz:.1f}".replace(".", "p")
            head = rf"$s_z$ = {sz:.1f} mm" if var == "sz" else rf"$\nu$ = {nu}"
            cases.append((r, c, key, head, f"({next(L)})"))
    return (cases, stem, [h for r_, c_, k_, h, l_ in cases if r_ == 0],
            [("Fusiform", "$d_e = d_c$ = 20%"), ("Saccular", "$d_e = d_c$ = 37%")])


VAR = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] in ("sz", "nu", "fus", "sac", "lad") else "sz"
CASES, STEM, COL_HEADERS, ROW_LABELS = cases_of(VAR)
# per-row colour bar: (lo, hi, n bands, labelled ticks, tick format, colour map, below colour)
# Each row: one complete ColorBrewer diverging map, white (position 0.5) at J = 1, which is
# a band edge. Band k on the side of 1 (k = 0 next to 1) sits at map position
# 0.5 +/- (A + (B - A) * ((k + 0.5) / n_side)**G) / 2, so no band is white and consecutive
# bands have nearly equal CIE Lab differences (checked in main). Fusiform J < 0.5 (tube
# ends) has its own band (0 - 0.5). Each bar is drawn from the smallest to the largest J of
# its row, so it shows only colours and values that occur on the tubes.
# A small: the bands next to J = 1 are nearly white, as in Fig. 4; colour grows away from 1.
DIV = {0: ("PuOr_r", 0.06, 0.97, 0.80),  # purple below 1, orange above
       1: ("BrBG", 0.06, 0.97, 0.80),    # brown below 1, teal above
       2: ("PiYG_r", 0.06, 0.97, 0.80)}  # green below 1, magenta above (profile figures, author 2026-10-04)
DE_MIN = {0: 6.0, 1: 8.0}                # min Delta E76 between consecutive bands


def row_cmap(r, below, n_up, rep_up=1):
    """ListedColormap of the row. below: distance from white (0-1) of each band under J = 1,
    low to high J; n_up bands above 1 (each repeated rep_up times, for a finer uniform grid)."""
    name, a, b, g = DIV[r]
    f = lambda k, n: a + (b - a) * ((k + 0.5) / n) ** g
    pos = [0.5 - d / 2 for d in below]
    for k in range(n_up):
        pos += [0.5 + f(k, n_up) / 2] * rep_up
    return matplotlib.colors.ListedColormap([matplotlib.colormaps[name](x) for x in pos],
                                            name=name + "_J1")


def lab(rgb):
    """sRGB (0-1) -> CIE Lab (D65)."""
    c = np.asarray(rgb, float)[..., :3]
    c = np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
    xyz = c @ M.T / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > (6 / 29) ** 3, np.cbrt(xyz), xyz / (3 * (6 / 29) ** 2) + 4 / 29)
    return np.stack([116 * f[..., 1] - 16, 500 * (f[..., 0] - f[..., 1]),
                     200 * (f[..., 1] - f[..., 2])], -1)


def delta_e(a, b):
    return float(np.linalg.norm(lab(a) - lab(b)))


# bands under J = 1 are 0.25 wide in both rows; the one next to 1 is nearly white
_CM = {0: row_cmap(0, [0.95, 0.90, 0.50, 0.14], 11, rep_up=2),   # grid of 0.25: 4 + 22 = 26
       1: row_cmap(1, [0.50, 0.14], 8)}
# (lo, hi, n bands, interior tick step, tick format, colour map, unused)
ROW_BARS = [(0.0, 6.5, 26, 1.0, "%.1f", _CM[0], None),
            (0.5, 3.0, 10, 0.5, "%.1f", _CM[1], None)]
# radius ratio, one bar for both rows: bands 0.1 wide, white at 1 (a band edge), each side on its own
# half of the diverging map. Nodal values span 0.607 - 2.573 (sac nu = 7), so the bands run from 0.6 to
# 2.6: 4 below 1 and 16 above
# PiYG_r in every figure, Fig. 6 included (author, 2026-10-05: same colours in Figs. 6, 7 and 8)
LAM_BAR = (0.6, 2.6, 20, 0.2, "%.1f", row_cmap(2, [0.85, 0.62, 0.40, 0.16], 16), None)
ROW_BARS = [LAM_BAR] * 2
LABEL_CLEAR = 40                         # target clearance (px) label block - mesh, per row
EDGE_W = 1.0                             # edge line width (px)
LIFT = 0.01                              # mm: edge lines moved toward the camera (no z-fight)
HEX27_EDGES = [(0, 8, 1), (1, 9, 2), (2, 10, 3), (3, 11, 0), (4, 12, 5), (5, 13, 6),
               (6, 14, 7), (7, 15, 4), (0, 16, 4), (1, 17, 5), (2, 18, 6), (3, 19, 7)]
GAP = 80                                 # white gap between the two rows (px, ~7% of H)
NROW, NCOL = (1 if VAR in LADDER else 2), 4
PANEL = (418, 531)                       # data panel size (px)
WINDOW = (NCOL * PANEL[0], NROW * PANEL[1])
LM, BARW = (40 if VAR in LADDER else 120), 150                       # left margin (row labels), right bar column
MARGIN = 0.02                            # min free fraction at each panel border
OUT = os.path.join(HERE, "ladder-lam-s1-common.png" if VAR == "lad" else
                   f"ladder-{VAR}-lam-s1.png" if VAR in LADDER else f"profile-{VAR}-lam-s1.png")


def node_lam(nodes_txt):
    """Per node (0-based, as conn): r / r_o, radial position at the end of the insult over the value of
    the block *Time = 1, the nodal ratio of lam_ratio_range()."""
    blocks, cur = {}, None
    for line in open(nodes_txt):
        s = line.strip()
        if s.startswith("*Time"):
            cur = blocks.setdefault(float(s.split("=")[1]), {})
        elif s and not s.startswith("*"):
            q = s.split(",")
            cur[int(q[0])] = np.hypot(float(q[1]), float(q[2]))
    hom, fin = blocks[1.0], blocks[max(blocks)]
    r = np.zeros(max(hom))
    for i in hom:
        r[i - 1] = fin[i] / hom[i]
    return r


def all_range():
    """Smallest and largest nodal lambda_r/lambda_ro over every panel of the four figures."""
    lo, hi = np.inf, -np.inf
    for var in ("sz", "nu", "fus", "sac"):
        cases, stem, *_ = cases_of(var)
        for *_, key, _h, _l in cases:
            p = os.path.join(HERE, stem[key])
            L = m.lam_ratio_range(p + "-NODES.txt")
            lo, hi = min(lo, L[0]), max(hi, L[1])
    return lo, hi


def load():
    data = {}
    for r, c, tag, d, letter in CASES:
        p = os.path.join(HERE, STEM[tag])
        conn = m.parse_connectivity(p + ".feb")
        pts = m.deformed_points(p + "-NODES.txt")
        J, sed = m.elem_fields(p + "-ALL.txt", conn.shape[0])
        L = node_lam(p + "-NODES.txt")                    # nodal, coloured as point data
        assert len(L) == len(pts), (tag, len(L), len(pts))
        # element edges only (hex27 boundaries, 3-node polylines), lifted toward the camera
        edges = {}
        for e in conn:
            for a_, mid, b_ in HEX27_EDGES:
                key = (min(e[a_], e[b_]), max(e[a_], e[b_]))
                edges.setdefault(key, (e[a_], e[mid], e[b_]))
        lines = pv.PolyData()
        vd = np.asarray(m.VIEWDIR, float) / np.linalg.norm(m.VIEWDIR)
        lines.points = pts + LIFT * vd
        lines.lines = np.hstack([[3, *ids] for ids in edges.values()]).astype(np.int64)
        assert lines.n_verts == 0 and lines.n_lines == len(edges)
        data_lines = lines
        grid = m.build_grid(conn, pts, L[conn].mean(axis=1), sed)   # geometry: 8 sub-hexes per hex27
        del grid.cell_data["J"]                           # element value replaced by the nodal one
        grid.point_data["J"] = L
        used = np.unique(conn)
        data[tag] = dict(grid=grid, J=(J.min(), J.max()), L=(L[used].min(), L[used].max()),
                         Lall=L,
                         lines=data_lines,
                         conn=conn, pts=pts,
                         lam=m.lam_ratio_range(p + "-NODES.txt"))
    return data


def render(data, scale, zoom, letters, fname, probes=None):
    """probes = [(key, tag, world point)]: returns also their pixel (x, y) in the screenshot."""
    campos = [tuple(np.add(m.FOCAL, m.VIEWDIR)), m.FOCAL, (0.0, 0.0, 1.0)]
    pl = pv.Plotter(shape=(NROW, NCOL), window_size=WINDOW, off_screen=True, border=False)
    reset_scales = {}
    for r, c, tag, d, letter in CASES:
        lo, hi = ROW_BARS[r][:2]
        pl.subplot(r, c)
        pl.add_mesh(data[tag]["grid"], scalars="J", cmap=ROW_BARS[r][5], clim=(lo, hi),
                    n_colors=ROW_BARS[r][2], show_edges=False,
                    lighting=False, show_scalar_bar=False)          # J < lo -> lowest band
        pl.add_mesh(data[tag]["lines"], color="black", line_width=EDGE_W, lighting=False)
        if letters:
            pl.add_text(letter, font_size=15, color="black", position="lower_right")
        pl.camera_position = campos
        pl.camera.parallel_projection = True
        pl.reset_camera()
        pl.camera.Roll(m.ROLL)
        reset_scales[tag] = pl.camera.parallel_scale
        if scale is not None:
            pl.camera.parallel_scale = scale      # same mm-per-pixel in every panel
        pl.camera.zoom(zoom)
    pl.render()
    pix = {}
    rc = {tag: (r, c) for r, c, tag, *_ in CASES}
    for key, tag, pt in (probes or []):
        pl.subplot(*rc[tag])
        ren = pl.renderer
        ren.SetWorldPoint(*map(float, pt), 1.0)
        ren.WorldToDisplay()
        x, y, _ = ren.GetDisplayPoint()
        pix[key] = (x, WINDOW[1] - y)                     # display origin is bottom-left
    pl.screenshot(fname)
    pl.close()
    return (reset_scales, pix) if probes else reset_scales


def panel_extents(png):
    """Per panel: (left, right, top, bottom) free pixels between mesh and panel border."""
    a = np.asarray(Image.open(png).convert("RGB")).astype(int)
    out = {}
    for r, c, tag, *_ in CASES:
        x0, y0 = c * PANEL[0], r * PANEL[1]
        sub = a[y0:y0 + PANEL[1], x0:x0 + PANEL[0]]
        ys, xs = np.nonzero((sub != 255).any(axis=2))
        out[tag] = (xs.min(), PANEL[0] - 1 - xs.max(), ys.min(), PANEL[1] - 1 - ys.max())
    return out


STRIPS = [0] * NROW                      # extra strip above each row's tube area (set in main)


def row_y(r):
    """Top of row r (strip + tube area) on the final canvas; rows separated by GAP."""
    return sum(STRIPS[k] + PANEL[1] + GAP for k in range(r))


def label_block(tag, d, data):
    """The three mathtext lines stacked into one RGBA image (4 px between lines)."""
    lines = label_lines(tag, d, data)                     # Fig. 6: the d_e = d_c line in bold (author 2026-10-04)
    ims = [m.mathtext_image(s_, m.LABEL_PX, bold=(VAR == "lad" and i == 0)) for i, s_ in enumerate(lines)]
    blk = Image.new("RGBA", (max(t.width for t in ims),
                             sum(t.height for t in ims) + 4 * (len(ims) - 1)), (255, 255, 255, 0))
    yy = 0
    for t in ims:
        blk.alpha_composite(t, (0, yy))                   # keeps glyph alpha (paste would thin it)
        yy += t.height + 4
    return blk


def label_lines(tag, d, data):
    jl, jh = data[tag]["J"]
    ll, lh = data[tag]["lam"]
    return ([d] if VAR == "lad" else []) + [f"$J$ = {jl:.2f} – {jh:.2f}",
            rf"$\lambda_r/\lambda_{{ro}}$ = {ll:.2f} – {lh:.2f}"]


def draw_bar(canvas, r, lo, hi, n, step, fmt, cmap, below, vmin, vmax):
    """Discrete strip of row r, drawn from vmin to vmax (the smallest and largest J on the
    tubes of the row): only the bands that occur, the end ones cut at vmin / vmax. Both ends
    labelled with their value; interior ticks every `step`, dropped when too close to an end.
    Returns (x centre, function J -> y pixel)."""
    draw = ImageDraw.Draw(canvas)
    cols = (cmap(np.linspace(0, 1, n))[:, :3] * 255).round().astype(int)   # as pyvista LUT
    x0 = LM + WINDOW[0] + 14
    sw, sh = 46, int(0.80 * PANEL[1])
    title = m.mathtext_image(r"$\lambda_r/\lambda_{ro}$", m.LABEL_PX)
    T, B = row_y(0) + STRIPS[0], row_y(NROW - 1) + STRIPS[NROW - 1] + PANEL[1]   # tube area
    sh = int((0.70 if NROW == 1 else 0.60) * (B - T))     # centred, 60 % of the two rows (Fig. 4)
    ytop = int((T + B) / 2 - sh / 2)
    yof = lambda v: ytop + sh * (vmax - v) / (vmax - vmin)
    w = (hi - lo) / n
    for i in range(n):
        b0, b1 = max(lo + i * w, vmin), min(lo + (i + 1) * w, vmax)
        if b1 > b0:
            draw.rectangle([x0, int(round(yof(b1))), x0 + sw, int(round(yof(b0)))], fill=tuple(cols[i]))
    # s1/ (3 Oct 2026): an end band thinner than 3 px (vmin = 0.248 sits 0.002 below the band edge 0.25) is
    # covered by its neighbour; it is drawn again over the last 3 px of that end, so that every colour on the
    # tubes also appears in the bar
    occ = [i for i in range(n) if min(lo + (i + 1) * w, vmax) > max(lo + i * w, vmin)]
    for i, end in ((occ[0], vmin), (occ[-1], vmax)):
        b0, b1 = max(lo + i * w, vmin), min(lo + (i + 1) * w, vmax)
        if abs(int(round(yof(b0))) - int(round(yof(b1)))) < 3:
            ye = int(round(yof(end)))
            y0, y1 = (ye - 3, ye) if end == vmin else (ye, ye + 3)
            draw.rectangle([x0, y0, x0 + sw, y1], fill=tuple(cols[i]))
    canvas.paste(title, (int(x0 + sw / 2 - title.width / 2), int(ytop - title.height - 48)), title)   # clear of the end value
    ends = [(vmin, fmt % vmin), (vmax, fmt % vmax)]                 # ends rounded to 0.1 (author, 2026-10-05)
    inner = [(v, fmt % v) for v in np.arange(np.ceil(vmin / step) * step, vmax, step)
             if min(abs(yof(v) - yof(vmin)), abs(yof(v) - yof(vmax))) >= 1.5 * m.LABEL_PX]   # 2.40 hit 2.44
    n_all = len(np.arange(np.ceil(vmin / step) * step, vmax, step))
    print("bar ticks: %d of %d interior kept" % (len(inner), n_all))
    for v, txt in inner:
        t = m.mathtext_image(txt, m.LABEL_PX)
        canvas.paste(t, (x0 + sw + 6, int(round(yof(v) - t.height / 2))), t)
    # the two end values sit just outside the bar ends, so they never crowd the round ticks
    t = m.mathtext_image(ends[1][1], m.LABEL_PX)
    canvas.paste(t, (x0 + sw + 6, int(round(yof(vmax) - t.height + 4))), t)
    t = m.mathtext_image(ends[0][1], m.LABEL_PX)
    canvas.paste(t, (x0 + sw + 6, int(round(yof(vmin) - 4))), t)
    return x0 + sw // 2, yof


def main():
    data = load()
    for r, c, tag, d, _ in CASES:                         # values against the row limits
        lo, hi = ROW_BARS[r][:2]
        jl, jh = data[tag]["L"]
        assert lo <= jl and jh <= hi, (tag, jl, jh, lo, hi)
    _all = all_range()
    # one bar, common to the four figures, its ends rounded outward to 0.1 (author, 2026-10-05: 2.57 -> 2.6)
    _bar = (np.floor(_all[0] * 10 + 1e-9) / 10, np.ceil(_all[1] * 10 - 1e-9) / 10)
    assert _bar[0] <= _all[0] and _all[1] <= _bar[1], (_bar, _all)
    RANGE = {r: _bar for r in range(NROW)}
    print("assert: all J inside the band limits; bar ranges", RANGE)
    blocks = {tag: label_block(tag, d, data) for _, _, tag, d, _ in CASES}

    import tempfile
    tmp = os.path.join(tempfile.gettempdir(), "profile-matrix-calib.png")   # scratch, removed below
    s = render(data, None, 1.0, False, tmp)              # reset scales, per mesh
    S = max(s.values())                                   # most dilated tube's reset scale
    print("reset parallel_scale per panel:", {k: round(v, 3) for k, v in s.items()}, "-> S =", round(S, 3))
    render(data, S, 1.0, False, tmp)                      # common scale, zoom 1
    ext = panel_extents(tmp)
    # zoom shrinks/grows about the panel centre: free = half - Z*(half - free1)
    zmax = []
    for tag, (fl, fr, ft, fb) in ext.items():
        hw, hh = PANEL[0] / 2, PANEL[1] / 2
        for free, half, dim in ((fl, hw, PANEL[0]), (fr, hw, PANEL[0]), (ft, hh, PANEL[1]), (fb, hh, PANEL[1])):
            zmax.append((half - MARGIN * dim - 1) / (half - free))
    zoom = np.floor(min(zmax) * 100) / 100
    render(data, S, zoom, False, tmp)
    ext = panel_extents(tmp)
    for tag, (fl, fr, ft, fb) in ext.items():
        ok = min(fl, fr) >= MARGIN * PANEL[0] and min(ft, fb) >= MARGIN * PANEL[1]
        assert ok, (tag, fl, fr, ft, fb)
    print("zoom = %.2f; margin check passed (min free px L/R/T/B):" % zoom,
          [min(e[i] for e in ext.values()) for i in range(4)])
    mesh_mask = (np.asarray(Image.open(tmp).convert("RGB")).astype(int) != 255).any(axis=2)
    os.remove(tmp)

    # label placement: per row, the lowest common offset p (label top relative to the tube
    # area top) such that every label pixel of the row is >= LABEL_CLEAR px from any mesh
    # pixel, searching downward from fully above the tube; the strip above the row is what
    # the label needs above the tube area (6 px top pad), possibly 0.
    from scipy.ndimage import distance_transform_edt
    PAD = 400
    label_p = []
    for r in range(NROW):
        rowmask = np.vstack([np.zeros((PAD, WINDOW[0]), bool),
                             mesh_mask[r * PANEL[1]:(r + 1) * PANEL[1]]])
        dist = {}                                         # distance to own panel's mesh only
        for c in range(NCOL):
            own = np.zeros_like(rowmask)
            own[:, c * PANEL[0]:(c + 1) * PANEL[0]] = rowmask[:, c * PANEL[0]:(c + 1) * PANEL[0]]
            dist[c] = distance_transform_edt(~own)
        hb = max(blocks[t].height for rr, _, t, *_ in CASES if rr == r)

        def clear_at(p):
            for rr, c, tag, *_ in CASES:
                if rr != r:
                    continue
                a = np.asarray(blocks[tag])[:, :, 3] > 0
                x0, y0 = c * PANEL[0] + 6, PAD + p
                if (dist[c][y0:y0 + a.shape[0], x0:x0 + a.shape[1]][a] < LABEL_CLEAR).any():
                    return False
            return True
        p = -hb - LABEL_CLEAR - 6
        assert clear_at(p)
        while p + 1 <= PANEL[1] - hb and clear_at(p + 1):
            p += 1
        label_p.append(p)
        STRIPS[r] = max(0, 6 - p)
    print("label offset p per row:", label_p, "strip px per row:", STRIPS)

    # pixel checks in panels (d) and (h): the visible outer face-centre node (face towards the camera)
    # with the highest nodal value, and the one with the lowest (contracted tube end, which must show
    # the lowest band); that node is projected to the image and its own value is the expected colour
    vd = np.asarray(m.VIEWDIR) / np.linalg.norm(m.VIEWDIR)
    probes, probe_J = [], {}
    for r, tag in ((r_, t_) for r_, c_, t_, *_ in CASES if c_ == NCOL - 1):   # last panel of each row
        conn, pts, J = data[tag]["conn"], data[tag]["pts"], data[tag]["Lall"]
        vis = []
        for e in range(conn.shape[0]):
            fn = conn[e, 20:26]                           # the six face-centre nodes
            k = fn[np.argmax(np.hypot(pts[fn, 0], pts[fn, 1]))]
            o = pts[k]
            rad = np.array([o[0], o[1], 0.0]) / np.hypot(o[0], o[1])
            if rad @ vd > 0.8:
                vis.append((J[k], o))
        for key, (Je, o) in (("max", max(vis, key=lambda v: v[0])), ("min", min(vis, key=lambda v: v[0]))):
            probes.append(((tag, key), tag, o))
            probe_J[(tag, key)] = (r, Je)
    _, pix = render(data, S, zoom, True, OUT, probes)
    img = Image.open(OUT).convert("RGB")
    canvas = Image.new("RGB", (LM + WINDOW[0] + BARW, row_y(NROW) - GAP), "white")
    mask = np.zeros((canvas.height, WINDOW[0]), bool)
    for r in range(NROW):                                 # rows pasted with a white gap
        y0 = row_y(r) + STRIPS[r]                          # tube area below the strip
        canvas.paste(img.crop((0, r * PANEL[1], WINDOW[0], (r + 1) * PANEL[1])), (LM, y0))
        mask[y0:y0 + PANEL[1]] = mesh_mask[r * PANEL[1]:(r + 1) * PANEL[1]]
    mesh_mask = mask

    # per-panel labels, left-aligned, one common vertical offset per row
    clear = {}
    for r, c, tag, d, _ in CASES:
        own = np.zeros_like(mesh_mask)
        own[:, c * PANEL[0]:(c + 1) * PANEL[0]] = mesh_mask[:, c * PANEL[0]:(c + 1) * PANEL[0]]
        dist = distance_transform_edt(~own)
        x, y = c * PANEL[0] + 6, row_y(r) + STRIPS[r] + label_p[r]
        blk = blocks[tag]
        canvas.paste(blk, (LM + x, y), blk)
        a = np.asarray(blk)[:, :, 3] > 0
        dmin = dist[y:y + a.shape[0], x:x + a.shape[1]][a].min()
        assert dmin >= LABEL_CLEAR - 1, (tag, dmin)        # clearance to own panel's mesh
        clear[tag] = round(float(dmin), 1)
        print(tag, "J = %.2f - %.2f, lam_r/lam_ro nodes %.2f - %.2f, elements %.3f - %.3f"
              % (*data[tag]["J"], *data[tag]["lam"], *data[tag]["L"]))
    for r in range(NROW):                                 # closest label of each row: 35-45 px
        cl = min(v for (rr, _, t, *_), v in zip(CASES, clear.values()) if rr == r)
        assert 35 <= cl <= 45, (r, cl)
    print("assert: label-mesh clearance (own panel) px:", clear)

    # row labels, rotated, in the left margin (style of Fig. 4): morphology, then d_e = d_c
    HEAD_PX = 36                                          # bold, larger than the panel labels (author 2026-10-04)
    for r, lines in enumerate(ROW_LABELS):
        if lines is None:                                 # one-row figure: morphology in the caption
            continue
        ims = [m.mathtext_image(t, HEAD_PX, bold=True) for t in lines]
        tile = Image.new("RGBA", (max(t.width for t in ims), sum(t.height for t in ims) + 6), (255, 255, 255, 0))
        yy = 0
        for t in ims:
            tile.alpha_composite(t, ((tile.width - t.width) // 2, yy))
            yy += t.height + 6
        rot = tile.rotate(90, expand=True)
        canvas.paste(rot, (int(LM / 2 - rot.width / 2), int(row_y(r) + STRIPS[r] + PANEL[1] / 2 - rot.height / 2)), rot)

    bars = [draw_bar(canvas, 0, *ROW_BARS[0], *RANGE[0])] * NROW   # one bar over both rows

    a = np.asarray(canvas).astype(int)
    for key, (r, Jp) in probe_J.items():
        lo, hi, n = ROW_BARS[r][:3]
        x, y = pix[key]
        cx, cy = int(round(LM + x)), int(round(row_y(r) + STRIPS[r] + y - r * PANEL[1]))
        win = a[cy - 3:cy + 4, cx - 3:cx + 4].reshape(-1, 3)   # face colour = most frequent
        vals, cnt = np.unique(win, axis=0, return_counts=True)   # (edge lines are dense here)
        mesh = vals[np.argmax(cnt)]
        band = min(n - 1, max(0, int((Jp - lo) / (hi - lo) * n)))   # clamped like VTK
        bx, yof = bars[r]
        wb = (hi - lo) / n
        jm = (max(lo + band * wb, RANGE[r][0]) + min(lo + (band + 1) * wb, RANGE[r][1])) / 2
        bar = a[int(round(yof(jm))), bx]
        diff = int(np.abs(mesh - bar).max())
        # rABC-7: the saccular bulge can hide the probed face behind elements of a higher band; the pixel
        # must then show the colour of a band of the bar one step above (never a colour absent from it)
        cols = (ROW_BARS[r][5](np.linspace(0, 1, n))[:, :3] * 255).round().astype(int)
        near = int(np.argmin(np.abs(cols - mesh).max(axis=1)))
        dnear = int(np.abs(cols[near] - mesh).max())
        print(f"pixel check {key}: lam = {Jp:.3f} (band {band}), mesh {mesh.astype(int).tolist()}"
              f" at ({cx},{cy}), bar {bar.tolist()}, max diff {diff}; nearest bar band {near} (diff {dnear})")
        # profile figure: sac nu = 6 hides the probed face behind elements two bands higher
        assert dnear <= 3 and (near == band or (key[1] == "max" and near > band)), (key, mesh, bar, near)

    heads = [m.mathtext_image(h, HEAD_PX, bold=True) for h in COL_HEADERS]
    TM = max(t.height for t in heads) + 24 if heads else 0   # lad: no column headers
    full = Image.new("RGB", (canvas.width, canvas.height + TM), "white")
    full.paste(canvas, (0, TM))
    for c, t in enumerate(heads):
        full.paste(t, (int(LM + (c + 0.5) * PANEL[0] - t.width / 2), int((TM - t.height) / 2)), t)
    canvas = full
    canvas.save(OUT)
    pdf = OUT.replace(".png", ".pdf")
    canvas.save(pdf, "PDF", resolution=200.0)
    W, H = canvas.size
    print("wrote", OUT, pdf, "size %d x %d px = 170 x %.1f mm; label em = %.2f pt at 170 mm"
          % (W, H, 170.0 * H / W, m.LABEL_PX * 170.0 / W / 0.3528))


if __name__ == "__main__":
    main()
