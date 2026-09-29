import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import matplotlib.ticker as mtick
import matplotlib.colors as mcolors 
import matplotlib.patheffects as pe

def radial_stretch(x,u,y,v):
	return np.sqrt(x**2+y**2)/(np.sqrt((x-u)**2+(y-v)**2))

# Helper function to generate identical formatting for all 6 subplots
J_BAND = 4.0   # upper bound of the admissible states (paper, Section 5.5)
FS = 26      # text size: each panel (two maps) prints ~8 cm wide, so 26 pt here is ~8 pt on paper


def create_contour_plot(ax, x, y, z, lims, cbar_label, levels=15, cmap='viridis_r', use_scatter=False, smooth=False, crange=None, band=None):

	# Establecer los límites ANTES de dibujar
	x_max, y_max, z_max = lims
	ax.set_xlim([0, x_max])
	ax.set_ylim([0, y_max])
	# major ticks from 0 to the limit, the limit included, so both ends carry a label
	_step = 0.1 if x_max <= 0.7 + 1e-9 else 0.2
	ax.set_xticks(np.arange(0, x_max + 1e-9, _step))
	ax.set_yticks(np.arange(0, y_max + 1e-9, _step))

	# Calculate absolute max value for annotation
	z_max_actual = np.max(z)
	# Optional cubic (C1) refinement of the SAME data -> smooth contours (no new data).
	_cargs = (x, y, z)
	if smooth and not use_scatter:
		try:
			from matplotlib.tri import Triangulation, UniformTriRefiner
			_tri, _zc = UniformTriRefiner(Triangulation(x, y)).refine_field(z, subdiv=3)
			_cargs = (_tri, _zc)
		except Exception as _e:
			print(f'  cubic refine failed ({_e}); linear fallback')
	
	# converged states with J > 4, in light grey: the triangles of the computed grid (all
	# converged points, neighbours only) with at least one vertex above J = 4. The coloured
	# fill below covers the triangles whose vertices are all admissible, so the grey band
	# starts where the colour ends; points without a converged solution stay white.
	if band is not None:
		from matplotlib.tri import Triangulation
		bx, by, bJ = band
		tri = Triangulation(bx, by)
		t = tri.triangles
		edge = np.max(np.hypot(bx[t] - bx[np.roll(t, 1, axis=1)], by[t] - by[np.roll(t, 1, axis=1)]), axis=1)
		step = 1.0 / 19.0
		tri.set_mask((edge > 1.01 * np.sqrt(2) * step) | np.all(bJ[t] <= J_BAND, axis=1))
		ax.tricontourf(tri, np.zeros(len(bx)), levels=[-1, 1], colors=['0.82'], zorder=0)

	if use_scatter:
		contour_filled = ax.scatter(x, y, c=z, cmap=cmap, marker='o') 
		
		# Colorbar ticks uniformes desde el min al max de 'levels'

	else:
		# 1. Create a buttery smooth background gradient (150 levels) over the colour range
		# (crange = the colour bar), which is independent of the contour lines drawn
		c0, c1 = crange if crange is not None else (levels[0], levels[-1])
		fill_levels = np.linspace(c0, c1, 150)
			
		contour_filled = ax.tricontourf(*_cargs, levels=fill_levels, 
							cmap=cmap, alpha=0.8)
		
							
	# thin lines, no halo: near the frontier the levels crowd into a narrow strip, and
	# halos + labels there turned into a dark smudge at the printed size
	# few lines, and every line carries its value. Near the frontier the lines crowd (J grows
	# steeply there), so the labels are spread out: start with neighbouring levels in
	# different thirds of their lines
	# (level i in third i mod 3, counted from the d_c axis), then refined: each label in turn
	# moves to the point of its line (away from the frame) farthest from all the other labels
	sets, cands, pos = [], [], []
	for i, lev in enumerate(levels):
		cl = ax.tricontour(*_cargs, levels=[lev], colors='0.15', linewidths=0.5)
		segs = [sg for sg in cl.allsegs[0] if len(sg) > 1]
		if not segs:
			continue
		sg = max(segs, key=lambda q: np.sum(np.hypot(*np.diff(q, axis=0).T)))
		if sg[0, 0] > sg[-1, 0]:
			sg = sg[::-1]                     # start at the d_c-axis end
		arc = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(sg, axis=0).T))])
		# away from the ends and from the frame (half a label from each axis)
		ok = ((arc >= 0.05 * arc[-1]) & (arc <= 0.95 * arc[-1])
		      & (sg[:, 0] > 0.07 * x_max) & (sg[:, 1] > 0.06 * y_max))
		cand = sg[ok] if ok.any() else sg[[len(sg) // 2]]
		f0, f1 = ((0.1, 0.37), (0.37, 0.63), (0.63, 0.9))[i % 3]
		band = np.flatnonzero((arc[ok] >= f0 * arc[-1]) & (arc[ok] <= f1 * arc[-1])) if ok.any() else [0]
		sets.append(cl)
		cands.append(cand / np.array([x_max, y_max]))
		pos.append(band[len(band) // 2] if len(band) else len(cand) // 2)
	# distance in label sizes: the lines run along the diagonal, and a label is ~0.19 of the
	# axis long along its line and ~0.09 high across it (measured on the rendered panels;
	# 0.22 and 0.11 below include a clear margin around each label,
	# halo included); two labels are clear of each other when max(|du|, |dv|) >= 1
	def _lab(q):
		return np.stack([(q[..., 0] - q[..., 1]) / (np.sqrt(2) * 0.22),
		                 (q[..., 0] + q[..., 1]) / (np.sqrt(2) * 0.11)], axis=-1)
	L = [_lab(c) for c in cands]
	def _spread(ps):
		for _ in range(10):
			for n in range(len(ps)):
				others = np.array([L[m][ps[m]] for m in range(len(ps)) if m != n])
				if len(others):
					d = np.min(np.max(np.abs(L[n][:, None, :] - others[None]), axis=2), axis=1)
					ps[n] = int(np.argmax(d))
		q = np.array([L[m][ps[m]] for m in range(len(ps))])
		dm = min((np.max(np.abs(q[m] - q[n])) for m in range(len(ps)) for n in range(m)), default=np.inf)
		return ps, dm
	# several starts (the three-thirds start, its mirror, all at the middle, and 40 random
	# ones with a fixed seed, so the figure is reproducible), keep the best spread
	rng = np.random.default_rng(0)
	starts = [list(pos), [len(c) - 1 - k for c, k in zip(cands, pos)], [len(c) // 2 for c in cands]]
	starts += [[int(rng.integers(len(c))) for c in cands] for _ in range(40)]
	pos, _dm = max((_spread(list(st)) for st in starts), key=lambda r: r[1])
	print(f'  label spacing ({cbar_label}): {_dm:.2f} label sizes')
	texts = []
	for cl, cand, k in zip(sets, cands, pos):
		new = ax.clabel(cl, manual=[tuple(cand[k] * np.array([x_max, y_max]))],
		                inline=True, inline_spacing=2, fontsize=FS - 3, fmt='%.1f')
		if not new:
			print(f'  UNLABELLED line {cl.levels[0]} ({cbar_label})')
		texts += new
	print(f'  lines drawn ({cbar_label}): {[float(c.levels[0]) for c in sets]}')
	for txt in texts:
		txt.set_path_effects([pe.withStroke(linewidth=3, foreground='white')])
		# OBTENEMOS LA POSICIÓN DEL TEXTO
		x_pos, y_pos = txt.get_position()
		
		# MÁRGENES DE SEGURIDAD (8% del tamaño del gráfico)
		margin_x = x_max * 0.1
		margin_y = y_max * 0.1
		
		# REPOSICIONAMIENTO INTELIGENTE: Si está en los bordes, cambiamos su anclaje
		if x_pos < margin_x:
			txt.set_horizontalalignment('left')   # Empuja a la derecha
		elif x_pos > (x_max - margin_x):
			txt.set_horizontalalignment('right')  # Empuja a la izquierda
				
	# no colour bar and no title per map: the two scales are common to all maps and drawn
	# once per lesion type (param_cbar-{axis,asym}.png), left map = lambda_r/lambda_ro, right map = J
	
	# Standardize Axes
	# no grid
	# tick numbers in percent, without the sign (the unit is in the axis label)
	pct = mtick.FuncFormatter(lambda v, _: f"{100 * v:.0f}")
	ax.xaxis.set_major_formatter(pct)
	ax.yaxis.set_major_formatter(pct)
	
	# Labels: the symbols only (their names are in the caption)
	ax.set_xlabel(r'$d_e$ (%)', fontsize=FS)
	ax.set_ylabel(r'$d_c$ (%)', fontsize=FS)
	ax.tick_params(labelsize=FS - 4)
	ax.set_aspect('equal')
	
	# Clean up chart borders (remove top and right spines for a modern scientific look)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)


parser = argparse.ArgumentParser()
# Removed 'type' positional argument to allow batch looping
parser.add_argument('--scatter', action='store_true', help='Use scatter plot instead of contours to view point density')
args = parser.parse_args()

# Define the folder where your FEBio file lives
work_dir = os.path.dirname(os.path.abspath(__file__)) + "/"

# List of all subfolders to process
'''
sub_folds = [
	'nu=1_alph=ON_zod=1',
	'nu=2_alph=ON_zod=1',
	'nu=3_alph=ON_zod=1',
	'nu=4_alph=ON_zod=1',
	'nu=5_alph=ON_zod=1',
	'nu=6_alph=ON_zod=1',
	#'nu=1_alph=OFF_zod=1',
	#'nu=2_alph=OFF_zod=1',
	#'nu=3_alph=OFF_zod=1',
	'nu=4_alph=OFF_zod=1',
	'nu=5_alph=OFF_zod=1',
	'nu=6_alph=OFF_zod=1',
	'nu=2_alph=ON_zod=1.2',
	'nu=2_alph=ON_zod=1.4',
	'nu=2_alph=ON_zod=1.6',
	'nu=2_alph=ON_zod=1.8',
	'nu=2_alph=ON_zod=2',
	'nu=2_alph=ON_zod=1_unit-frac',
]


sub_folds = [
	'nu=6_alph=ON_zod=1',
	'nu=6_alph=OFF_zod=1',
	'nu=5_alph=ON_zod=1.2',
	'nu=5_alph=ON_zod=1.4',
	'nu=5_alph=ON_zod=1.6',
	'nu=5_alph=ON_zod=1.8',
	'nu=5_alph=ON_zod=2',
	'nu=5_alph=ON_zod=2.2',
	'nu=5_alph=ON_zod=2.4',
	'nu=5_alph=ON_zod=2.6',
	'nu=5_alph=ON_zod=2.8',
	'nu=5_alph=ON_zod=3',
	'nu=5_alph=ON_zod=3.2',
	'nu=5_alph=ON_zod=4',
]
''';

sub_folds = [
'nu=1_alph=ON_zod=1',
'nu=2_alph=ON_zod=1',
'nu=2_alph=ON_zod=1.2',
'nu=2_alph=ON_zod=1.4',
'nu=2_alph=ON_zod=1.6',
'nu=2_alph=ON_zod=1.8',
'nu=2_alph=ON_zod=1_unit-frac',
'nu=2_alph=ON_zod=2',
'nu=2_alph=ON_zod=2.2',
'nu=2_alph=ON_zod=2.4',
'nu=2_alph=ON_zod=2.6',
'nu=2_alph=ON_zod=2.8',
'nu=2_alph=ON_zod=3',
'nu=2_alph=ON_zod=3.2',
'nu=2_alph=ON_zod=3.4',
'nu=2_alph=ON_zod=3.6',
'nu=2_alph=ON_zod=3.8',
'nu=2_alph=ON_zod=4',
'nu=2_alph=ON_zod=5',
'nu=3_alph=ON_zod=1',
'nu=3_alph=ON_zod=1.2',
'nu=3_alph=ON_zod=1.4',
'nu=3_alph=ON_zod=1.6',
'nu=3_alph=ON_zod=1.8',
'nu=3_alph=ON_zod=2',
'nu=3_alph=ON_zod=2.2',
'nu=3_alph=ON_zod=2.4',
'nu=3_alph=ON_zod=2.6',
'nu=3_alph=ON_zod=2.8',
'nu=3_alph=ON_zod=3',
'nu=4_alph=OFF_zod=1',
'nu=4_alph=ON_zod=1',
'nu=4_alph=ON_zod=1.2',
'nu=4_alph=ON_zod=1.4',
'nu=4_alph=ON_zod=1.6',
'nu=4_alph=ON_zod=1.8',
'nu=4_alph=ON_zod=2',
'nu=4_alph=ON_zod=2.2',
'nu=4_alph=ON_zod=2.4',
'nu=4_alph=ON_zod=2.6',
'nu=4_alph=ON_zod=2.8',
'nu=4_alph=ON_zod=3',
'nu=5_alph=OFF_zod=1',
'nu=5_alph=ON_zod=1',
'nu=5_alph=ON_zod=1.2',
'nu=5_alph=ON_zod=1.4',
'nu=5_alph=ON_zod=1.6',
'nu=5_alph=ON_zod=1.8',
'nu=5_alph=ON_zod=2',
'nu=5_alph=ON_zod=2.2',
'nu=5_alph=ON_zod=2.4',
'nu=5_alph=ON_zod=2.6',
'nu=5_alph=ON_zod=2.8',
'nu=5_alph=ON_zod=3',
'nu=5_alph=ON_zod=3.2',
'nu=5_alph=ON_zod=3.4',
'nu=5_alph=ON_zod=4',
'nu=6_alph=OFF_zod=1',
'nu=6_alph=ON_zod=1'
]


# Axisymmetric (axis) and Asymmetric (asym)
types = ['axis', 'asym']

# Define custom levels for J to match your previous setup
# contour lines (each one labelled) and colour ranges (the colour bars). The lambda_r range
# follows the data of the admissible states (0.75 <= J <= 4): max 1.91 fusiform, 3.06
# saccular over the 15 paper configurations -> 1-2 and 1-3.5; J reaches 4 in both
# (3.5 lies on the frontier itself). The fusiform admissible triangle is small (it ends
# at d_e + d_c ~ 35-60% on 0-70% axes), and four labelled lines do not fit in it at the
# printed label size: three lambda_r lines and two J lines there
LEVELS_J = {'asym': [1.5, 2.0, 2.5, 3.0], 'axis': [1.5, 2.5]}   # fusiform J=1.5 line ~1 label long
LEVELS_LR = {'asym': [1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0], 'axis': [1.2, 1.5, 1.8]}
CRANGE_J = (1.0, 4.0)
CRANGE_LR = {'axis': (1.0, 2.0), 'asym': (1.0, 3.5)}

# BATCH LOOP
sub_folds = ['nu=3_alph=ON_zod=1','nu=3_alph=ON_zod=1.2','nu=3_alph=ON_zod=1.4','nu=3_alph=ON_zod=1.6','nu=3_alph=ON_zod=1.8','nu=3_alph=ON_zod=2','nu=3_alph=ON_zod=2.2','nu=3_alph=ON_zod=2.4','nu=3_alph=ON_zod=2.6','nu=3_alph=ON_zod=2.8','nu=1_alph=ON_zod=1','nu=2_alph=ON_zod=1','nu=4_alph=ON_zod=1','nu=5_alph=ON_zod=1','nu=6_alph=ON_zod=1']  # paper figure set
for sub_fold in sub_folds:
	for type_ in types:
		print(f"\n==========================================")
		print(f"Processing: Folder = {sub_fold} | Type = {type_}")
		print(f"==========================================")

		file_EA = os.path.join(work_dir, sub_fold, f"summary_{type_}_elem_A.csv")
		file_EB = os.path.join(work_dir, sub_fold, f"summary_{type_}_elem_B.csv")
		file_NA = os.path.join(work_dir, sub_fold, f"summary_{type_}_node_A.csv")
		file_NB = os.path.join(work_dir, sub_fold, f"summary_{type_}_node_B.csv")

		# Check if files exist to prevent crashes
		if not os.path.exists(file_EA):
			print(f"WARNING: Files for {type_} in {sub_fold} not found. Skipping...")
			continue

		data_EA = np.array(np.loadtxt(file_EA, delimiter=',', comments='#'))
		data_EB = np.array(np.loadtxt(file_EB, delimiter=',', comments='#'))
		data_NA = np.array(np.loadtxt(file_NA, delimiter=',', comments='#'))
		data_NB = np.array(np.loadtxt(file_NB, delimiter=',', comments='#'))

		plot_mu = data_EA[:, 0]
		plot_cc = data_EA[:, 1]

		# Degradation (0 to 1)
		X = 1. - plot_mu
		Y = 1. - plot_cc

		# RETRIEVE ELEMENT DATA 
		s1A = data_EA[:, 3] # stress x 
		s1B = data_EB[:, 3] 
		s2A = data_EA[:, 4] # stress y 
		s2B = data_EB[:, 4] 
		s3A = data_EA[:, 5] # stress z 
		s3B = data_EB[:, 5] 
		JA = data_EA[:, 6] # VOLUME CHANGE (J)
		JB = data_EB[:, 6] 

		# RETRIEVE NODE DATA
		xA = data_NA[:, 3] # pos x 
		xB = data_NB[:, 3] 
		uxA = data_NA[:, 4] # disp x 
		uxB = data_NB[:, 4] 
		yA = data_NA[:, 5] # pos y 
		yB = data_NB[:, 5] 
		uyA = data_NA[:, 6] # disp y 
		uyB = data_NB[:, 6] 

		# RADIAL STRETCH Coordinates based on type
		if 'axis' in type_:
			xo = -2.31774149752e-15
			uo = 0.0
			yo =  0.605691058705
			vo = 0.0
		elif 'asym' in type_:
			xo = 0.640459751827
			uo = 0.0 #-0.0380821381617
			yo = 0.101438859344
			vo = 0.0 #-0.0060316181391
			
		lro_axis = radial_stretch(xo,uo,yo,vo)

		lr_axis_A  = radial_stretch(xA, uxA, yA, uyA)
		lr_ratio_A = lr_axis_A/lro_axis 

		lr_axis_B  = radial_stretch(xB, uxB, yB, uyB)
		lr_ratio_B = lr_axis_B/lro_axis 

		# PRESSURE
		press_A = (s1A + s2A + s3A)/3.
		press_B = (s1B + s2B + s3B)/3.

		# ==========================================
		# FILTER UNREALISTIC/UNSTABLE DATA
		# ==========================================
		eps = 1e-02
		# admissible states of the paper: converged, J <= 4 (no converged state has J < 1);
		# converged states with J > 4 are drawn as a grey band (create_contour_plot)
		J_max_0 = J_BAND
		J_max = J_max_0 + eps  

		lr_max = 3.
		mu_max = 1. 

		# Distinct y-axis constraints based on type
		cc_max = 1.0 
		# the fusiform admissible region never passes d_e = 63%, d_c = 47% over the 15 paper
		# configurations: both axes stop at 70%, the smallest round tick above it (every
		# axis ends on a labelled tick), the same in every fusiform figure; panels square
		if 'axis' in type_:
			mu_max, cc_max = 0.7, 0.7

		xy_max_J  = (mu_max, cc_max, J_max_0)   
		xy_max_lr = (mu_max, cc_max, lr_max)    

		# Create distinct masks for A
		mask_A = (JA <= J_max)

		# Clean Data for Point A
		X_A, Y_A = X[mask_A], Y[mask_A]
		press_A_clean = press_A[mask_A]
		lr_A_clean = lr_ratio_A[mask_A]
		JA_clean = JA[mask_A]

		# ==========================================
		# PLOTTING (point A only)
		# ==========================================
		print(f"Data stats for {type_} in {sub_fold}:")
		print('max J_A =', np.max(JA_clean) if len(JA_clean) > 0 else "N/A")
		print('min J_A =', np.min(JA_clean) if len(JA_clean) > 0 else "N/A")
		print('max lr_A =', np.max(lr_A_clean) if len(lr_A_clean) > 0 else "N/A")
		print('min lr_A =', np.min(lr_A_clean) if len(lr_A_clean) > 0 else "N/A")

		cmap_J  = 'viridis_r'
		cmap_lr = 'cividis_r'

		# Create a 2x1 grid
		fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10, 4.8))
		(ax2, ax3) = axes

		# --- ROW 1: POINT A ---
		create_contour_plot(ax2, X_A, Y_A, lr_A_clean, xy_max_lr, 
					r'$\lambda_{r}/\lambda_{ro}$', 
					levels=[l for l in LEVELS_LR[type_] if l < CRANGE_LR[type_][1]], 
					cmap=cmap_lr, smooth=True,
					use_scatter=args.scatter, crange=CRANGE_LR[type_], band=(X, Y, JA))
					
		create_contour_plot(ax3, X_A, Y_A, JA_clean,   xy_max_J,  
					r'$J$', 
					levels=LEVELS_J[type_], 
					cmap=cmap_J,
					use_scatter=args.scatter, crange=CRANGE_J, band=(X, Y, JA))

		plt.tight_layout(pad=0.2)
		
		# Save automatically without pausing
		plot_name = 'param_' + type_
		save_filename = f"{plot_name}-{sub_fold}.png"
		
		fig.savefig(save_filename, bbox_inches='tight', pad_inches=0.02, dpi=400)
		print(f"Saved figure: {save_filename}")

		# Close the figure to free memory for the next loop iteration
		plt.close(fig)

# the two colour scales, drawn once per lesion type and placed above each figure. Same
# text size as a map panel and 1.5 times its width (room for the ticks), so that, included at
# the scale of the panels (v2_build.py), bar text and map text print at the same size;
# bars ~2.5 mm tall.
import matplotlib as mpl
for t in ('axis', 'asym'):
	fig, cax = plt.subplots(1, 2, figsize=(15, 0.5), gridspec_kw=dict(wspace=0.25))
	for a, (cm, (c0, c1), lab) in zip(cax, (('cividis_r', CRANGE_LR[t], r'$\lambda_{r}/\lambda_{ro}$'),
	                                        ('viridis_r', CRANGE_J, r'$J$'))):
		cb = mpl.colorbar.ColorbarBase(a, cmap=plt.get_cmap(cm), orientation='horizontal',
		                               norm=mpl.colors.Normalize(c0, c1),
		                               ticks=np.arange(c0, c1 + 1e-9, 0.2 if c1 <= 2 else 0.5))
		cb.solids.set_alpha(0.8)
		cb.ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
		cb.ax.tick_params(labelsize=FS - 4)
		cb.set_label(lab, fontsize=FS, labelpad=4)      # centred under the bar
	fig.savefig(f'param_cbar-{t}.png', bbox_inches='tight', pad_inches=0.02, dpi=400)
	plt.close(fig)
print("\nBatch processing complete!")
