# Figures, tables and numbers: input, result, script

Paths are relative to `aneurysm-knockout/` unless they start with `paper/` or `semi-analytical/`.
Figure files are named as in the article source. Every script writes its figure to `paper/`.

## Figures

| figure file | inputs (`.feb`) | result files read | script |
|---|---|---|---|
| `mbe-scheme.pdf` | none | none | `paper/mbe-scheme.tex` (standalone TikZ, `pdflatex`) |
| `tension-4fib-all-again.pdf` | `uniax-fractions/uniax-MBE.feb`, `uniax-fractions/uniax-full-sload{1,10,100}-ko.feb` | `uniax-fractions/{MBE,sp1,sp10,sp100}_c.txt` (FE), `uniax-fractions/{MBE,sp1,sp10,sp100}_a.txt` (semi-analytical) | `uniax-fractions/make_fig_uniax_4fib.py` |
| `uniax-mbe-homeostasis.pdf` | `uniax-fractions/uniax-MBE8.feb`, `uniax-fractions/uniax-PASS8.feb` | `uniax-fractions/{elem1,node7}_{MBE8,PASS8}.txt` | `uniax-fractions/make_fig_homeostasis.py` |
| `matrix_J.pdf`, `matrix_energy.pdf` | `param-x/matrix-{fus,sac}-{el,elcc}.feb` | `param-x/matrix-*-NODES.txt` (deformed coordinates), `param-x/matrix-*-ALL.txt` (J, strain energy) | `param-x/render_matrix.py` |
| `aneurysm-4fib-axis.pdf` (fusiform) | `full-KG/apex/rate-fusiform-{MBE,full-sload10,full-sload100}.feb`, `full-KG/apex/ds-aneu/fus_ds001.feb` | `full-KG/apex/{elem,node}-axisym-{mbe,sp10,sp100}.txt`, `full-KG/apex/ds-aneu/{elem,node}-ds001.txt` | `full-KG/make_fig_aneurysm_4fib.py` |
| `aneurysm-4fib-asym.pdf` (saccular) | `full-KG/apex/rate-saccular-{MBE,full-sload10,full-sload100}.feb`, `full-KG/apex/ds-aneu/sac_ds001.feb` | `full-KG/apex/{elem,node}-axis-{FINALmbe,sp10,sp100}.txt`, `full-KG/apex/ds-aneu/{elem,node}-axis-FINALfk-SDS.txt` | `full-KG/make_fig_aneurysm_4fib.py` |
| `nu-sensitivity.pdf` | parametric sweep, see below | `param-x/MASTER_PARAMETRIC_STUDY.csv` | `param-x/make_fig_nu_sensitivity.py` |
| `fig-param-s/{axis,asym}-nu=*.png`, `fig-param-s/cbar-*.png` (Online Resource 1) | parametric sweep, see below | `param-x/nu=*_alph=ON_zod=*/summary_{axis,asym}_{elem,node}_{A,B}.csv` | `param-x/plot_febio_many-4.py`, then `paper/shrink_maps.py` |

The rate-comparison runs with `s_load = 1` day are the `ds-aneu/*_ds001` inputs. The
decks in `full-KG/apex/` log the lesion apex: element 10 and node 20.

`plot_febio_many-4.py` writes `param_<axis|asym>-<folder>.png` and `param_cbar-<axis|asym>.png`
to `param-x/`. The article uses reduced copies. To regenerate them, copy each image to
`paper/fig-param/` without the `param_` prefix and run `paper/shrink_maps.py`, which writes
`paper/fig-param-s/` at one third of the width with a 256-colour palette. The
full-resolution `paper/fig-param/` is not included.

## Parametric sweep

| item | file |
|---|---|
| templates | `param-x/TAA-axis-MBE-BASE.feb` (fusiform), `param-x/TAA-asym-MBE-BASE.feb` (saccular) |
| configurations in the article | 15 folders `param-x/nu=<nu>_alph=ON_zod=<s_z>`: `nu` = 1 to 6 at `s_z` = 1 mm; `s_z` = 1.0 to 2.8 mm at `nu` = 3 |
| sweep of one configuration | `param-x/run_febio_many.py` (20 x 20 grid in `d_e`, `d_c`), `param-x/densify_runner.py` (39 x 39 grid: the 20 x 20 points plus the midpoints) |
| saccular rerun with the apex probe | `param-x/rerun_asym_apex.py` |
| collected table | `param-x/rebuild_master.py` → `param-x/MASTER_PARAMETRIC_STUDY.csv` |

`MASTER_PARAMETRIC_STUDY.csv` also holds rows of exploratory configurations that are
not in the article. `paper/make_numbers.py` and `param-x/make_fig_nu_sensitivity.py`
select the configurations of the article from it.

## Tables and numbers in the text

| item | inputs (`.feb`) | result files read | script |
|---|---|---|---|
| convergence tables, CPU time, step and Newton-iteration counts | `timing/conv-uniax-{MBE,HYP}.feb`, `timing/u_ds0.5.feb`, `timing/conv-aneu-{MBE,HYP}.feb`, `timing/a_ds1.0.feb` | `timing/*.log`, `timing/timing.csv` | `timing/timing.py` (one thread, repeated runs) |
| time-step study, uniaxial | `uniax-fractions/ds-study/u_ds*.feb` | `uniax-fractions/ds-study/elem1_ds*.txt` | none |
| time-step study, aneurysm | `full-KG/ds-aneu/fus_ds{01,005,001}.feb`, `full-KG/ds-study/a_ds1.0.feb` | `full-KG/ds-aneu/elem-ds*.txt`, `full-KG/ds-study/elem-axisym-sp100_a_ds1.0.txt` | none |
| MBE step count, fusiform rate comparison | `full-KG/rate-fusiform-MBE.feb` | `full-KG/rate-fusiform-MBE.log` | none |
| volumetric penalty `k_J` | `kJ-study/{uniax,aneu}_kJ1e{3..8}.feb` | `kJ-study/uniax_kJ*_elem1_MBE.txt`, `kJ-study/aneu_kJ*_matrix-fus-el-{eA,eB}.txt` | `kJ-study/run_kJ.py` |
| radial mesh refinement | `mesh-study/{fus,sac}_nr{1,2,3}.feb` | `mesh-study/{fus,sac}_nr*_matrix-*-el-ALL.txt` | `mesh-study/run_clinical.py`, `mesh-study/gen_mesh.py` |
| uniaxial homeostasis, semi-analytical agreement | as for `tension-4fib-all-again.pdf` and `uniax-mbe-homeostasis.pdf` | same | `semi-analytical/uniax_4fib.py` |
| field-map values | `param-x/matrix-*.feb` | `param-x/matrix-*-{ALL,NODES,eA,nA}.txt` | none |

`paper/make_numbers.py` reads the result files listed in `paper/numbers-inputs.txt`
and writes every number of the text as a LaTeX macro in `paper/numbers.tex`.

`kJ-study/run_kJ.py` writes the `kJ-study/*.feb` from `param-x/matrix-fus-el.feb` and
`uniax-fractions/uniax-MBE.feb`, changing only `k_J`, and runs them.

## Result-file names

| file name contains | written by |
|---|---|
| `-axisym-` | a fusiform (axisymmetric) input |
| `-axis-`, `-FINALfk-SDS` | a saccular input |
| `elem1_`, `node7_` | a uniaxial coupon input |

A `.feb` names its result files in its `<logfile>` section. The file name
`node-axis-min_sp1.txt` is written by three saccular inputs in `full-KG/apex/`, and
no script reads it.
