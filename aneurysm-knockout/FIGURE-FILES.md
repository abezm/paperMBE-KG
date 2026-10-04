# Figures, tables and numbers: input, result, script

Paths are relative to `aneurysm-knockout/s1/` unless they start with `paper/`, `semi-analytical/`
or `aneurysm-knockout/`. Figure files are named as in the article source. The scripts of Figures 2, 3,
9, 10 and 11 write the figure to `paper/`; `render_matrix_iso.py` (Figures 4, 5) and
`render_profile_matrix.py` (Figures 6-8) write it next to the script.

## Figures

| figure | figure file | inputs (`.feb`) | result files read | script |
|---|---|---|---|---|
| 1 | `mbe-scheme.pdf` | none | none | `paper/mbe-scheme.tex` (standalone TikZ, `pdflatex`) |
| 2 | `tension-4fib-s1.pdf` | `alp1/uniax-lo1/uniax-MBE.feb`, `alp1/uniax-lo1/uniax-full-sload{1,100}-ko.feb`, `alp1/uniax-lo1/to5/*.feb` (`s_load` = 10 days) | `alp1/uniax-lo1/prep/{MBE,sp1,sp10,sp100}_c.txt` (FE), `alp1/uniax-lo1/prep/{MBE,sp1,sp10,sp100}_a.txt` (semi-analytical) | `alp1/uniax-lo1/prep/make_fig_uniax_4fib_s1.py` |
| 3 | `uniax-mbe-homeostasis-s1.pdf` | `alp1/uniax-lo1/uniax-MBE8.feb`, `uniax-fractions/uniax-PASS8.feb` | `alp1/uniax-lo1/{elem1,node7}_MBE8.txt`, `uniax-fractions/{elem1,node7}_PASS8.txt` | `alp1/uniax-lo1/make_fig_homeostasis_s1.py` |
| 4, 5 | `matrix_J_iso-s1.pdf`, `matrix_energy_iso-s1.pdf` | `param-x/matrix-{fus,sac}-{el,elcc}.feb` | `param-x/matrix-*-NODES.txt` (deformed coordinates), `param-x/matrix-*-ALL.txt` (J, strain energy) | `param-x/render_matrix_iso.py` |
| 6 | `ladder-lam-s1-common.pdf` | `severity-fig/lad-{fus,sac}-d*.feb` (8 lesions) | `severity-fig/lad-*-{NODES,ALL}.txt` | `profile-fig/render_profile_matrix.py lad` |
| 7 | `profile-sz-lam-s1.pdf` | `profile-fig/prof-{fus,sac}-nu5-sz{1p2,1p8,2p4}.feb`, `severity-fig/lad-fus-d20p000-v5.feb`, `severity-fig/lad-sac-d37p000-lev2.feb` | `profile-fig/prof-*-{NODES,ALL}.txt`, `severity-fig/lad-*-{NODES,ALL}.txt` | `profile-fig/render_profile_matrix.py sz` |
| 8 | `profile-nu-lam-s1.pdf` | `profile-fig/prof-{fus,sac}-nu{1,3,7}-sz3p0.feb`, the two `severity-fig` lesions of Figure 7 | as Figure 7 | `profile-fig/render_profile_matrix.py nu` |
| 9 | `aneurysm-4fib-axis-s1.pdf` (fusiform) | `alp2/full-KG/apex/rate-fusiform-MBE.feb`, `alp2/full-KG/apex/ds-aneu-t1/fus_ds001.feb` (`s_load` = 1 day), `alp2/full-KG/apex/to5/rate-fusiform-full-sload{10,100}.feb` | `alp2/full-KG/apex/{node,elem}-axisym-mbe.txt`, `alp2/full-KG/apex/ds-aneu-t1/{node,elem}-ds001.txt`, `alp2/full-KG/apex/to5/{node,elem}-axisym-sp{10,100}.txt` | `alp2/full-KG/make_fig_aneurysm_4fib_s1.py` |
| 10 | `aneurysm-4fib-asym-s1.pdf` (saccular) | `alp2/full-KG/apex/rate-saccular-MBE.feb`, `alp2/full-KG/apex/ds-aneu-t1/sac_ds001.feb` (`s_load` = 1 day), `alp2/full-KG/apex/rate-saccular-full-sload10.feb`, `alp2/full-KG/apex/to5/rate-saccular-full-sload100.feb` | `alp2/full-KG/apex/{node,elem}-axis-FINALmbe.txt`, `alp2/full-KG/apex/ds-aneu-t1/{node,elem}-axis-FINALfk-SDS.txt`, `alp2/full-KG/apex/{node,elem}-axis-sp10.txt`, `alp2/full-KG/apex/to5/{node,elem}-axis-sp100.txt` | `alp2/full-KG/make_fig_aneurysm_4fib_s1.py` |
| 11 | `nu-sensitivity-s1.pdf` | parametric sweep, see below | `param-x/MASTER_PARAMETRIC_STUDY.csv` | `param-x/make_fig_nu_sensitivity_s1.py` |
| S1-S4 | `fig-param-s1b-s/{axis,asym}-nu=*.png`, `fig-param-s1b-s/cbar-*.png` (Online Resource 1) | parametric sweep, see below | `param-x/nu=*_alph=ON_zod=*/summary_{axis,asym}_{elem,node}_{A,B}.csv` | `param-x/plot_maps_s1.py`, then `paper/shrink_maps.py` |

The decks in `alp2/full-KG/apex/` log the lesion apex: element 10 and node 20.

`plot_maps_s1.py` writes `param_<axis|asym>-<folder>.png` and `param_cbar-<axis|asym>.png` to
`param-x/maps-s1b/`. `paper/shrink_maps.py ../aneurysm-knockout/s1/param-x/maps-s1b fig-param-s1b-s`
writes the copies of the article: one third of the width, 256-colour palette, prefix `param_` removed.

`profile-fig/make_decks.py` writes the `prof-*.feb` from the two `severity-fig` lesions of Figure 7,
changing only `z_spread`, `z_nu` and `t_nu`. `profile-fig/run_all.sh` runs them.

## Parametric sweep

| item | file |
|---|---|
| templates | `param-x/TAA-axis-MBE-BASE.feb` (fusiform), `param-x/TAA-asym-MBE-BASE.feb` (saccular) |
| profiles of the article | `PAPER_PROFILES` in `param-x/run_maps.py`: `s_z` = 1.0, 1.4, 1.8, 2.2, 2.6, 3.0 mm at `nu` = 5; `nu` = 1 to 6 at `s_z` = 3 mm |
| sweep of one profile | `param-x/run_maps.py` (20 x 20 grid in `d_e`, `d_c`, then the midpoints), with `aneurysm-knockout/param-x/densify_runner.py` |
| retries at the edge of the admissible region | `param-x/edge-retry/run_edge_s1.py`, `param-x/edge-retry/collect.py` → `param-x/edge-retry/results.json` |
| collected table | `param-x/rebuild_master_s1.py` → `param-x/MASTER_PARAMETRIC_STUDY.csv` |

`MASTER_PARAMETRIC_STUDY.csv` also holds the profiles at `nu` = 3 and the profile at `s_z` = 2.8 mm,
which are computed and not shown. `paper/make_numbers.py` and `param-x/make_fig_nu_sensitivity_s1.py`
select the profiles of the article from it.

## Tables and numbers in the text

| item | inputs (`.feb`) | result files read | script |
|---|---|---|---|
| convergence table and CPU time, uniaxial | `alp1/uniax-lo1/timing-stop/conv-uniax-{MBE,HYP}.feb`, `alp1/uniax-lo1/timing-stop/u_ds0.5.feb`, `alp1/uniax-lo1/timing-stop/find/u_ds0.5-long.feb` (stop time) | `alp1/uniax-lo1/timing-stop/*.log`, `alp1/uniax-lo1/timing-stop/timing.csv` | `alp1/uniax-lo1/timing-stop/run_stop.py` |
| convergence table, fusiform lesion | `alp2/timing/conv-aneu-MBE.feb`, `alp2/timing/a_ds1.0.feb` | `alp2/timing/*.log` | none |
| convergence table, saccular lesion | `alp2/speedup-fair-sac/loaded/T-{full-ds5,mbe-1p10}-sac-r1.feb` | `alp2/speedup-fair-sac/loaded/T-*-sac-r1.log` | `alp2/speedup-fair-sac/loaded/timing_sac.py` |
| step and Newton-iteration counts, CPU time, fusiform lesion | `alp2/speedup-fair/*.feb`, `alp2/speedup-fair/stop/*.feb`, `alp2/speedup-fair/stop/find/*.feb` (stop time) | `alp2/speedup-fair/timing_fair2.csv`, `alp2/speedup-fair/stop/timing_stop.csv`, `alp2/speedup-fair/stop/T-*_elem-fus-sl100-2st-ds*.txt` | `alp2/speedup-fair/timing_fair2.py`, `alp2/speedup-fair/stop/timing_stop.py`, `alp2/make_stop.py` |
| step and Newton-iteration counts, CPU time, saccular lesion | `alp2/speedup-fair-sac/{clean,loaded}/*.feb`, `alp2/speedup-fair-sac/find/*.feb` (stop time) | `alp2/speedup-fair-sac/loaded/timing_sac.csv` (steps, iterations), `alp2/speedup-fair-sac/clean/timing_sac.csv` (CPU time) | `alp2/speedup-fair-sac/setup_sac.py`, `alp2/speedup-fair-sac/clean/timing_sac.py` |
| time-step study, uniaxial | `alp1/uniax-lo1/ds-study/u_ds*.feb` | `alp1/uniax-lo1/ds-study/elem1_ds*.txt`, `alp1/uniax-lo1/ds-study/u_ds1.0t1.log` | none |
| time-step study, lesions | `alp2/full-KG/apex/ds-aneu-t1/{fus_ds01,fus_ds005,fus_ds001,sac_ds01}.feb` | `alp2/full-KG/apex/ds-aneu-t1/elem-*.txt` | none |
| MBE step count, fusiform rate comparison | `alp2/full-KG/apex/rate-fusiform-MBE.feb` | `alp2/full-KG/apex/rate-fusiform-MBE.log` | none |
| rate-comparison lesions with elastin degradation only (`cc_in` = 1) | `alp1/full-KG/apex/rate-{fusiform,saccular}-MBE.feb` | `alp1/full-KG/apex/{node,elem}-*.txt` | none |
| volumetric penalty `k_J` | `alp1/uniax-lo1/kJ-study/*.feb`, `kJ-study/aneu_kJ1e{3..8}.feb` | `alp1/uniax-lo1/kJ-study/uniax_kJ*_elem1_MBE.txt`, `kJ-study/aneu_kJ*_matrix-fus-el-{eA,eB}.txt` | none |
| radial mesh refinement | `mesh-study/{fus,sac}_nr{1,2,3}.feb` | `mesh-study/{fus,sac}_nr*_matrix-*-el-ALL.txt` | none |
| axial mesh refinement, segment length | `mesh-study2/{gate,z2,L30}/*.feb` | `mesh-study2/summary.csv` | `mesh-study2/gen.py`, `mesh-study2/run.sh`, `mesh-study2/summary.py` |
| lesions with equal losses | `severity-fig/lad-*.feb` | `severity-fig/lad-*-{ALL,NODES}.txt` | `severity-fig/make_ladder_csv.py` |
| uniaxial homeostasis, semi-analytical agreement | as for Figures 2 and 3 | same | `semi-analytical/uniax_4fib.py` |
| field-map values | `param-x/matrix-*.feb` | `param-x/matrix-*-{ALL,NODES,eA,nA}.txt` | none |
| tangent check against central finite differences | `aneurysm-knockout/tancheck/*/*.feb` (`fd_check` = 1) | `aneurysm-knockout/tancheck/*/*.log` | `aneurysm-knockout/tancheck/setup.py`, `aneurysm-knockout/tancheck/setup_r8.py` |

`paper/make_numbers.py` reads the result files listed in `paper/numbers-inputs.txt`
and writes every number of the text as a LaTeX macro in `paper/numbers.tex`.
It also reads `aneurysm-knockout/severity-fig/ladder-numbers.json`, from runs before `s1/`, for macros
that the article does not print.

The inputs in `aneurysm-knockout/tancheck/` run only with FEBio rebuilt with `febio-source/tancheck-fd.patch`.
The check writes its finite-difference comparison to the FEBio log.

## Result-file names

| file name contains | written by |
|---|---|
| `-axisym-` | a fusiform (axisymmetric) input |
| `-axis-`, `-FINALfk-SDS` | a saccular input |
| `elem1_`, `node7_` | a uniaxial coupon input |

A `.feb` names its result files in its `<logfile>` section.
