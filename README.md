# MBE-KG: mechanobiologically-equilibrated finite kinematic growth

Source code, finite-element inputs, results and post-processing scripts for

> A. Bezmalinovic, C. García-Herrera, D. Celentano, M. Latorre.
> *A mechanobiologically-equilibrated finite kinematic growth model for soft-tissues:
> formulation and fast, rate-independent, finite-element implementation.*

MBE-KG adds mechanobiological equilibrium to classical kinematic growth (KG) as an
algebraic constraint, solved together with mechanical equilibrium. The grown
state follows from the constraint directly, with no time integration of the growth
transient. The model is implemented as a material in a source-modified FEBio 4.12.

## Contents

All simulations of the article are in `aneurysm-knockout/s1/`.

| directory | contents |
|---|---|
| `febio-source/` | FEBio patch, the `KinematicGrowth` material source, optional finite-difference tangent check, build notes |
| `semi-analytical/` | `uniax_4fib.py`: semi-analytical solution of the uniaxial coupon, independent of FEBio |
| `aneurysm-knockout/s1/alp1/uniax-lo1/` | uniaxial coupon: equilibrated and rate-dependent runs, homeostasis check, time-step study (`ds-study/`), volumetric penalty (`kJ-study/`), convergence and CPU time (`timing-stop/`) |
| `aneurysm-knockout/s1/uniax-fractions/` | uniaxial coupon without growth (`uniax-PASS8.feb`), homeostasis check |
| `aneurysm-knockout/s1/param-x/` | field maps of the four lesions, parametric study of the degradation profile, retries at the edge of the admissible region (`edge-retry/`) |
| `aneurysm-knockout/s1/severity-fig/` | lesions with equal losses `d_e = d_c` |
| `aneurysm-knockout/s1/profile-fig/` | lesions with other spread parameters `s_z` and exponential decay parameters `nu` |
| `aneurysm-knockout/s1/alp2/` | rate-dependent vs equilibrated comparison of the fusiform and saccular lesions (`full-KG/`), convergence and CPU times (`timing/`, `speedup-fair/`, `speedup-fair-sac/`) |
| `aneurysm-knockout/s1/alp1/full-KG/apex/` | equilibrated lesions of the rate comparison with elastin degradation only (`cc_in` = 1) |
| `aneurysm-knockout/s1/kJ-study/` | volumetric penalty `k_J`, lesion |
| `aneurysm-knockout/s1/mesh-study/`, `aneurysm-knockout/s1/mesh-study2/` | radial mesh refinement; axial refinement and segment length |
| `aneurysm-knockout/tancheck/` | tangent check against central finite differences (needs `febio-source/tancheck-fd.patch`) |
| `paper/` | `make_numbers.py` and its output `numbers.tex`; the figure files of the article |

[`aneurysm-knockout/FIGURE-FILES.md`](aneurysm-knockout/FIGURE-FILES.md) lists, for
every figure and table, the input files, the result files and the script.

## FEBio build

Every input declares `<material type="KinematicGrowth">`. Released FEBio does not
contain this material, so the inputs need FEBio 4.12 rebuilt with the patch in
`febio-source/`:

```bash
git clone --branch v4.12 https://github.com/febiosoftware/FEBio.git
cd FEBio
git apply /path/to/febio-source/mbe-kg-febio-4.12.patch
git apply /path/to/febio-source/mkl2025-compat.patch   # only with oneAPI MKL 2025 or later
git apply /path/to/febio-source/tancheck-fd.patch      # only for the inputs in aneurysm-knockout/tancheck/
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_MKL=ON -DUSE_ZLIB=ON -DUSE_STATIC_STDLIBS=ON \
         -DMKLROOT=/opt/intel/oneapi/mkl -DMKL_INC=/opt/intel/oneapi/mkl/latest/include \
         -DMKL_LIB_DIR=/opt/intel/oneapi/mkl/latest/lib/intel64 \
         -DMKL_OMP_LIB=/opt/intel/oneapi/compiler/latest/lib/libiomp5.so
make -j
export FEBIO=$PWD/bin/febio4
```

`USE_MKL=ON` is required. No input declares a `<linear_solver>`, so FEBio uses its
default solver, PARDISO, which is part of MKL.

`febio-source/README.md` lists the modified files.

## Model switches in the inputs

| tag | value | meaning |
|---|---|---|
| `MBE` | `1` | equilibrated growth (MBE-KG) |
| `MBE` | `0` | rate-dependent growth (full KG), rate `k_parameter` = K_σ k [1/day] |
| `gr_start` | beyond the end time | growth never starts: hyperelastic response |
| `lo_opt` | `1` | F_e = F_o after `gr_start` |
| `alp_update` | `1` | diagonal collagen angle updated during growth |
| `ec_degrad` | `1` | `mu_in`, `cc_in` and the degradation profile applied after `gr_start` |
| `mu_in`, `cc_in` | < 1 | elastin and collagen degradation at the lesion centre |
| `z_spread`, `z_nu`, `t_nu`, `t_pif` | | axial and circumferential degradation profile |
| `asymm_aneurysm` | `0` / `1` | fusiform (axisymmetric) / saccular lesion |

`lo_opt`, `alp_update` and `ec_degrad` are `1` in every input, except `uniax-fractions/uniax-PASS8.feb` (no growth)
and `alp1/uniax-lo1/timing-stop/conv-uniax-HYP.feb` (hyperelastic), where `lo_opt` and `ec_degrad` are `0`.

Material sets (paths relative to `aneurysm-knockout/s1/`):

| inputs | `Ge_t`, `Ge_z` | `Gm`, `Gc` | `f_e`, `f_c` |
|---|---|---|---|
| `alp1/uniax-lo1/`, `uniax-fractions/` | 1.0, 1.0 | 1.0, 1.0 | 0.34, 0.33 |
| `param-x/`, `severity-fig/`, `profile-fig/`, `kJ-study/`, `mesh-study/`, `mesh-study2/` | 1.9, 1.62 | 1.2, 1.25 | 0.34, 0.33 |
| `alp2/`, `alp1/full-KG/` | 1.0, 1.0 | 1.2, 1.25 | 0.45, 0.45 |

The third set is used where the rate-dependent model runs: with the second set the
rate-dependent model loses convergence by element inversion once the lesion dilates
the wall.

## Running

Each `.feb` writes its result files to the directory it is run from:

```bash
cd aneurysm-knockout/s1/alp1/uniax-lo1
$FEBIO -i uniax-MBE.feb -o uniax-MBE.log
```

The run and timing scripts (`s1/profile-fig/run_all.sh`, `s1/alp2/run_all.sh`, `s1/param-x/run_maps.py`,
`s1/param-x/edge-retry/run_edge_s1.py`, `s1/alp1/uniax-lo1/timing-stop/run_stop.py`,
`s1/alp2/speedup-fair/timing_fair2.py`, `s1/alp2/speedup-fair/stop/timing_stop.py`,
`s1/alp2/speedup-fair-sac/clean/timing_sac.py`) read the FEBio executable from the environment variable
`FEBIO`, with `febio4` on the `PATH` as the default.

The homeostatic volumetric stress `sigma_vo` is stored at the last time step with
`t <= gr_start`. Every input therefore carries a must-point at `gr_start` in its
time-stepper load curve. Without it the stored stress belongs to a partly loaded
state, and the run still terminates normally.

## Numbers and figures without FEBio

All result files read by the article are included, so the numbers and the figures are
regenerated from them without running FEBio. Python 3.10 with numpy, scipy,
matplotlib, Pillow and pyvista:

```bash
cd paper && python3 make_numbers.py                                   # numbers.tex, numbers-inputs.txt
cd aneurysm-knockout/s1/alp1/uniax-lo1 && python3 prep/make_fig_uniax_4fib_s1.py && python3 make_fig_homeostasis_s1.py
cd aneurysm-knockout/s1/param-x && python3 render_matrix_iso.py && python3 make_fig_nu_sensitivity_s1.py
cd aneurysm-knockout/s1/profile-fig && for v in lad sz nu; do python3 render_profile_matrix.py $v; done
cd aneurysm-knockout/s1/alp2/full-KG && python3 make_fig_aneurysm_4fib_s1.py
cd aneurysm-knockout/s1/param-x && python3 plot_maps_s1.py \
  && cd ../../../paper && python3 shrink_maps.py ../aneurysm-knockout/s1/param-x/maps-s1b fig-param-s1b-s
```

Three scripts do not overwrite an existing figure in `paper/`: `make_fig_uniax_4fib_s1.py` needs
`--overwrite`, `make_fig_homeostasis_s1.py` needs `--overwrite-v2`, and `make_fig_nu_sensitivity_s1.py`
needs `paper/nu-sensitivity-s1.pdf` removed first.

`numbers.tex` holds every number of the article as a LaTeX macro. `numbers-inputs.txt`
lists the result files that `make_numbers.py` reads.

## License

MIT, see `LICENSE`. The files in `febio-source/` modify FEBio, which is distributed
under the MIT license by the University of Utah, Columbia University and others.
