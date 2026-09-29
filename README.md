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

| directory | contents |
|---|---|
| `febio-source/` | FEBio patch, the `KinematicGrowth` material source, build notes |
| `semi-analytical/` | `uniax_4fib.py`: semi-analytical solution of the uniaxial coupon, independent of FEBio |
| `aneurysm-knockout/uniax-fractions/` | uniaxial coupon: equilibrated and rate-dependent runs, homeostasis check, time-step study (`ds-study/`) |
| `aneurysm-knockout/full-KG/` | rate-dependent vs equilibrated comparison for the fusiform and saccular lesions; probe at the lesion apex (`apex/`) |
| `aneurysm-knockout/param-x/` | aneurysm field maps and the parametric study of the degradation profile |
| `aneurysm-knockout/timing/` | convergence tables and CPU times |
| `aneurysm-knockout/kJ-study/` | volumetric penalty `k_J` |
| `aneurysm-knockout/mesh-study/` | radial mesh refinement |
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
| `mu_in`, `cc_in` | < 1 | elastin and collagen degradation at the lesion centre |
| `z_spread`, `z_nu`, `t_nu`, `t_pif` | | axial and circumferential degradation profile |
| `asymm_aneurysm` | `0` / `1` | fusiform (axisymmetric) / saccular lesion |

Material sets:

| inputs | `Ge_t`, `Ge_z` | `Gm`, `Gc` | `f_e`, `f_c` |
|---|---|---|---|
| `uniax-fractions/`, `timing/conv-uniax-*`, `kJ-study/uniax_*` | 1.0, 1.0 | 1.0, 1.0 | 0.34, 0.33 |
| `param-x/`, `mesh-study/`, `kJ-study/aneu_*` | 1.9, 1.62 | 1.2, 1.25 | 0.34, 0.33 |
| `full-KG/`, `timing/conv-aneu-*`, `timing/a_ds1.0` | 1.0, 1.0 | 1.2, 1.25 | 0.45, 0.45 |

The third set is used where the rate-dependent model runs: with the second set the
rate-dependent model loses convergence by element inversion once the lesion dilates
the wall.

## Running

Each `.feb` writes its result files to the directory it is run from:

```bash
cd aneurysm-knockout/uniax-fractions
$FEBIO -i uniax-MBE.feb -o uniax-MBE.log
```

Study scripts (`timing/timing.py`, `kJ-study/run_kJ.py`, `mesh-study/run_clinical.py`,
`param-x/densify_runner.py`, `param-x/rerun_asym_apex.py`, `full-KG/apex/run_all.sh`)
read the FEBio executable from the environment variable `FEBIO`, with `febio4` on
the `PATH` as the default.

The homeostatic volumetric stress `sigma_vo` is stored at the last time step with
`t <= gr_start`. Every input therefore carries a must-point at `gr_start` in its
time-stepper load curve. Without it the stored stress belongs to a partly loaded
state, and the run still terminates normally.

## Numbers and figures without FEBio

All result files of the article are included, so the numbers and the figures are
regenerated from them without running FEBio. Python 3.10 with numpy, scipy,
matplotlib, Pillow and pyvista:

```bash
cd paper && python3 make_numbers.py                  # numbers.tex, numbers-inputs.txt
cd aneurysm-knockout/uniax-fractions && python3 make_fig_uniax_4fib.py && python3 make_fig_homeostasis.py
cd aneurysm-knockout/full-KG && python3 make_fig_aneurysm_4fib.py
cd aneurysm-knockout/param-x && python3 make_fig_nu_sensitivity.py && python3 render_matrix.py
```

`numbers.tex` holds every number of the article as a LaTeX macro. `numbers-inputs.txt`
lists the result files that `make_numbers.py` reads.

## License

MIT, see `LICENSE`. The files in `febio-source/` modify FEBio, which is distributed
under the MIT license by the University of Utah, Columbia University and others.
