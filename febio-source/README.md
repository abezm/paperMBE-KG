# FEBio source changes

`mbe-kg-febio-4.12.patch` applies to FEBio tag `v4.12`
(`https://github.com/febiosoftware/FEBio`, commit `9ea80f8`).

| file | change |
|---|---|
| `FEBioMech/FEKinematicGrowthZ.{h,cpp}` | new: the `KinematicGrowth` material (stress, tangent, equilibrated and rate-dependent growth) |
| `FEBioMech/FEBioMechModule.cpp` | registration of the material, of three plot variables and of nine log variables |
| `FEBioMech/FEBioMechPlot.{h,cpp}` | plot variables `elastin damage`, `collagen damage`, `current strain energy density` |
| `FEBioMech/FEBioMechData.{h,cpp}` | element log variables `Fgxx` to `Fgzy` (growth deformation gradient) |
| `FEBioMech/FEElasticMaterialPoint.{h,cpp}` | material-point variable `m_svo`: homeostatic volumetric stress |
| `FEBioMech/FESolidMaterial.{h,cpp}` | `TangentMm`: tangent with minor symmetries only |
| `FECore/mat3d.hpp` | dyadic products `dyad1`, `dyad1s` |

`FEKinematicGrowthZ.h` and `FEKinematicGrowthZ.cpp` are also included as plain files.

Stage I (`t <= gr_start`) of `FEKinematicGrowthZ.cpp`: F_g = 1 and S = J S_e + penalty term. Fibre angle
update (`alp_update` = 1): the stretches enter as F F_o⁻¹, and the tangent includes the variation of the
stress magnitude and of the volumetric multiplier with the angle.

`tancheck-fd.patch` (optional, applied after `mbe-kg-febio-4.12.patch`) adds a check of the tangent
against central finite differences, for the inputs in `aneurysm-knockout/tancheck/`. The material
parameters `fd_check`, `fd_time`, `fd_elem`, `fd_gp` and `fd_eps` select the time, element, integration
point and relative perturbation. Three checks are written to the FEBio log: material tangent (`FDCHK`),
element stiffness (`EFDCHK`) and global Jacobian along the Newton direction (`GFDCHK`). Files:
`FEBioMech/FEKinematicGrowthZ.{h,cpp}`, `FEBioMech/FEElasticSolidDomain.{h,cpp}`,
`FEBioMech/FESolidSolver2.cpp`. The patch adds code only.

`mkl2025-compat.patch` replaces the sparse triangular solve `mkl_dcsrtrsv`, which
oneAPI MKL 2025 removed, in `NumCore/ILU0_Preconditioner.cpp`,
`NumCore/ILUT_Preconditioner.cpp` and `NumCore/IncompleteCholesky.cpp`. FEBio v4.12
does not compile against MKL 2025 or later without it. The inputs of this repository
use PARDISO and do not call these preconditioners.
