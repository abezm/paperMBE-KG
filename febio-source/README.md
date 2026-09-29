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

`mkl2025-compat.patch` replaces the sparse triangular solve `mkl_dcsrtrsv`, which
oneAPI MKL 2025 removed, in `NumCore/ILU0_Preconditioner.cpp`,
`NumCore/ILUT_Preconditioner.cpp` and `NumCore/IncompleteCholesky.cpp`. FEBio v4.12
does not compile against MKL 2025 or later without it. The inputs of this repository
use PARDISO and do not call these preconditioners.
