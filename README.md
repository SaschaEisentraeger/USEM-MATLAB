# USEM: unsymmetric spectral element method (2D elasticity)

[![tests](https://github.com/SaschaEisentraeger/USEM-MATLAB/actions/workflows/tests.yml/badge.svg)](https://github.com/SaschaEisentraeger/USEM-MATLAB/actions/workflows/tests.yml)

MATLAB code for **high-order unsymmetric (Petrov–Galerkin) spectral elements** for linear
plane elasticity: static, transient, modal and frequency-domain analysis.

The element combines two sets of shape functions:

- **parametric test functions:** Lagrange polynomials of degree $p$ on Gauss–Lobatto–Legendre
  nodes, which keep the mesh conforming;
- **metric trial functions:** Lagrange polynomials in corotational **skew coordinates**, which
  are complete in Cartesian coordinates and frame-invariant.

The result reproduces every displacement field of total degree $\le p$ **exactly on
arbitrarily distorted meshes** (angular and curved-edge distortion). It generalizes the
unsymmetric 8-node element of Rajendran & Liew (2003) and Rajendran (2010) to arbitrary
polynomial degree.

## Features

- Quadrilateral spectral elements of arbitrary degree $p_\xi$, $p_\eta$ (GLL, equidistant or
  Chebyshev nodes), lifted from 4-node or 8-node ABAQUS meshes.
- Unsymmetric stiffness and mass matrices: consistent, GLL nodal quadrature (diagonal),
  row-sum and HRZ lumping.
- Isotropic and orthotropic materials, plane stress and plane strain.
- Loads: nodal forces, pressure and traction (also spatially varying), body forces, gravity.
  Dirichlet conditions: homogeneous, prescribed, time-dependent.
- Analyses: static, transient (more than 20 explicit, implicit and higher-order time
  integrators, including the adaptive hybrid explicit method of Soares 2020), modal,
  harmonic frequency domain.
- Contour plots of displacements, strains and stresses.

## Quick start

```matlab
cd code
Start_Sim          % static Lee & Bathe cubic patch test on a distorted mesh, p = 3
```

The model is defined in an ABAQUS input file (`code/InputFilesAbaqus`). The parameters
(degree, analysis type, mass matrix, time integration, …) are set at the top of
`Start_Sim.m`. Requirements: MATLAB R2018b or newer; the Parallel Computing Toolbox is
optional (parallel element assembly).

## Documentation

The documentation is written in LaTeX with TikZ figures (`docs/`):

| Chapter | Content |
|---|---|
| Theory | Petrov–Galerkin weak form, test and trial spaces, skew coordinates, element matrices |
| Getting started | requirements, running a simulation, parameters, output |
| ABAQUS input files | supported keywords, conventions, data structures, Lee & Bathe example |
| Code reference | all routines by module, third-party code |
| Analyses | static, dynamic (time integrators), modal, frequency domain |
| Test suite | unit and integration tests, helpers, fixtures |
| Known issues | remaining limitations and fixed bugs |
| Literature | annotated bibliography (`literature/references.bib`) |

Build the PDF (TeX Live or MiKTeX with biber):

```bash
cd docs
latexmk -pdf usem-documentation.tex
```

## Tests

```matlab
cd tests
run_tests          % about 630 tests, roughly 2 minutes, base MATLAB only
```

The suite checks, among other things, the completeness and frame invariance of the metric
trial functions, rigid-body modes and rank of the unsymmetric element matrices, exact patch
tests on distorted meshes for $p = 1,\dots,4$, the Lee & Bathe benchmarks against the
analytical solution, the ABAQUS readers, and the accuracy and convergence order of the time
integrators. GitHub Actions runs it on every push.

## Studies

- `studies/thick_cylinder/`: distortion sensitivity and convergence of the unsymmetric
  spectral element for non-polynomial fields (thick cylinder of Rajendran 2010, Sec. 5.9,
  $p=2,\dots,8$, gmsh meshes, manufactured solutions). Report:
  `studies/thick_cylinder/report/report.pdf`.
- `studies/rajendran_panel/`: reproduction of the nodal rotation patch test and the
  convergence test of Rajendran (2010), Secs. 5.5-5.7 ($\theta\ge90^\circ$, folded elements),
  9-node and $p\le8$. Report: `studies/rajendran_panel/report/report.pdf`.

## Repository layout

```
code/
  Start_Sim.m              main script (settings → input → analysis)
  InputFilesAbaqus/        example ABAQUS input files (Lee & Bathe patch tests)
  subroutines/             element technology, assembly, BCs, solvers, plots, ...
docs/                      documentation (LaTeX, TikZ figures)
tests/                     test suite (matlab.unittest)
literature/references.bib  bibliography (the PDFs are not distributed)
```

## Key references

- S. Rajendran, K. M. Liew (2003). A novel unsymmetric 8-node plane element immune to mesh
  distortion under a quadratic displacement field. *IJNME* 58(11):1713–1748.
  [doi:10.1002/nme.836](https://doi.org/10.1002/nme.836)
- S. Rajendran (2010). A technique to develop mesh-distortion immune finite elements.
  *CMAME* 199:1044–1063. [doi:10.1016/j.cma.2009.11.017](https://doi.org/10.1016/j.cma.2009.11.017)
- Q. Xie, K. Y. Sze, Y. X. Zhou (2016). Modified and Trefftz unsymmetric finite element
  models. *IJMMD* 12(1):53–70. [doi:10.1007/s10999-014-9289-3](https://doi.org/10.1007/s10999-014-9289-3)
- N.-S. Lee, K.-J. Bathe (1993). Effects of element distortions on the performance of
  isoparametric elements. *IJNME* 36(20):3553–3576.
  [doi:10.1002/nme.1620362009](https://doi.org/10.1002/nme.1620362009)

## Licence

MIT, see [LICENSE](LICENSE). Third-party routines keep their own licences (chapter "Code
reference" of the documentation). The authors of the individual routines are given in the
file headers (S. Duczek / S. Eisenträger, Otto von Guericke University Magdeburg).
