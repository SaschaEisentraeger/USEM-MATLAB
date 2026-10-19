# USEM: unsymmetric spectral element method (2D elasticity)

MATLAB code for **high-order unsymmetric (Petrov–Galerkin) spectral elements** for linear plane
elasticity. It covers static, transient, modal and frequency-domain analysis.

The element combines two sets of shape functions:
- **parametric test functions:** Lagrange polynomials of degree $p$ on GLL nodes, which keep
  the mesh conforming;
- **metric trial functions:** Lagrange polynomials in corotational **skew coordinates**,
  which are complete in Cartesian coordinates and frame-invariant.

The result is a spectral element that reproduces every displacement field of total degree
$\le p$ **exactly on arbitrarily distorted meshes** (angular and curved-edge distortion). It
generalizes the unsymmetric 8-node element of Rajendran & Liew (2003) and Rajendran (2010) to
arbitrary polynomial degree.

## Features

- Quadrilateral spectral elements of arbitrary degree $p_\xi$, $p_\eta$ (GLL, equidistant or
  Chebyshev nodes), lifted from 4-node or 8-node ABAQUS meshes.
- Unsymmetric stiffness and mass matrices: consistent, GLL nodal quadrature (diagonal),
  row-sum and HRZ lumping.
- Isotropic and orthotropic materials; plane stress and plane strain.
- Loads: nodal forces, pressure and traction (also spatially varying), body forces and
  gravity. Dirichlet conditions: homogeneous, prescribed, time-dependent.
- Analyses: static, transient (more than 20 explicit, implicit and higher-order time
  integrators), modal, harmonic frequency domain.
- Post-processing: displacement, strain and stress contour plots.

## Quick start

```matlab
cd code
Start_Sim          % static Lee & Bathe cubic patch test on a distorted mesh, p = 3
```

The model is defined in an ABAQUS input file (`code/InputFilesAbaqus`). The parameters
(degree, analysis type, mass matrix, time integration, …) are set at the top of `Start_Sim.m`.
The Parallel Computing Toolbox is optional (parallel element assembly), see
[Getting started](docs/02-getting-started.md).

## Repository layout

```
code/
  Start_Sim.m              main script (settings → input → analysis)
  InputFilesAbaqus/        example ABAQUS input files (Lee & Bathe patch tests)
  subroutines/             element technology, assembly, BCs, solvers, plots, ...
docs/                      documentation (see below)
tests/                     test suite (matlab.unittest)
literature/                annotated bibliography (README.md, references.bib); PDFs are not distributed
```

## Documentation

1. [Theory: the unsymmetric spectral element method](docs/01-theory.md)
2. [Getting started](docs/02-getting-started.md)
3. [ABAQUS input files](docs/03-input-files.md)
4. [Code reference](docs/04-code-reference.md)
5. [Analysis types and time integration](docs/05-analyses.md)
6. [Test suite](docs/06-testing.md)
- [Known issues](docs/known-issues.md)
- [Literature](literature/README.md)

## Tests

```matlab
cd tests
run_tests             % about 620 tests, roughly 2 minutes, only base MATLAB required
```

Among other things, the suite checks the following:
- completeness and frame invariance of the metric trial functions;
- rigid-body modes and rank of the unsymmetric element matrices;
- exact patch tests on distorted meshes for $p = 1\dots4$, together with a control
  experiment showing that the classical element fails;
- the Lee & Bathe benchmarks against the analytical solution;
- reading of the ABAQUS input files;
- accuracy and convergence order of the time integrators.

Tests of open bugs are tagged `KnownIssue` and excluded by default (there are none at the
moment). Remaining limitations are listed in [docs/known-issues.md](docs/known-issues.md). GitHub Actions runs the suite on every push
(`.github/workflows/tests.yml`).

## Key references

- S. Rajendran, K. M. Liew (2003). A novel unsymmetric 8-node plane element immune to mesh
  distortion under a quadratic displacement field. *IJNME* 58(11):1713–1748.
  [doi:10.1002/nme.836](https://doi.org/10.1002/nme.836)
- S. Rajendran (2010). A technique to develop mesh-distortion immune finite elements.
  *CMAME* 199:1044–1063. [doi:10.1016/j.cma.2009.11.017](https://doi.org/10.1016/j.cma.2009.11.017)
- Q. Xie, K. Y. Sze, Y. X. Zhou (2016). Modified and Trefftz unsymmetric finite element models.
  *IJMMD* 12(1):53–70. [doi:10.1007/s10999-014-9289-3](https://doi.org/10.1007/s10999-014-9289-3)
- N.-S. Lee, K.-J. Bathe (1993). Effects of element distortions on the performance of
  isoparametric elements. *IJNME* 36(20):3553–3576.
  [doi:10.1002/nme.1620362009](https://doi.org/10.1002/nme.1620362009)

The full bibliography is in [literature/references.bib](literature/references.bib).

## Authors

See the file headers for the authors of the individual routines (S. Duczek / S. Eisenträger,
Otto von Guericke University Magdeburg). Third-party code and its licences are listed in
[docs/04-code-reference.md](docs/04-code-reference.md#48-third-party-code).

## License

MIT, see [LICENSE](LICENSE). Third-party routines keep their own licences, see
[docs/04-code-reference.md](docs/04-code-reference.md#48-third-party-code).
