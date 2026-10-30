# USEM: unsymmetric spectral elements with hybrid trial space (paper code)

[![tests](https://github.com/SaschaEisentraeger/USEM-MATLAB/actions/workflows/tests.yml/badge.svg?branch=Minimal_Code_Paper)](https://github.com/SaschaEisentraeger/USEM-MATLAB/actions/workflows/tests.yml)

Minimal MATLAB code to reproduce all results of the paper *Unsymmetric spectral elements
with hybrid trial space* (`manuscript/`). This branch contains only what the studies of the
paper use. The general-purpose code (analysis driver, time integration, modal and
frequency-domain analyses, ABAQUS input with loads and boundary conditions, visualization)
is on the branch `main`.

The unsymmetric element combines two sets of shape functions:

- **parametric test functions:** Lagrange polynomials of degree $p$ on Gauss-Lobatto-Legendre
  nodes, which keep the mesh conforming;
- **metric trial functions:** Lagrange polynomials in corotational **skew coordinates**
  (trial space $Q_p(s,t)$), or the **hybrid trial space** $T_p$ of the paper: complete
  polynomials of degree $p$ in skew coordinates, complemented by Legendre products in the
  reference coordinates.

The result reproduces every displacement field of total degree $\le p$ exactly on
arbitrarily distorted meshes.

## Requirements

MATLAB R2026a was used for all results. The thick-cylinder study needs the Symbolic Math
Toolbox; all other studies and the tests need base MATLAB only. The run scripts need `bash`,
the post-processing Python 3 (standard library), and the manuscript a LaTeX distribution
with `pgfplots`.

## Reproducing the paper

```bash
cd manuscript/studies
bash run_calibration.sh   # Gauss rules of all studies (data/quadrature)
bash run_reruns.sh        # mesh-level studies, thick cylinder
bash run_element.sh       # element-level studies, panel, timing
bash run_postprocess.sh   # plot data and tables (data/, generated/)
```

Then build the paper:

```bash
cd manuscript
latexmk -pdf "Unsymmetric Spectral Elements with Hybrid Trial Space.tex"
```

Chapter "Reproducing the results of the paper" of the documentation lists which script
produces which result.

## Documentation

The documentation is written in LaTeX with TikZ figures (`docs/`):

| Chapter | Content |
|---|---|
| Theory | Petrov-Galerkin weak form, test and trial spaces, skew coordinates, stiffness matrix |
| Metric basis | bases of the metric Vandermonde matrix, hybrid trial space |
| Reproducing the paper | requirements, run scripts, scripts and results |
| Code reference | all routines by module |
| Test suite | unit and integration tests, helpers |
| Known limitations | open limitations of the code |
| Literature | annotated bibliography (`literature/references.bib`) |

```bash
cd docs
latexmk -pdf usem-documentation.tex
```

## Tests

```matlab
cd tests
run_tests          % base MATLAB only
```

The suite checks the completeness and frame invariance of the metric trial functions,
rigid-body modes and rank of the unsymmetric stiffness matrices, and exact patch tests on
distorted meshes for $p = 1,\dots,4$. GitHub Actions runs it on every push.

## Repository layout

```
code/subroutines/          element technology: shape functions, stiffness matrix, topology
tests/                     test suite (matlab.unittest) and helpers, also used by the studies
studies/thick_cylinder/    thick cylinder of Rajendran (2010): scripts, meshes, raw results
studies/rajendran_panel/   panel of Rajendran (2010): scripts and meshes
manuscript/                paper, study scripts (studies/), data (data/) and tables (generated/)
docs/                      documentation (LaTeX, TikZ figures)
literature/references.bib  bibliography of the documentation
```

## Key references

- S. Rajendran, K. M. Liew (2003). A novel unsymmetric 8-node plane element immune to mesh
  distortion under a quadratic displacement field. *IJNME* 58(11):1713–1748.
  [doi:10.1002/nme.836](https://doi.org/10.1002/nme.836)
- S. Rajendran (2010). A technique to develop mesh-distortion immune finite elements.
  *CMAME* 199:1044–1063. [doi:10.1016/j.cma.2009.11.017](https://doi.org/10.1016/j.cma.2009.11.017)
- Q. Xie, K. Y. Sze, Y. X. Zhou (2016). Modified and Trefftz unsymmetric finite element
  models. *IJMMD* 12(1):53–70. [doi:10.1007/s10999-014-9289-3](https://doi.org/10.1007/s10999-014-9289-3)

## Contact

Sascha Eisenträger (corresponding author), Otto von Guericke University Magdeburg,
[sascha.eisentraeger@ovgu.de](mailto:sascha.eisentraeger@ovgu.de)

## Licence

MIT, see [LICENSE](LICENSE). Third-party routines keep their own licences (chapter "Code
reference" of the documentation). The authors of the individual routines are given in the
file headers (S. Duczek / S. Eisenträger, Otto von Guericke University Magdeburg).
