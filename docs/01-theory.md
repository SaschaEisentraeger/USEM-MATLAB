# 1. Theory: the unsymmetric spectral element method

This chapter describes the formulation **as it is implemented** in `code/subroutines`. The
references are listed in [`literature/README.md`](../literature/README.md).

## 1.1 Problem

The code solves linear plane elasticity (plane stress or plane strain) in a domain $\Omega$:

$$
\nabla\cdot\boldsymbol\sigma + \mathbf b = \rho\,\ddot{\mathbf u},\qquad
\boldsymbol\sigma = \mathbf C\,\boldsymbol\varepsilon(\mathbf u),\qquad
\boldsymbol\varepsilon = [\varepsilon_{x},\ \varepsilon_{y},\ \gamma_{xy}]^{\mathsf T},
$$

with $\mathbf u = \bar{\mathbf u}$ on $\Gamma_u$ and $\boldsymbol\sigma\mathbf n = \bar{\mathbf t}$
on $\Gamma_t$. Here $\mathbf C$ is the $3\times3$ plane stress or plane strain matrix, which
is obtained from the 3D Hooke matrix (`PlaneStressConstitutiveMatrices`,
`PlaneStrainConstitutiveMatrices`). The thickness is 1.

## 1.2 Petrov–Galerkin weak form

The classical (Bubnov–Galerkin) finite or spectral element method uses the same
isoparametric shape functions for the test and the trial space. The unsymmetric method
uses two different spaces (Rajendran 2010):

| Space | Functions | Property |
|---|---|---|
| Test $\mathbf w_h \in \mathcal W_h$ | **parametric** shape functions $N_i(\xi,\eta)$ | inter-element compatible ($C^0$-conforming) |
| Trial $\mathbf u_h \in \mathcal U_h$ | **metric** shape functions $M_j(x,y)$ | complete in Cartesian coordinates |

The discrete problem is to find $\mathbf u_h\in\mathcal U_h$ such that for all
$\mathbf w_h\in\mathcal W_h$ (with $\mathbf w_h = \mathbf 0$ on $\Gamma_u$)

$$
\int_\Omega \boldsymbol\varepsilon(\mathbf w_h)^{\mathsf T}\mathbf C\,\boldsymbol\varepsilon(\mathbf u_h)\,\mathrm d\Omega
+ \int_\Omega \rho\,\mathbf w_h^{\mathsf T}\ddot{\mathbf u}_h\,\mathrm d\Omega
= \int_\Omega \mathbf w_h^{\mathsf T}\mathbf b\,\mathrm d\Omega
+ \int_{\Gamma_t} \mathbf w_h^{\mathsf T}\bar{\mathbf t}\,\mathrm d\Gamma .
$$

Both sets of shape functions are interpolatory at the same nodes, so the nodal displacements
are the unknowns in both cases. The resulting element matrices are **unsymmetric**.

### Exactness on distorted meshes

Let the exact solution $\mathbf u$ be a polynomial of total degree $\le p$ in $(x,y)$. The
metric interpolant of $\mathbf u$ is then $\mathbf u$ itself on every element, so the strains
and stresses are exact. Because the test functions are conforming, integration by parts
gives

$$\int_\Omega \boldsymbol\varepsilon(\mathbf w_h)^{\mathsf T}\boldsymbol\sigma(\mathbf u)\,\mathrm d\Omega
= \int_\Omega \mathbf w_h^{\mathsf T}\mathbf b\,\mathrm d\Omega
+ \int_{\Gamma_t}\mathbf w_h^{\mathsf T}\bar{\mathbf t}\,\mathrm d\Gamma .$$

The exact nodal values therefore satisfy the discrete equations, and the discrete solution
is exact, **independently of the element shape**. The theorem only holds if the quadrature is
exact. This is the distortion immunity of Rajendran's unsymmetric elements, here for
arbitrary $p$. The classical isoparametric element loses it as soon as the geometry map is
non-affine and $x^p\notin\mathcal W_h$: for example, quadratic fields on elements with curved
edges (Lee & Bathe 1993).

The metric functions need not be continuous across element edges; compatibility is provided
by the test functions.

## 1.3 Parametric (test) shape functions

`ShapeFunctions2D_SEM(xi, eta, p)` builds tensor products of 1D Lagrange polynomials through
$p+1$ nodes. Node distributions (`p.ndist`, `nodal_distribution.m`):

| `ndist` | Nodes |
|---|---|
| 1 | equidistant |
| 2 | Gauss–Lobatto–Legendre (GLL), the **default** (spectral elements) |
| 3 | Gauss–Lobatto–Chebyshev |

The geometry is isoparametric: the element nodes are generated from the 4-node (bilinear) or
8-node (serendipity) input mesh by `pSpectralElement_NodesTopology`. The geometry of an
element is therefore at most quadratic, whatever $p$ is.

### Local node numbering

For degree $p$, the element has $(p+1)^2$ nodes in this order: 4 corners, then the interior
nodes of edge 1 ($\eta=-1$), edge 2 ($\xi=+1$), edge 3 ($\eta=+1$) and edge 4 ($\xi=-1$), each
in the direction of increasing $\xi$ or $\eta$, then the interior nodes row by row ($\xi$
fastest). Example for $p = 3$:

```
 eta = +1   4 ---- 9 ---- 10 ---- 3
            |                     |
           12     15      16      8
            |                     |
           11     13      14      7
            |                     |
 eta = -1   1 ---- 5 ----- 6 ---- 2
          xi = -1               xi = +1
```

Global degrees of freedom are interleaved: node $n$ has $u_x$ at DOF $2n-1$ and $u_y$ at
DOF $2n$ (`Topology_DoF_2D`).

## 1.4 Metric (trial) shape functions in skew coordinates

`ShapeFunctions2D_Metric(xQ, yQ, xN, yN, p)` builds functions that are interpolatory at the
element nodes $(x_k, y_k)$ in physical space:

1. **Skew coordinates** (`CorotationalSkewCoordinates`, Xie, Sze & Zhou 2016). With the element
   centre $\mathbf x_0 = \mathbf x(0,0)$ and the Jacobian there,
   $\mathbf J_0 = \partial(x,y)/\partial(\xi,\eta)|_{(0,0)}$:

   $$\begin{bmatrix} s\\ t\end{bmatrix} = \mathbf J_0^{-\mathsf T}\,(\mathbf x - \mathbf x_0).$$

   $(s,t)$ is an affine, non-dimensional coordinate system that moves and rotates with the
   element. For a parallelogram it coincides with $(\xi,\eta)$.
2. **Basis.** $\phi_k(s,t) = s^a t^b$ with $0\le a\le p_\xi$ and $0\le b\le p_\eta$ (tensor-product
   space $Q_p$).
3. **Interpolation.** With the generalized Vandermonde matrix $V_{kl} = \phi_k(s_l,t_l)$ at the
   nodes, $\mathbf M(\mathbf x) = \mathbf V^{-1}\boldsymbol\phi(s(\mathbf x),t(\mathbf x))$, so
   $M_j(\mathbf x_k)=\delta_{jk}$.
4. **Derivatives.** $\nabla_{\mathbf x} M = \mathbf J_0^{-1}\nabla_{(s,t)}M$.

Properties (all covered by `tests/unit/MetricShapeFunctionsTest.m`):

- **Completeness.** $(s,t)$ is affine in $(x,y)$, so $Q_p(s,t)\supseteq P_p(x,y)$. Every
  polynomial of total degree $\le p$ is reproduced exactly on any element shape.
- **Frame invariance.** Rigid motions and uniform scaling of the element do not change the
  function values. The older variant `method = 0` (monomials in centroid-shifted global
  coordinates, Rajendran & Liew 2003) is complete but rotation-dependent for $p\ge2$
  (Ooi et al. 2008).
- **Affine elements.** On a parallelogram, metric = parametric, so the element is identical
  to the classical spectral element and $\mathbf K_e$ is symmetric.

## 1.5 Element matrices

With $J = \partial(x,y)/\partial(\xi,\eta)$ in the convention
$\mathbf J = \begin{bmatrix} x_{,\xi} & y_{,\xi}\\ x_{,\eta} & y_{,\eta}\end{bmatrix}$
(`Jacobi_Matrix`):

**Stiffness** (`Ke_2D`):

$$\mathbf K_e = \int_{\hat\Omega}\mathbf B_L^{\mathsf T}\mathbf C\,\mathbf B_R\,\det\mathbf J\,
\mathrm d\xi\,\mathrm d\eta
= \sum_g \tilde{\mathbf B}_L^{\mathsf T}(\xi_g)\,\mathbf C\,\mathbf B_R(\mathbf x_g)\,w_g,
\qquad
\tilde{\mathbf B}_L = \det\mathbf J\,\mathbf B_L .$$

$\tilde{\mathbf B}_L$ is built with the **adjugate** $\operatorname{adj}\mathbf J = \det\mathbf J\,\mathbf J^{-1}$,
so $\det\mathbf J$ never appears in a denominator. This is Rajendran's third condition for
distortion immunity. The rule assumes $\det\mathbf J > 0$ (no inverted elements).
$\mathbf B_R$ contains the derivatives of the metric functions evaluated at the images
$\mathbf x_g = \mathbf x(\xi_g,\eta_g)$ of the Gauss points.
Quadrature: $2p_\xi \times 2p_\eta$ Gauss–Legendre points.

Consequences:

- Rigid-body modes lie in both the right and the left null space of $\mathbf K_e$, and
  $\operatorname{rank}\mathbf K_e = n_\text{dof}-3$.
- $\mathbf K_e$ is unsymmetric unless the element is a parallelogram.

**Mass** (`Me_2D`, and `KMe_2D` together with the stiffness):

$$\mathbf M_e = \int_{\hat\Omega}\rho\,\mathbf N^{\mathsf T}\mathbf M\,\det\mathbf J\,\mathrm d\xi\,\mathrm d\eta ,$$

where $\mathbf N$ holds the parametric and $\mathbf M$ the metric functions.

| `mass_matrix` / `lumping` | Mass matrix |
|---|---|
| 1 | consistent (unsymmetric) |
| 2 | nodal quadrature with GLL points (diagonal, because both $N_i$ and $M_j$ are interpolatory at the nodes); requires `ndist = 2` |
| 3 | row-sum lumping |
| 4 | HRZ lumping (diagonal scaling) |

`Me_2D` uses $2p$ Gauss points per direction, `KMe_2D` $p+1$ (see
[known issues](known-issues.md#ki-2)).

**Loads.** Edge loads (pressure, traction) and body loads are integrated with the
*parametric* (test) functions (`CreatePress`, `CreateTraction`, `CreateBload`,
`CreateGravity`), which makes them consistent with the weak form above.

**Post-processing.** Displacements, strains and stresses inside the elements are evaluated
with the *metric* (trial) functions (`Calc_results`, `Results_Stress_Strain`).

## 1.6 Practical implications of unsymmetric matrices

- **Linear solvers.** $\mathbf K$ is unsymmetric, so use LU, i.e. MATLAB's `\` or
  `decomposition`, and not Cholesky.
- **Eigenvalues.** The generalized eigenproblem $\mathbf K\boldsymbol\phi = \omega^2\mathbf M\boldsymbol\phi$
  can have complex eigenvalues (Ooi 2006).
- **Time integration.** Energy is not conserved, and schemes that are unconditionally stable
  for symmetric systems (e.g. trapezoidal Newmark) can become unstable (Wang & Hillman 2022;
  Zähringer & Betsch 2026). Check dynamic results critically.
