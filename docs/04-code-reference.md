# 4. Code reference

All routines are in `code/subroutines`. `Start_Sim.m` adds this folder with all subfolders
to the path (`genpath`). Scripts (no `function` header) operate on the base workspace.

## 4.1 Program flow

```
Start_Sim.m                         settings, parallel pool
 ├─ Adjust_Set_names_input_file_ABAQUS, Read_*_ABAQUS     input   (chapter 3)
 ├─ pSpectralElement_NodesTopology                        mesh of degree p
 ├─ pSpectralElement_BCs                                  edge sets of degree p
 ├─ FE_Mesh                                               plot (optional)
 └─ call_Static | call_Dynamic | call_Modal | call_FrequencyDomain        (chapter 5)
       ├─ HookeMatrix3d → PlaneStress/PlaneStrainConstitutiveMatrices
       ├─ AssemblyFEM → Ke_2D | KMe_2D (parfor over elements)
       ├─ DirichletBCs_homogeneous / _nonhomogeneous(_dynamic)
       ├─ NeumannBCs → CreateForce, CreatePress, CreateTraction, CreateBload, CreateGravity
       ├─ solver (\, eigs, Solver_Transient → time integrator)
       └─ Contour_plot → Calc_results, Results_Stress_Strain
```

## 4.2 Element technology

| Routine | Purpose |
|---|---|
| `ElementMatrices/Ke_2D(p, C, x, y [, nIPxi, nIPeta])` | unsymmetric element stiffness $\int \tilde{\mathbf B}_L^{\mathsf T}\mathbf C\mathbf B_R$ (default $2p$ Gauss points per direction) |
| `ElementMatrices/Me_2D(p, rho, x, y, lumping [, nIPxi, nIPeta])` | unsymmetric element mass, lumping 1–4 |
| `ElementMatrices/KMe_2D(p, C, rho, x, y, lumping)` | stiffness and mass in one loop (used by the dynamic, modal and frequency-domain analyses) |
| `ShapeFunctions/2D/ShapeFunctions2D_SEM(xi, eta, p)` | parametric (test) shape functions and $\partial/\partial\xi$, $\partial/\partial\eta$; returns `[nNodes x nPoints]` arrays |
| `ShapeFunctions/2D/ShapeFunctions2D_Metric(xQ, yQ, xN, yN, p [, method])` | metric (trial) shape functions and $\partial/\partial x$, $\partial/\partial y$; `method` 1 = skew coordinates (default), 0 = global coordinates |
| `ShapeFunctions/2D/CorotationalSkewCoordinates(xQ, yQ, xN, yN, p)` | skew coordinates $(s,t)$ and $\mathbf J_0$ |
| `ShapeFunctions/2D/ShapeFunctions2D_4nodes`, `_8nodes` | bilinear and serendipity functions (geometry of the input mesh) |
| `ShapeFunctions/2D/ShapeFunctions2D_Metric_Tria` | unfinished metric functions for triangles ($p\le2$, not used) |
| `ShapeFunctions/1D/LagrangePolynomials(x)` | 1D Lagrange polynomials (barycentric form) as function handles |
| `ShapeFunctions/1D/ShapeFunctions1D_SEM(p, xi, ndist)` | 1D shape functions, order `[xi=-1, xi=+1, interior]` |
| `GeometryMapping/Jacobi_Matrix(p, [xi eta], x, y)` | $\mathbf J = [x_{,\xi}\ y_{,\xi};\ x_{,\eta}\ y_{,\eta}]$ |
| `GeometryMapping/Local2GlobalCoordinatesQuad(p, P, [x y])` | $\mathbf x(\xi,\eta)$ |
| `GeometryMapping/Global2LocalCoordinatesQuad` | inverse map with `fsolve` (Optimization Toolbox) |
| `Quadrature/lgwt(n, a, b)` | Gauss–Legendre points and weights (G. von Winckel) |
| `Quadrature/lglnodes(N)` | $N+1$ Gauss–Lobatto–Legendre points and weights (G. von Winckel) |
| `NodalDistribution/nodal_distribution(p, ndist)` | supporting points on $[-1,1]$ |
| `Lumping/RowSum_Method`, `HRZ_Method` | lumping of an element mass matrix |

## 4.3 Mesh and topology

| Routine | Purpose |
|---|---|
| `Topology/pSpectralElement_NodesTopology(T, x, y, p, p_org)` | creates the edge and interior nodes of degree $p$ from the 4-node ($p_\text{org}=1$) or 8-node ($p_\text{org}=2$) mesh; shared edges get the same nodes, with their order flipped for opposite orientation |
| `Topology/Topology_Edge_2D(T)` | global edge numbers `[nElem x 4]` (edges 1-2, 2-3, 3-4, 4-1) |
| `Topology/ConnectivityEdgeElement(TE)` | elements and local edge numbers of every edge |
| `Topology/ConnectivityNodeElement(T)` | elements of every node |
| `Topology/Boundary_Edges(TE)` | boundary edges `[element, local edge, edge]` |
| `Topology/Topology_DoF_2D(T)` | DOF connectivity, node $n$ → DOFs $2n-1$, $2n$ |
| `MeshDistortion/distort_mesh(X, T, 'Q4'|'Q8', ...)` | random distortion of interior nodes with a $\det\mathbf J$ safeguard (stand-alone) |
| `MeshDistortion/distort_single_element` | distortion of the element closest to the centre |

## 4.4 Material

Voigt notation `[11 22 33 12 23 13]` (engineering shear strains).

| Routine | Purpose |
|---|---|
| `HookeMatrix3d(E, nu)` | isotropic $6\times6$ stiffness |
| `HookeMatrix3d_ortho(E, G, nu)` | orthotropic $6\times6$ stiffness |
| `HookeMatrix3d_transiso(EL, ET, GLT, nuLT, nuTT)` | transversely isotropic (fibre direction 1) |
| `RotationHooke(C, phi, 'x'|'y'|'z')` | rotation of a $6\times6$ stiffness (angle in rad) |
| `PlaneStressConstitutiveMatrices(C)` | static condensation of $\sigma_{33}=0$, giving a $3\times3$ matrix |
| `PlaneStrainConstitutiveMatrices(C)` | in-plane sub-matrix, $3\times3$ |
| `HookeMatrix2d(E, nu, antype)` | closed-form $3\times3$ matrices (1 plane stress, 2 plane strain) |

## 4.5 Assembly, boundary conditions and loads

| Routine | Purpose |
|---|---|
| `Assembly/AssemblyFEM(T, p, x, y, nElem, nDoF, C, rho, nMat, mass_matrix, elsets, alpha, beta, zeta, flag_comp)` | sparse global $\mathbf K$, $\mathbf M$, $\mathbf R = \alpha\mathbf M + \beta\mathbf K$ (only for elements with damping). `flag_comp`: 0 = $\mathbf K$; 1 = $\mathbf K$, $\mathbf M$; 2 = $\mathbf K$, $\mathbf M$ and the critical time step $\Delta t = 0.99\cdot 2/\omega_{\max}$ from element eigenvalues. Optional 6th output `ElemData`: element matrices and DOF connectivity (adaptive hybrid explicit method) |
| `BoundaryConditions/Dirichlet/DirichletBCs_homogeneous(bc)` | constrained DOFs `[DOF_C, NOI]` (`fixed`, `UxUy 0`, `XSYMM`/`Ux 0`, `YSYMM`/`Uy 0`; `noi` → nodes of interest) |
| `BoundaryConditions/Dirichlet/DirichletBCs_nonhomogeneous(K, F, nDoF, bc, x, y)` | prescribed values by elimination: $\mathbf F \mathrel{-}= \mathbf K_{:,b}\bar{\mathbf u}_b$, rows and columns of $b$ set to zero, diagonal kept, $F_b = K_{bb}\bar u_b$ |
| `BoundaryConditions/Dirichlet/DirichletBCs_nonhomogeneous_dynamic(K, C, M, ...)` | lifting of constant and time-dependent prescribed displacements; returns load parts multiplied by $s(t)$, $\dot s(t)$, $\ddot s(t)$ |
| `BoundaryConditions/Dirichlet/SpatialDisplacementFunction` | polynomial equilibrium fields for higher-order patch tests ($p$ and $\nu = 0.3$ hard-coded) |
| `BoundaryConditions/Dirichlet/HigherOrderDisplacementField.m` | script (Symbolic Math Toolbox) that generates those fields |
| `BoundaryConditions/Neumann/NeumannBCs(p, cload, dsload, body_force, gravity, MAT, nDoF, T, x, y, nNodes [, flag_mRHS])` | load vector(s) `[F, F_dt]` (constant and time-dependent part); with `flag_mRHS` one column per load |
| `BoundaryConditions/Neumann/CreateForce`, `CreatePress`, `CreateTraction`, `CreateBload`, `CreateGravity` | individual loads (edge integrals with $p+1$ Gauss points, area integrals with $(p+1)^2$) |
| `BoundaryConditions/Neumann/SpatialLoadingFunction` | spatially varying edge loads (Lee & Bathe cubic patch test) |
| `BoundaryConditions/pSpectralElement_BCs(dsload, bc, T, p, p_org)` | adds the high-order edge nodes to the node sets of loads and boundary conditions that have an element set |

## 4.6 Analyses and time integration

See [chapter 5](05-analyses.md).

## 4.7 Post-processing and utilities

| Routine | Purpose |
|---|---|
| `Visualization/Contour_plot(p, x, y, U, scale, fields, T, C, elsets, flag_ave, alpha)` | contour plots on the deformed mesh (`'ux'`, `'uy'`, `'mag'`, `'Mises'`, `'sigma_x'`, `'sigma_y'`, `'tau_xy'`, `'epsilon_x'`, `'epsilon_y'`, `'gamma_xy'`); evaluation with the metric functions on $(20p+1)^2$ points per element |
| `Visualization/Contour_plot/Calc_results`, `Results_Stress_Strain`, `stress_averaging` | field evaluation and nodal averaging |
| `Visualization/Discretization/FE_Mesh` | mesh plot with loads and boundary conditions |
| `Visualization/Colormap` | diverging colour map (Moreland), `maxdistcolor` (S. Cobeldick, BSD licence) |
| `Misc/strain_energy(K, u)` | $\tfrac12\mathbf u^{\mathsf T}\mathbf K\mathbf u$ |
| `Misc/RotateFeMatrix`, `Misc/quadmesh` | utilities |
| `MMS_2D/` | method of manufactured solutions (Symbolic Math Toolbox): generates displacement and body-load functions for convergence studies; not connected to `Start_Sim` |

## 4.8 Third-party code

| Code | Author | Licence |
|---|---|---|
| `Quadrature/lgwt.m`, `Quadrature/lglnodes.m` | Greg von Winckel (MATLAB File Exchange) | BSD (File Exchange); no licence file in the repository |
| `Visualization/Colormap/MaxDistColor/*` | Stephen Cobeldick | BSD-3, `license.txt` |
| `Visualization/Colormap/Diverging/diverging_map.m` | after Kenneth Moreland, by Andy Stein | none given |
| `Misc/quadmesh.m` | Allan P. Engsig-Karup (based on MathWorks' `trimesh`) | none given |
| `TimeIntegration/HigherOrder/Implicit/Pade/MatrixExponential/expm2.m`, `expm3.m` | I. Kuprov (Spinach); Bader, Blanes & Casas (with MathWorks code) | see file headers |
