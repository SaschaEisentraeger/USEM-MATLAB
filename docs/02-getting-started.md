# 2. Getting started

## 2.1 Requirements

| Product | Needed for |
|---|---|
| MATLAB R2018b or newer (developed and tested with R2026a) | everything |
| Parallel Computing Toolbox | `Start_Sim.m` opens a parallel pool unconditionally (lines 33–41). Without the toolbox, delete these lines; the `parfor` loops then run serially. |
| Symbolic Math Toolbox | `MMS_2D` (manufactured solutions), `HigherOrderDisplacementField.m`, Padé helpers |
| Signal Processing Toolbox | `FFT_Signal` with `SignalType = 5` (`findpeaks`) |
| Optimization Toolbox | `Global2LocalCoordinatesQuad` (`fsolve`), not used in the main flow |
| Deep Learning Toolbox | `combvec` in the Q8 branch of `distort_mesh` / `distort_single_element` |

The test suite needs only base MATLAB.

## 2.2 Running a simulation

1. Change the MATLAB current folder to `code/`. All paths in `Start_Sim.m` are relative to it.
2. Choose an input file in `Start_Sim.m` (`filename = ...`). The ABAQUS input files in
   [`code/InputFilesAbaqus`](../code/InputFilesAbaqus) are found through the MATLAB path:

   | File | Test |
   |---|---|
   | `PatchTestLeeQuadratic_{undistorted,angular,midedge}.inp` | Lee & Bathe quadratic patch test (pure bending) |
   | `PatchTestLeeCubic_{undistorted,angular,midedge}.inp` | Lee & Bathe cubic patch test |
   | `PatchTestLeeCubic_Rajendran_adj.inp` | cubic test on five distorted 8-node elements |

3. Set the parameters (next section) and run `Start_Sim`.

```matlab
cd code
Start_Sim
```

`Start_Sim` clears the workspace, sets `dbstop if error`, deletes `subroutines/Temp/*.mat`,
reads the input file, lifts the mesh to degree $p$ and then calls one of the analysis scripts
(`call_Static`, `call_Dynamic`, `call_Modal`, `call_FrequencyDomain`). They are scripts, not
functions, so all variables stay in the base workspace.

## 2.3 Parameters in `Start_Sim.m`

| Variable | Meaning |
|---|---|
| `filename` | ABAQUS input file (`.inp`, without path) |
| `analysis_type` | `'Static_analysis'`, `'Dynamic_analysis'`, `'Modal_analysis'`, `'FrequencyDomain_analysis'` |
| `saveWS` | save the workspace to `subroutines/Results/...` |
| `plot_mesh`, `plot_disp`, `plot_stress`, `plot_strain` | plots |
| `antype` | 1 = plane stress, 2 = plane strain |
| `p.pxi`, `p.peta` | polynomial degree in $\xi$ and $\eta$. If the local coordinate systems of neighbouring elements are not aligned, both are set to the minimum. |
| `p.ndist` | node distribution (1 equidistant, 2 GLL, 3 GLC) |
| `mass_matrix` | 1 consistent, 2 nodal quadrature (GLL), 3 row-sum, 4 HRZ |
| `flag_mRHS` | dynamic analysis: one right-hand side per load case. If true, only the central difference (`Method < 300`) or HHT-α (`Method >= 300`) solvers are used. |
| `omega_vec`, `zeta_vec` | Rayleigh damping: two angular frequencies [rad/s] and damping ratios |
| `numEigen`, `EVP_Solver.*` | modal analysis (`eigs` options) |
| `SignalType`, `ExciteParam` | time signal of time-dependent loads, see [chapter 5](05-analyses.md#53-excitation-signals) |
| `save_var` | dynamic analysis: 0 = save only the nodes of interest (DOF 6 in the input file), 1 = all DOFs |
| `TimeIntegration.Method` | time integration scheme, see [chapter 5](05-analyses.md#52-dynamic-analysis-call_dynamic) |
| `TimeIntegration.Order`, `.rho`, `.dtImp`, `.nInterp`, `.SolverSave` | options of the time integrators |
| `fHarmonic` | excitation frequency of the frequency-domain analysis [Hz] |

## 2.4 Output

| Analysis | Results file |
|---|---|
| static | `subroutines/Results/Static/<name>_Static/<name>_ndof<N>.mat` |
| dynamic | `subroutines/Results/Dynamic/<name>_<Method>/<name>_ndof<N>_<mass>.mat` |
| modal | `subroutines/Results/Modal/<name>_Modal/<name>_ndof<N>_<mass>.mat` |
| frequency domain | `subroutines/Results/Frequency_Domain/<name>_Frequency_Domain/...` |

These folders are excluded from version control.

The most important variables in a saved workspace:

| Variable | Content |
|---|---|
| `Topology_Node`, `xCoor`, `yCoor` | spectral element mesh |
| `K`, `M`, `R` | stiffness, mass, damping (sparse, after the Dirichlet conditions) |
| `F` | load vector |
| `u`, `ux`, `uy` | static displacements |
| `EF`, `U` | eigenfrequencies and mode shapes (modal analysis) |
| `t_comp`, `y_comp` | time and displacement history (dynamic analysis) |

## 2.5 Using the building blocks without `Start_Sim`

The test helpers show how to use the routines directly:

```matlab
addpath(genpath('code/subroutines')); addpath('tests/helpers');

% 3 x 2 mesh of 8-node elements with curved edges, lifted to p = 4
[nodes, topo] = usemTestMesh(3, 2, 3, 2, 2);
nodes = usemDistortMesh(nodes, topo, 'combined', 0.3);
model = usemBuildModel(nodes, topo, 4);          % GLL nodes

% unsymmetric stiffness matrix (plane stress, E = 1, nu = 0.3)
K = usemAssemble(model, HookeMatrix2d(1, 0.3, 1));

% displacement patch test with the field u = (Re z^4, -Im z^4)
field = @(x, y) deal(real((x + 1i*y).^4), -imag((x + 1i*y).^4));
[u, uExact, err] = usemDirichletSolve(K, model, field);   % err ~ 1e-15
```

`usemReadAbaqus(inpFile, p)` runs the input part of `Start_Sim.m` (reading and lifting to
degree $p$) and returns the model as a struct.
