# Known issues and fixed bugs

Tests that document an open bug get the tag `KnownIssue`. They are excluded from the default
run and fail until the bug is fixed:

```matlab
cd tests
run_tests('KnownIssues', true)   % include them
```

There are currently **no open `KnownIssue` tests**. The bugs found while writing the test
suite are fixed (see below). The remaining limitations are listed first.

## Remaining limitations

- **Time integration method 109** (`AdaptiveHybridExplicitMethod`, Soares 2020) is not
  available. The routine works element by element and expects element matrices, a cell-array
  DOF topology and triangular global matrices, i.e. it was written for another code.
  `Solver_Transient` stops with a clear error message for this method.
- **Default route of the dynamic analysis.** With `flag_mRHS = true`, `Solver_Transient` ignores
  the method code. It runs the central difference method (codes below 300) or HHT-α (300 and
  above). See [chapter 5](05-analyses.md).
- **ABAQUS readers** expect the structure written by ABAQUS/CAE:
  - boundary-condition data lines must contain `Set`;
  - sets must be defined in the assembly;
  - node labels must be consecutive from 1;
  - in 8-node meshes the corner nodes must be numbered first.

  Violations of the last three now raise errors instead of producing wrong data or endless
  loops. See [chapter 3](03-input-files.md).
- **Modal and frequency-domain analyses** only honour homogeneous Dirichlet conditions. The
  frequency-domain analysis ignores damping.
- `CentralDifferenceMethod_te`, `Solver_Transient_symrcm` and the second variants of the Padé
  integrators (`PadeTimeIntegrator2`, `PadeTimeIntegrator_ND2`) are not called by the program
  and are not tested.
- `ShapeFunctions2D_Metric_Tria` (triangles) is unfinished ($p\le2$) and not used.
- The result workspaces that came with the code (`subroutines/Results/Static`, not under version
  control) were partly created with an earlier version of the load routines (opposite load sign
  in the quadratic Lee & Bathe cases). The regression tests therefore use the analytical
  solution (`LeeBathePatchTest`).
- **Unsymmetric matrices in dynamics.** Eigenvalues can be complex, and energy is not
  conserved; see [chapter 1](01-theory.md#16-practical-implications-of-unsymmetric-matrices).

## Fixed bugs

| Area | Bug | Fix | Test |
|---|---|---|---|
| Nodal distributions | Mock-Chebyshev nodes (`ndist = 4`) were reversed, which gave indefinite stiffness matrices for $p\ge3$ | Mock and "fake" node code removed; `nodal_distribution` rejects unknown distributions | `NodalDistributionTest` |
| Path | Copies of core functions in `NodalDistribution/FakeNodes/src` shadowed the originals | Folder removed | – |
| Geometry | `Global2LocalCoordinatesQuad` failed on every call (`InverseMapping` used `N'` instead of `N`) | Corrected | `GeometryMappingTest/globalToLocalInvertsTheMapping` |
| Mass matrix | `KMe_2D` (dynamic, modal, frequency analyses) integrated the unsymmetric mass with $p+1$ Gauss points, `Me_2D` with $2p$; on distorted elements the entries differed by 1–10 % | `KMe_2D` uses $2p$ points | `ElementMatricesTest/combinedRoutineMatchesSeparateMassOnDistortedElements` |
| Lumping | `lumping = 2` with equidistant or GLC nodes crashed instead of falling back to HRZ / row-sum lumping | Quadrature selection restructured (`Me_2D`, `KMe_2D`) | `ElementMatricesTest/nodalQuadratureFallsBackForNonGLLNodes` |
| Sets | Every node, element or material set with three entries was taken for an ABAQUS `generate` triplet (`DirichletBCs_homogeneous`, `NeumannBCs`, `pSpectralElement_BCs`, `AssemblyFEM`, `Contour_plot`); e.g. the $p = 2$ edge set {2, 8, 3} became {2, 5, 8}, and a material on three elements got the wrong elements | Sets are used as given; the readers expand `generate` (now also for `*Cload`) | `DirichletBCTest`, `AbaqusReaderTest/concentratedLoadOnGeneratedSet` |
| Loads | `NeumannBCs(..., false)` crashed (`nLC` undefined); the body-force loop passed all body forces at once; the gravity loop tested `body_force(i)` | Corrected | `LoadsTest` |
| Loads | `CreateBload`/`CreateGravity` failed for element subsets; `CreateGravity` used the density of material 1 everywhere | Loops run over the element set; density per material | `LoadsTest` |
| Loads | `SpatialLoadingFunction`: constant part of the cubic-patch-test traction not scaled by the amplitude | Corrected | `LoadsTest/spatialLoadIsScaledByAmplitude` |
| Readers | `Read_boundary_condition_ABAQUS` looped forever if a set was missing; the U2 branch never read the next line | Error message; next line read | `AbaqusReaderTest` |
| Readers | Set names matched by prefix (`Set-1` also matched `Set-10`) | Exact matching (`isSetHeader`) in all readers | `AbaqusReaderTest/setNamesAreMatchedExactly` |
| Readers | `Adjust_Set_names_input_file_ABAQUS` used the file content as format string (`%`, `\`) and a hard-coded `\` separator | `fprintf(fout,'%s',...)`, `fullfile` | `AbaqusReaderTest/percentSignsAreCopiedLiterally` |
| Readers | Node labels ignored; corner nodes of 8-node meshes silently assumed first | Both checked, with error messages | `AbaqusReaderTest`, `TopologyTest/cornerNodesMustBeNumberedFirst` |
| Readers | `GRAV` with three direction components gave NaN; `Read_step_ABAQUS` overwrote all step times with 0 | Corrected | `AbaqusReaderTest/gravityWithThreeComponents` |
| Time integration | `rho_BatheMethod` (306) could not run (`Fdt` vs. `F_dt`) | Corrected | `TimeIntegrationTest`, `SolverTransientTest` |
| Time integration | `KimMethod2Stage` (202): wrong call, outputs `[U,T]`, `Minv*` instead of `Minv.*`, time step not adjusted to the simulation time | Corrected | `TimeIntegrationTest`, `SolverTransientTest` |
| Time integration | Padé method (403): Fortran solver missing; `K+K'` doubled the full matrices; undefined variable with non-zero initial conditions; damped solvers called `TimeIntgCoeff` with the Padé code instead of the order | MATLAB solvers used if the Fortran solver is missing; all errors corrected (verified convergence orders) | `SolverTransientTest` |
| Time integration | `Solver_Transient`: `save_var = 1` stored only one DOF | `1:numel(ActiveDoF)` | `SolverTransientTest/allDegreesOfFreedomAreSavedOnRequest` |
| Time integration | `CentralDifferenceMethod` (and `_mRHS`): velocity and acceleration stored one step late | Stored at the correct time; final values extrapolated | `TimeIntegrationTest/centralDifferenceVelocityAndAcceleration` |
| Modal analysis | Damped case: eigenfrequencies in rad/s instead of Hz, state-space eigenvectors assigned to the mode shapes (dimension error); random start vector of wrong length | Corrected | – |
| Dynamic analysis | `NOI(1:3)` failed for fewer than three nodes; the export block at the end ("For Jing") always ran | Corrected; export only when it is possible | – |
| `Start_Sim` | Required the Parallel Computing Toolbox | Pool only opened if the toolbox is available | – |
| Mesh distortion | `distort_single_element` used the $x$ range twice for the characteristic length | Corrected | – |

`MirbagheriMethod` (104) converges with first order. This is by design, not a bug: the factor
$\eta = 1.025$ on the acceleration term of the Noh–Bathe predictor adds dissipation against
spurious oscillations, and any $\eta \ne 1$ makes the scheme first-order accurate.
