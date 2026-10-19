# Known issues

The issues below were found while writing the documentation and the test suite. The
scientific code was **not** changed. Issues marked *(test)* have a test with the tag
`KnownIssue`. These tests fail until the issue is fixed and are excluded from the default
run:

```matlab
cd tests
run_tests('KnownIssues', true)   % include them
```

After fixing an issue, remove the `KnownIssue` tag from its test.

## Element technology

<a id="ki-1"></a>
**KI-1 `Global2LocalCoordinatesQuad` always fails** *(test: `GeometryMappingTest/globalToLocalInvertsTheMapping`)*
- `InverseMapping` computes `N'*x` with an already transposed `N`, which is a dimension
  error. The routine is not used by the main program.
- Fix: `F = [N*x - P(:,1), N*y - P(:,2)]`.

<a id="ki-2"></a>
**KI-2 Different mass-matrix quadrature in `KMe_2D` and `Me_2D`** *(test: `ElementMatricesTest/combinedRoutineMatchesSeparateMassOnDistortedElements`)*
- `KMe_2D` integrates the consistent mass with $p+1$ Gauss points per direction, `Me_2D` with
  $2p$. The dynamic, modal and frequency-domain analyses use `KMe_2D`.
- On distorted elements the metric functions are not polynomials of degree $p$ in
  $(\xi,\eta)$, so $p+1$ points under-integrate $\mathbf N^{\mathsf T}\mathbf M\det\mathbf J$.
  The entries differ by 1–10 %; the total mass is exact in both.
- Decide which rule is intended.

<a id="ki-3"></a>
**KI-3 Nodal-quadrature lumping with non-GLL nodes crashes** *(test: `ElementMatricesTest/nodalQuadratureFallsBackForNonGLLNodes`)*
- `lumping = 2` with `ndist = 1` or `3` is supposed to fall back to HRZ or row-sum lumping.
  The branch sets the new lumping flag but never defines the quadrature points
  (`Me_2D` line 85, `KMe_2D`), giving "Unrecognized function or variable 'xi'".
- In `Me_2D`, `pxi`/`peta` are also undefined when the quadrature orders are passed as
  arguments.

## Boundary conditions and loads

<a id="ki-4"></a>
**KI-4 Three-node sets are mistaken for "generate" triplets** *(test: `DirichletBCTest/threeNodeSetIsNotMistakenForGenerateTriplet`)*
- `DirichletBCs_homogeneous` interprets every node set with three entries `[a; b; c]` as
  `a:c:b` if that range ends at `b`. For example, `{2, 8, 3}` becomes `{2, 5, 8}`. Edge sets of
  $p = 2$ elements have three nodes, so wrong DOFs can be constrained.
- Fix: expand `generate` sets in the reader only (it already does that).

<a id="ki-5"></a>
**KI-5 Loads (`NeumannBCs`, `Create*`)**
- `NeumannBCs(..., flag_mRHS)` with `flag_mRHS = false` passed explicitly leaves `nLC`
  undefined.
- The body-force loop passes the whole `body_force` array to `CreateBload` instead of
  `body_force(i)`.
- The gravity loop tests `body_force(i).type` (error if there is no body force).
- `CreateBload` and `CreateGravity` index `elset(i)` with the element number, so they fail
  for element subsets.
- `CreateGravity` always uses the density of material 1.
- `SpatialLoadingFunction`: the constant `+6` in the traction of the cubic patch test is not
  scaled by the amplitude.

## Input files

<a id="ki-6"></a>
**KI-6 ABAQUS readers**
- `Read_boundary_condition_ABAQUS` loops forever when a boundary-condition set is not found.
  This happens, for example, when the sets are defined inside the part (before
  `*End Instance`) and are therefore not renamed to `PickedSet-`. The `Uy` branch
  (line ~192) can also loop without reading a new line.
- Boundary conditions are only counted if their data line contains `Set`.
- `Adjust_Set_names_input_file_ABAQUS` writes the file with `fprintf(fout, ContNew)`, which
  treats `%` and `\` in the input as format characters. Use `fprintf(fout, '%s', ContNew)`.
  It also uses a hard-coded `'\'` separator.
- `Read_mesh_ABAQUS` ignores node and element labels (they must be consecutive from 1).
  For 8-node meshes, `pSpectralElement_NodesTopology` additionally assumes that the corner
  nodes are numbered first.
- `Read_step_ABAQUS` overwrites the whole `tend` vector with 0 for non-dynamic steps.
- `Read_load_ABAQUS`: a `GRAV` line with three components gives NaN; `*Cload` sets do not
  support `generate`.
- The ABAQUS input files that `Start_Sim.m` refers to were missing. They were reconstructed
  in `code/InputFilesAbaqus` (see [chapter 3](03-input-files.md)).

## Time integration and analysis drivers

<a id="ki-7"></a>
**KI-7 `rho_BatheMethod` (code 306) cannot run** *(test: `TimeIntegrationTest/rhoBatheMethodRuns`)*
- The 5th argument is named `Fdt`, but the body uses `F_dt`.

<a id="ki-8"></a>
**KI-8 `Solver_Transient`**
- Line 47: `dof_poi = 1:size(ActiveDoF,2)` is always `1`, because `ActiveDoF` is a column. So
  `save_var = 1` stores a single DOF.
- Method 109: `FE_Matrices.mat` is never written, `Rayleigh` is undefined and the argument
  list does not match `AdaptiveHybridExplicitMethod`.
- Method 202: passes one argument too many to `KimMethod2Stage`, whose outputs are declared
  as `[U, T]`.
- Method 403: `PadeTimeIntegrator_ND_FORTRAN` is not in the repository, and `K + K'`,
  `M + M'`, `C + C'` double the full (not triangular) matrices while `F` is not doubled.
- With `flag_mRHS = true` (the default) the method code is ignored (see [chapter 5](05-analyses.md)).

<a id="ki-9"></a>
**KI-9 `MirbagheriMethod` (code 104) converges with first order** *(test: `TimeIntegrationTest/mirbagheriMethodIsSecondOrder`)*
- On the 2-DOF free-vibration problem, the error halves (instead of quartering) when $\Delta t$
  is halved. Compare with the published scheme.

<a id="ki-10"></a>
**KI-10 Modal analysis with damping and start vector**
- `call_Modal`: the damped (state-space) branch returns `EF` in rad/s instead of Hz, and
  assigns the $2n$-row eigenvectors to `U(ActiveDoF,:)` (dimension error).
- `EVP_Solver.v0 = true` creates a start vector of length `num_DoF` instead of
  `numel(ActiveDoF)`.

<a id="ki-11"></a>
**KI-11 Further observations**
- `CentralDifferenceMethod`: the stored velocity and acceleration of step $n+1$ are the
  values at $t_n$ (one step shifted).
- `call_Dynamic`: `NOI(1:3)` fails for fewer than three active nodes, and the post-processing
  block at the end ("For Jing") always runs.
- `Start_Sim` requires the Parallel Computing Toolbox (unconditional `parpool`).
- The result workspaces that came with the code (`subroutines/Results/Static`) were partly
  created with an earlier version of the load routines. The quadratic Lee & Bathe cases have
  the opposite load sign, and some $p\ge3$ runs had a uniform tension load. Only the cubic
  cases agree with the current code. The regression tests therefore check the analytical
  solution instead (`LeeBathePatchTest`).
- `distort_single_element` computes the characteristic length from the $x$ range twice.
- `Contour_plot` computes `ylim` from the $x$ extent.
