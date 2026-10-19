# 6. Test suite

The tests use the MATLAB unit test framework (`matlab.unittest`) and need only base MATLAB.

```matlab
cd tests
run_tests                       % unit + integration tests (about 2 minutes)
run_tests('unit')               % only tests/unit
run_tests('integration')        % only tests/integration
run_tests('KnownIssues', true)  % also the tests of open bugs (they fail)
run_tests('Coverage', true)     % Cobertura coverage report of code/subroutines
```

From the command line or in CI:

```bash
matlab -batch "cd tests; run_tests"
```

`run_tests` writes `tests/results/junit.xml` (and `coverage.xml`) and throws an error if a
test fails. It temporarily switches off the automatic creation of a parallel pool, so the
`parfor` loops run serially; pass `'Parallel', true` to allow a pool.

Every test class puts `code/subroutines` (with subfolders, as `Start_Sim.m` does) and
`tests/helpers` on the path with a `PathFixture`, so single classes can also be run with
`runtests`.

## 6.1 Unit tests (`tests/unit`)

| Class | What is checked |
|---|---|
| `QuadratureTest` | Gauss–Legendre and Gauss–Lobatto rules integrate polynomials of degree $2n-1$ exactly; interval mapping; symmetry |
| `NodalDistributionTest` | node distributions: ascending, symmetric, endpoints; GLL nodes are roots of $P_p'$; GLC nodes are Chebyshev extrema |
| `ShapeFunctionsTest` | 1D/2D Lagrange functions: Kronecker property, partition of unity, reproduction of polynomials and their derivatives, anisotropic degrees; bilinear and serendipity functions |
| `MetricShapeFunctionsTest` | metric functions: Kronecker property; **completeness in $(x,y)$ on distorted elements**; finite-difference derivatives; equality with the parametric functions on parallelograms; **frame invariance** (rotation, translation, scaling) and the frame dependence of `method = 0`; skew coordinates |
| `GeometryMappingTest` | Jacobian (constant for affine maps, finite differences, $\det\mathbf J>0$), node interpolation, element area |
| `MaterialTest` | Hooke matrices: symmetry, positive definiteness, Lamé constants, plane stress/strain reduction against closed forms, transverse isotropy and rotations |
| `ElementMatricesTest` | **rigid-body modes in the right and left null space**, rank $n-3$, symmetry on parallelograms (identical to the classical element), unsymmetry on distorted elements, constant-strain energy, rotation invariance $\mathbf K' = \mathbf T\mathbf K\mathbf T^{\mathsf T}$; mass: total mass, symmetry, diagonal GLL nodal quadrature, lumping |
| `TopologyTest` | DOF and edge connectivity, boundary edges, number and position of the high-order nodes on distorted 4- and 8-node meshes |
| `TimeIntegrationTest` | 18 integrators on a 2-DOF oscillator: free vibration, step load, **convergence order**, multiple right-hand sides |
| `ExcitationAndDampingTest` | excitation signals and their derivatives (finite differences), Rayleigh damping parameters |

## 6.2 Integration tests (`tests/integration`)

| Class | What is checked |
|---|---|
| `PatchTest` | displacement patch tests on regular, angular, curved and combined distorted meshes for $p = 1\dots4$, plane stress and plane strain: linear fields, all harmonic fields $(\operatorname{Re} z^k, -\operatorname{Im} z^k)$ with $k\le p$, and pure bending are reproduced **to machine precision**. Control experiment: the classical symmetric element fails on the same curved mesh. Also checks that fields of higher degree are not exact, and exponential p-convergence for $(e^x\cos y, -e^x\sin y)$. |
| `LeeBathePatchTest` | the Lee & Bathe patch tests (reference data in `tests/fixtures/reference`), repeated with the steps of `call_Static`: the quadratic test matches the analytical pure-bending solution for all meshes and $p\ge2$; the cubic test for $p = 3$ is mesh-independent and matches the stored solution |
| `AbaqusInputTest` | end to end: the input files in `code/InputFilesAbaqus` are read, lifted to degree $p$ (`pSpectralElement_*`) and compared with the reference models; static analysis from the input file |
| `DirichletBCTest` | homogeneous and prescribed displacements (elimination), uniaxial stretch on a distorted mesh |

## 6.3 Helpers and fixtures

| File | Purpose |
|---|---|
| `helpers/usemTestMesh` | structured 4- or 8-node mesh in ABAQUS node order |
| `helpers/usemDistortMesh` | deterministic angular, curved-edge or combined distortion of interior nodes |
| `helpers/usemTestElement` | single test elements (square, parallelogram, trapezoid, curved, general) |
| `helpers/usemBuildModel` | `pSpectralElement_NodesTopology` wrapper |
| `helpers/usemAssemble` | quiet `AssemblyFEM` wrapper |
| `helpers/usemDirichletSolve` | displacement-driven patch test |
| `helpers/usemKeSymmetric` | classical (Bubnov–Galerkin) element stiffness, for comparisons only |
| `helpers/usemReadAbaqus` | the input part of `Start_Sim.m` as a function |
| `fixtures/reference/*.mat` | model data and solutions of the Lee & Bathe tests, extracted from the result workspaces that came with the code (`make_reference_fixtures.m`) |
| `fixtures/make_abaqus_inputs.m` | reconstruction of the ABAQUS input files from the fixtures |

## 6.4 Writing new tests

- Put the new class in `tests/unit` or `tests/integration` and copy the `TestClassSetup`
  block (path fixtures) from an existing class.
- Use `evalc` around routines that print progress (`AssemblyFEM`, `NeumannBCs`,
  integrators), or the helpers that already do this.
- Use deterministic meshes (`usemDistortMesh`) or a seeded `RandStream`.
- Tag tests that document an unfixed bug with `TestTags = {'KnownIssue'}` and describe the
  bug in [known-issues.md](known-issues.md).
