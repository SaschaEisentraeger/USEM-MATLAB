# 5. Analysis types

`Start_Sim.m` selects the analysis with `analysis_type`. All four drivers
(`subroutines/AnalysisType/call_*.m`) start the same way:
- set up the $3\times3$ material matrix of every material: 2 constants give isotropic
  `HookeMatrix3d`, 9 constants give orthotropic;
- reduce it to plane stress or plane strain;
- assemble with `AssemblyFEM`.

## 5.1 Static analysis (`call_Static`)

1. Assemble $\mathbf K$ (`flag_comp = 0`).
2. Find the homogeneous Dirichlet DOFs.
3. Assemble the loads (`NeumannBCs`).
4. Apply the prescribed displacements (`DirichletBCs_nonhomogeneous`).
5. Solve $\mathbf K_{aa}\mathbf u_a = \mathbf F_a$ with `\` (sparse LU, since $\mathbf K$ is
   unsymmetric).
6. Save and plot (`Contour_plot`).

## 5.2 Dynamic analysis (`call_Dynamic`)

1. `AssemblyFEM(..., 2)`: $\mathbf K$, $\mathbf M$, Rayleigh $\mathbf R$ and the stable time
   step `dte` of the central difference method.
2. Assemble the loads with one column per load case if `flag_mRHS`. Prescribed displacements
   are handled by lifting (`DirichletBCs_nonhomogeneous_dynamic`).
3. Write the matrices to `subroutines/Temp/`, then call `Solver_Transient`. It reduces
   everything to the active DOFs, chooses the time step from
   `dt = [dtf, dte, TimeIntegration.dtImp]` and calls the integrator:
   - `dtf = 0.05/f_max` (20 samples per period of the highest excited frequency, `FFT_Signal`);
   - the simulation time is the time period of the ABAQUS `*Dynamic` step.
4. Only the nodes of interest are stored (`save_var = 0`).

**Important:** with `flag_mRHS = true` (the default) `Solver_Transient` ignores the method
code. It runs `CentralDifferenceMethod_mRHS` for codes below 300 and `HHTAlphaMethod_mRHS`
otherwise.

Time step selection inside the integrators:
- **Explicit:** `dt = min(dtmax)` after scaling `dte` by the stability factor of the method.
- **Implicit:** `dt = dtImp` if all three entries differ, otherwise `dtf`.
- **Both:** the step is then adjusted to `tsim/ceil(tsim/dt)`.

| Code | Routine | Remarks |
|---|---|---|
| 101, 601 | `CentralDifferenceMethod` | explicit, 2nd order |
| 102 | `NohBatheMethod` | explicit, $p = 0.54$ |
| 103 | `SoaresMethodExplicit` | |
| 104 | `MirbagheriMethod` | Noh–Bathe with the factor $\eta = 1.025$ on the acceleration term of the predictor (dissipation against spurious oscillations); first order by design |
| 105 | `KimMethodExplicit` | |
| 106, 107 | `ZhangXingMethod1`, `ZhangXingMethod2` | $\rho = 1$ and $\rho = 0$ |
| 108 | `ExplicitGeneralizedAlphaMethod` | $\rho_b = 1$ |
| 109 | `AdaptiveHybridExplicitMethod` | adaptive hybrid explicit method (Soares 2020), see [5.2.1](#521-adaptive-hybrid-explicit-method-109) |
| 201 | `KimExplicit_ho` | order 3 or 4 |
| 202 | `KimMethod2Stage` | two-stage explicit method |
| 301 | `NewmarkMethod` | trapezoidal rule $\beta = 1/4$, $\gamma = 1/2$; also the fallback for unknown codes |
| 302 | `HHTAlphaMethod` | $\alpha = -0.05$ |
| 303 | `GeneralizedAlphaMethod` | $\rho_\infty = 0.8$ |
| 304 | `BatheBaigMethod` | |
| 305 | `b1b2_BatheMethod` | $\beta_1 = 1/3$ |
| 306 | `rho_BatheMethod` | $\rho_\infty = 0$ |
| 307 | `SoaresMethodImplicit` | |
| 308, 309 | `KimMethodImplicit`, `KimMethodImplicitCollocation` | |
| 401 | `KimLagrange_ho` | order 1–5 |
| 402 | `SoaresMethodImplicit_ho` | |
| 403 | Padé-based integrator | order `Order` = 2–10; $\rho_\infty$ = `TimeIntegration.rho` (1: `PadeTimeIntegrator`, < 1: `PadeTimeIntegrator_ND`, or the Fortran solver if it is on the path) |
| 501 | `HPTIM_Newmark` | highly precise time integration |

Common calling convention:

```matlab
[T, U(, V, A)] = Method(K, M, C, F, F_dt, U0, V0, NumDoF, DOFsave, ActiveDoFs, ...
                        ExciteParam, SignalType, dtmax, tsim [, extra])
```

| Argument | Meaning |
|---|---|
| `K`, `M`, `C` | reduced (active-DOF) matrices |
| `F` | constant load (applied at $t = 0$) |
| `F_dt` | time-dependent load. Either a vector, or a struct with the fields `.Neumann`, `.K`, `.C` and `.M`, multiplied by $s(t)$, $s(t)$, $\dot s(t)$ and $\ddot s(t)$ respectively |
| `U0`, `V0` | initial displacement and velocity; the initial acceleration follows from the equation of motion |
| `NumDoF`, `DOFsave`, `ActiveDoFs` | total DOFs, stored (reduced) DOFs, map from reduced to global DOFs |
| `dtmax` | `[dtf, dte, dtImp]` |
| `tsim` | simulation time |
| `extra` | waitbar flag (central difference) or order (higher-order methods) |
| `T`, `U`, `V`, `A` | time vector and histories `[nT+1 x numel(DOFsave)]`. They are expanded to all DOFs if all active DOFs are stored; the `_mRHS` variants return 3-D arrays (time × DOF × load case) |

`tests/unit/TimeIntegrationTest.m` verifies the accuracy and the convergence order of the
schemes on a 2-DOF oscillator.

### 5.2.1 Adaptive hybrid explicit method (109)

Soares, D.: A novel time-marching formulation for wave propagation analysis: the adaptive
hybrid explicit method. *CMAME* 366 (2020) 113095. The implementation follows Table 1 of the
paper.

- **Elements.** Every element $e$ gets its maximum sampling frequency
  $\Omega_e = \omega_{\max,e}\,\Delta t$ from its own (unsymmetric) eigenproblem
  $\mathbf K_e\boldsymbol\phi = \omega^2\mathbf M_e\boldsymbol\phi$. Elements with
  $\Omega_e \le 2$ use the *regular* formulation (stability limit 2, one stiffness product
  per step); elements with $2 < \Omega_e \le 4$ use the *enhanced* formulation (stability
  limit 4, an additional correction of the displacements).
- **Time step.** $\Delta t = 2\,\Delta t_\text{CDM}$, twice the critical time step of the
  central difference method (`dte` from `AssemblyFEM`).
- **Adaptive dissipation.** Numerical damping is switched on only for elements in which a DOF
  oscillates, i.e. where its last two increments have opposite signs (Eq. 7). The parameters
  $\alpha_e$ (Eq. 8b, capped at 3) and $\beta_e$ (Eq. 12b) put the bifurcation of the
  amplification matrix at $\Omega_e$, so the stability limit is not reduced.
- **Damping.** Rayleigh damping is supported; the damping ratio of Eqs. (8b) and (12b) is
  computed from the mass-proportional part, $\zeta_e = \alpha_R\,\Delta t/2$.
- **Requirements.** The element matrices are needed: for method 109 `call_Dynamic` requests
  them from `AssemblyFEM` (6th output `ElemData`) and saves them to
  `subroutines/Temp/ElemData.mat`. Set `flag_mRHS = false`, otherwise the method code is
  ignored. Use lumped mass matrices (`mass_matrix = 2`, GLL nodal quadrature): then
  $\mathbf E = \mathbf M + \tfrac12\Delta t\,\mathbf C$ is diagonal and no linear solver is
  needed; consistent matrices work, but every step then solves with $\mathbf E$.

`tests/integration/AdaptiveHybridExplicitTest.m` checks second-order convergence, damping,
and stable, accurate results at twice the central-difference limit, where the central
difference method diverges.

## 5.3 Excitation signals

`ExcitationSignal(SignalType, ExciteParam, t)` returns a row vector, and
`ExcitationSignalDerivative(..., order)` its derivatives.

| `SignalType` | Signal | `ExciteParam` |
|---|---|---|
| 0 | none | `tE` |
| 1 | sine burst with Hann window: $\sin(2\pi f t)\sin^2(\pi f t/n)$ for $t < n/f$ | `[f, n]` |
| 2 | sine with Gaussian window | `f` |
| 3 | linear chirp | `[f0, f1, t1]` |
| 4 | harmonic $a\cos(2\pi f_0 t) + b\sin(2\pi f_1 t)$ | `[a, b, f0, f1]` |
| 5 | piecewise linear | table `[t, amplitude; ...; tE, 0]` |
| 6 | Ricker wavelet, delay $1/f$ | `f` |
| 7 | Gaussian impulses | `[a, b, c, d, n]` (amplitude, first peak, width, spacing, number) |

`RayleighDamping(omega, zeta, flag_vis)` returns $\alpha$ and $\beta$ such that
$\zeta(\omega) = \alpha/(2\omega) + \beta\omega/2$ matches the targets. It uses the closed
form for two frequencies and least squares otherwise.

## 5.4 Modal analysis (`call_Modal`)

1. Assemble $\mathbf K$ and $\mathbf M$, apply homogeneous constraints.
2. Undamped: `eigs(K_aa, M_aa, numEigen, EVP_Solver.type)`, giving $f = \sqrt{|\lambda|}/(2\pi)$
   in Hz. The full `eig` is used if `2*numEigen` exceeds the number of active DOFs.
3. With damping: state-space form of $\lambda^2\mathbf M\mathbf u + \lambda\mathbf R\mathbf u + \mathbf K\mathbf u = \mathbf 0$,
   giving $f = |\lambda|/(2\pi)$ in Hz and the displacement part of the eigenvectors.
4. The first mode is plotted.

Because $\mathbf K$ and $\mathbf M$ are unsymmetric on distorted meshes, eigenvalues can
become complex.

## 5.5 Frequency-domain analysis (`call_FrequencyDomain`)

Solves $(\mathbf K - \omega^2\mathbf M)\hat{\mathbf u} = \hat{\mathbf F}$ for
$\omega = 2\pi f_\text{Harmonic}$. Constant and time-dependent loads are treated as harmonic
loads with unit amplitude. Damping and prescribed displacements are not taken into account.
