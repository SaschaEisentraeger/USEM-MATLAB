# 3. ABAQUS input files

The model is defined in an ABAQUS/CAE input file (`*.inp`). The readers in
`subroutines/ReadInputFileAbaqus` are not a general ABAQUS parser. They expect the structure
that ABAQUS/CAE writes, including several of its comment lines. The safest approach is to
build the model in ABAQUS/CAE (a single part, a single instance, 2D planar, CPS4/CPE4 or
CPS8/CPE8 elements) and to start from one of the examples in
[`code/InputFilesAbaqus`](../code/InputFilesAbaqus).

## 3.1 Reading sequence (`Start_Sim.m`)

1. `Adjust_Set_names_input_file_ABAQUS(filename)` writes `<name>_adjusted.inp` to the current
   folder. In it:
   - every `Set-` after `End Instance` is renamed to `PickedSet-`;
   - the node/element set pairs named `Pressure` and `Traction` are numbered
     (`Pressure1`, `Pressure2`, …) in order of appearance.
2. `Read_mesh_ABAQUS`: nodes (x, y) and element connectivity.
3. `Read_material_ABAQUS`: materials and their element sets.
4. `Read_step_ABAQUS`: step names, types and the time period of `*Dynamic` steps.
5. `Read_load_ABAQUS`: concentrated forces, edge loads, body forces and gravity.
6. `Read_boundary_condition_ABAQUS`: Dirichlet conditions and nodes of interest.
7. The adjusted file is deleted.
8. `pSpectralElement_NodesTopology` and `pSpectralElement_BCs` lift the mesh and the
   edge-based sets to degree $p$.

## 3.2 Requirements and conventions

**Mesh**
- Only the first `*Node` block is read, and only x and y. Node labels must be numbered
  consecutively from 1 (checked by the reader).
- Element types ending in `4`/`4R` give 4-node elements, and types containing `8` give 8-node
  elements. One element type per model.
- For 8-node meshes the **corner nodes must be numbered before the mid-side nodes**
  (`pSpectralElement_NodesTopology` keeps nodes `1..nCorner` and stops with an error
  otherwise).
- The element block must be followed directly by `*Nset`, `*End Part` or `*End Instance`.

**Sections and materials**
- `*Solid Section, elset=<set>, material=<mat>`, preceded by a `** Section:` comment.
- A `** MATERIALS` header; the material block ends with the `** ----…` line that CAE writes.
- `*Density`, `*Elastic` (isotropic `E, nu`, `ENGINEERING CONSTANTS` or `ORTHO`) and `*Damping`.
  `*Damping` only switches Rayleigh damping on for the element set; the parameters come from
  `omega_vec` and `zeta_vec`.

**Sets**
- Sets used by boundary conditions and loads must be defined in the **assembly** (after
  `*End Instance`), as CAE does. Otherwise the boundary condition reader stops with the error
  "Node set ... was not found". Set names are matched exactly (`Set-1` is not `Set-10`).
- Node and element sets with `generate` are expanded by the readers (also for `*Cload`).
- Boundary-condition set names must contain `Set` (`Set-1`, `_PickedSet3`, …).
- A set on an **edge** must have a node set *and* an element set with the same name (CAE
  creates both when edges are picked). Only then can `pSpectralElement_BCs` add the
  high-order nodes of that edge. Point constraints only need a node set.
- Edge loads: define the node set and element set as `Pressure` (pressure loads) or
  `Traction` (surface tractions), in the same order as the loads in the step.

**Step, boundary conditions and loads** (only data after `** STEP:` is read)

| ABAQUS | Internal `type` | Notes |
|---|---|---|
| `*Boundary` `Set-1, ENCASTRE` / `PINNED` | `fixed` | $u_x = u_y = 0$ |
| `Set-1, XSYMM` / `YSYMM` | `XSYMM` / `YSYMM` | |
| `Set-1, 1, 1` | `Ux 0` | homogeneous |
| `Set-1, 1, 1, 0.001` | `Ux` (amplitude 0.001) | prescribed, applied by elimination |
| `Set-1, 1, 1` and `Set-1, 2, 2` (two lines, as CAE writes U1 + U2) | `UxUy 0` | |
| `*Boundary, amplitude=…` | `Ux dt`, `Uy dt` | time-dependent; time function from `SignalType` |
| `*Boundary, user` | `Ux fx`, `Uy fx` | spatially varying, `SpatialDisplacementFunction.m` |
| `Set-5, 6, 6` | `noi` | **node of interest** for the output of the dynamic analysis; not a constraint |
| `*Cload` `Set, 1, F` | `f_x` / `f_y` | the full value at every node of the set; `amplitude=` gives `f_x_dt` |
| `*Dsload` `Pressure, P, p` | `Pressure` | along the inward normal |
| `*Dsload` `Pressure, PNU, a` | `Pressure fx` | $p(x,y)$ from `SpatialLoadingFunction.m` |
| `*Dsload` `Traction, TRVEC, a, dx, dy, dz` | `Traction` | traction $a\,(d_x, d_y)$ |
| `*Dsload` `Traction, TRVECNU, a, dx, dy, dz` | `Traction fx` | $a\,(d_x, d_y)\,f(x,y)$ from `SpatialLoadingFunction.m` |
| `*Dload` with `BX`/`BY` (CAE comment `Body force`) | `Body_force` | |
| `*Dload` … `GRAV, g, dx, dy` | gravity | all elements |

`*Amplitude` tables are ignored. The time history of every time-dependent load is the signal
selected by `SignalType` in `Start_Sim.m`.

`SpatialLoadingFunction.m` contains the loads of the Lee & Bathe cubic patch test:
- `Pressure fx` gives $a \cdot 120\,y$;
- `Traction fx` gives $-a \cdot 6y^2 + 6$.

## 3.3 Data structures returned by the readers

| Variable | Layout |
|---|---|
| `nodes` | `[nNodes x 2]` coordinates |
| `topology` | `[nElem x 4]` or `[nElem x 8]` connectivity |
| `material_properties{i,1:5}` | density, dielectric (legacy), elastic constants (`[E nu]` or 9 orthotropic stiffness values), piezoelectric (legacy), damping |
| `material_elset{i}` | elements of material `i` |
| `cload(k)`, `dsload(k)` | struct arrays with the fields `set_name`, `amplitude`, `type`, `nset`, `elset`, `DirectionVector`, `step_name` (`cload` has no `elset` and no `DirectionVector`) |
| `body_force(k)`, `gravity(k)` | `amplitude = [bx; by]` or `[g; dx; dy]`, `elset` |
| `boundary_conditions(k)` | column struct array with the fields `set_name`, `amplitude`, `type`, `nset`, `elset`, `step_name`, **in reverse order of appearance** |

## 3.4 Example: Lee & Bathe quadratic patch test

The file `PatchTestLeeQuadratic_angular.inp` contains:
- a plate $[0,10]\times[-1,1]$ with two 4-node elements and angular distortion;
- $u_x = 0$ on the left edge (`Set-1`, node set *and* element set);
- $u_y = 0$ at $(0,-1)$ (`Set-2`);
- the pressure $p = 120\,y$ on the right edge (`Pressure`, `PNU`).

The exact solution is pure bending (plane stress, $\sigma_x = a y$):

$$u_x = \frac{a\,x\,y}{E},\qquad u_y = -\frac{a}{2E}\left(x^2 + \nu y^2\right) + \frac{a\,\nu}{2E}.$$

For $p\ge2$ the unsymmetric spectral element reproduces it to machine precision on all
three meshes (undistorted, angular, mid-edge). `tests/integration/LeeBathePatchTest.m` and
`tests/integration/AbaqusInputTest.m` check this.

The input files were reconstructed from the result workspaces that came with the code
(`tests/fixtures/make_abaqus_inputs.m`), and they reproduce the stored mesh, loads and
boundary conditions exactly.
