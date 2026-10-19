# Literature

This folder collects the literature on the unsymmetric finite element method (UFEM) that the
unsymmetric spectral element code (USEM) is built on. **The PDFs are not under version
control.** They are publisher-copyrighted or unpublished (see [`.gitignore`](../.gitignore)).
Use the DOIs below to get them, or ask the maintainer. All entries are available as BibTeX
in [`references.bib`](references.bib).

The references that the implementation follows most closely are:

| Code | Reference |
|---|---|
| Unsymmetric (Petrov–Galerkin) element, adjugate form of $\mathbf B_L$ (`Ke_2D`, `Me_2D`, `KMe_2D`) | Rajendran (2010); Rajendran & Liew (2003) |
| Metric trial functions in corotational skew coordinates (`ShapeFunctions2D_Metric`, `CorotationalSkewCoordinates`) | Xie, Sze & Zhou (2016); Yuan, Huang & Pian (1993); Ooi, Rajendran & Yeo (2008) |
| Lee & Bathe patch tests (`code/InputFilesAbaqus`, `SpatialLoadingFunction`) | Lee & Bathe (1993); Rajendran & Liew (2003) |
| Stability of time integration with unsymmetric matrices | Wang & Hillman (2022); Zähringer & Betsch (2026) |

See [`docs/01-theory.md`](../docs/01-theory.md) for the formulation as implemented.

Element acronyms used below:
PP = parametric test / parametric trial (classical isoparametric),
PM = parametric test / metric trial (unsymmetric, "US-"),
MM = metric / metric,
ATF = analytical trial functions.

---

## 1. Axisymmetric UFEM

| Reference | Key content | DOI / source |
|---|---|---|
| Zähringer F, Betsch P (2024). *Finite element formulations for axisymmetric solids: a new approach based on the Petrov–Galerkin method.* ECCOMAS Congress 2024, Lisbon. | 8-node **PG-Q8**: Serendipity test functions, metric trial functions in **skew coordinates**. It is exact for quadratic fields (circular plate in pure bending, cylindrical shell under body force) on distorted meshes. For the **thick-walled cylinder** (non-polynomial $u_r$) PG-Q8 is *not* generally better than Q8. A 4-node PG-EAS variant appears infeasible in axisymmetry. | proceedings, no DOI |
| Dou LZ, Huang YQ, Yang YF, Ju YL, Chen HB (2025). *A novel Petrov–Galerkin 4-node quadrilateral element with radial polynomial interpolation for linear elastic analysis.* IJNME 126(10):e70049. | **UQ4-RPIM**. Sec. 2.4 extends it to axisymmetry (r, z) with a modified polynomial basis for hollow bodies. Rotating disc and thick cylinder (ν = 0.3, 0.49). | [10.1002/nme.70049](https://doi.org/10.1002/nme.70049) |
| Huang Y, Huan Y, Chen H (2020). *An incompatible and unsymmetric four-node quadrilateral plane element with high numerical performance.* IJNME 121(15):3382–3396. | **IUQ4**. Sec. 3.7 extends it to axisymmetric problems, with incompatible test modes modified to pass the axisymmetric patch test. | [10.1002/nme.6363](https://doi.org/10.1002/nme.6363) |

## 2. Foundations of UFEM / Petrov–Galerkin

| Reference | Key content | DOI |
|---|---|---|
| Rajendran S, Liew KM (2003). IJNME 58(11):1713–1748. | **US-QUAD8**, the seminal paper. Serendipity (compatibility-enforcing) test functions combined with metric (completeness-enforcing) trial functions. Exact quadratic fields under angular, curved-edge and mid-node distortion. | [10.1002/nme.836](https://doi.org/10.1002/nme.836) |
| Rajendran S, Liew KM (2000). SEM 10(2):93–110. | Higher-order completeness conditions. Explains why distorted Q8 elements lose quadratic completeness. | [10.12989/sem.2000.10.2.093](https://doi.org/10.12989/sem.2000.10.2.093) |
| Rajendran S (2010). CMAME 199:1044–1063. | Three sufficient conditions for distortion immunity (compatible test, complete trial, no $\det\mathbf J$ in the denominator). PP/PM/MM naming. | [10.1016/j.cma.2009.11.017](https://doi.org/10.1016/j.cma.2009.11.017) |
| Ooi ET (2006). PhD thesis, NTU Singapore. | Comprehensive study of US-QUAD8/US-HEXA20: strong-form patch tests, GNL, mixed US-QUAD9/3, iterative corrections. Discusses frame dependence, P-matrix singularity and possible complex eigenvalues. | [10.32657/10356/5254](https://doi.org/10.32657/10356/5254) |
| Prathap G, Manju S, Senthilkumar V (2007). SEM 26(1):31–42. | PM formulation is not variationally correct, but recovers stresses as a best fit in metric space. | [10.12989/sem.2007.26.1.031](https://doi.org/10.12989/sem.2007.26.1.031) |
| Cen S, Wu CJ, Li Z, Shang Y, Li C (2019). Eng Comput 36(8):2811–2834. | Review of hybrid stress/displacement-function elements and the "improved UFEM" (US-ATF family). | [10.1108/EC-10-2018-0479](https://doi.org/10.1108/EC-10-2018-0479) |
| Chen H (2021). Master thesis, OVGU Magdeburg. | Abaqus UEL (Fortran77) implementation of PM-QUAD8 and PM-HEXA20. Cantilever and Kirsch-plate benchmarks under distortion. | thesis |

## 3. Skew coordinates and frame invariance (basis of the trial functions in this code)

| Reference | Key content | DOI |
|---|---|---|
| Xie Q, Sze KY, Zhou YX (2016). IJMMD 12(1):53–70. | **UQ8\***, a frame-invariant UQ8 whose metric trial functions are written in non-dimensional **skew coordinates**. Also introduces Trefftz-UFEM TQ4/TH8. | [10.1007/s10999-014-9289-3](https://doi.org/10.1007/s10999-014-9289-3) |
| Yuan KY, Huang YS, Pian THH (1993). IJNME 36(10):1747–1763. | Origin of the skew (oblique) coordinates defined at the element centre. | [10.1002/nme.1620361009](https://doi.org/10.1002/nme.1620361009) |
| Ooi ET, Rajendran S, Yeo JH (2008). CNME 24(11):1203–1217. | US-QUAD8 is frame-dependent for the cubic terms $x^2y$, $xy^2$, and its P-matrix can become singular. Remedies: an edge-aligned local frame and a small random rotation. | [10.1002/cnm.1026](https://doi.org/10.1002/cnm.1026) |
| Cowan T, Coombs WM (2014). CMAME 275:189–203. | Metric trial functions in a local convected frame, giving rotational invariance. Finite-deformation elasto-plasticity. | [10.1016/j.cma.2014.02.016](https://doi.org/10.1016/j.cma.2014.02.016) |
| Wisniewski K, Turska E (2008). IJNME 76(6):798–836. | Hellinger–Reissner 4-node elements with stresses in skew coordinates. | [10.1002/nme.2343](https://doi.org/10.1002/nme.2343) |

## 4. 2D plane elements

| Reference | Element | DOI |
|---|---|---|
| Rajendran S, Subramanian S (2004). SEM 17(6):767–788. Comparison of PP/MM/PM/MP 8-node elements. | PP, MM, PM, MP | [10.12989/sem.2004.17.6.767](https://doi.org/10.12989/sem.2004.17.6.767) |
| Liew KM, Rajendran S, Wang J (2006). CMAME 195:1207–1223. 6-node triangle, incl. free vibration with an unsymmetric mass matrix. | US-TRIA6 | [10.1016/j.cma.2005.04.012](https://doi.org/10.1016/j.cma.2005.04.012) |
| Cen S, Zhou GH, Fu XR (2012). IJNME 91(2):158–185. Airy-function analytical trial functions. | US-ATFQ8 | [10.1002/nme.4260](https://doi.org/10.1002/nme.4260) |
| Cen S, Zhou PL, Li CF, Wu CJ (2015). IJNME 103(7):469–500. "Breaking through" MacNeal's theorem. | US-ATFQ4 | [10.1002/nme.4899](https://doi.org/10.1002/nme.4899) |
| Shang Y, Ouyang W (2018). IJNME 113(10):1589–1606. Drilling DOFs, stress trial functions. | US-Q4 | [10.1002/nme.5711](https://doi.org/10.1002/nme.5711) |
| Shang Y, Cen S, Zhou MJ (2018). Acta Mech 229:5031–5049. Orthotropic, Airy stress solutions. | US-Q8 | [10.1007/s00707-018-2291-3](https://doi.org/10.1007/s00707-018-2291-3) |
| Huang YQ, Yang YF, Wang JZ, Liu XC, Chen HB (2022). IJNME 123(1):101–127. Wilson incompatible modes, explicit trial functions without inverting P. | UQ6, UH11 | [10.1002/nme.6849](https://doi.org/10.1002/nme.6849) |
| Yang YF, Huang YQ, Wang JZ, Liu XC, Chen HB (2023). IJNME 124(12):2731–2755. Cubic completeness through incompatible modes / EAS. | UQ10, UQ10/E4 | [10.1002/nme.7226](https://doi.org/10.1002/nme.7226) |

## 5. 3D solid and solid-shell elements

| Reference | Element | DOI |
|---|---|---|
| Ooi ET, Rajendran S, Yeo JH (2004). IJNME 60:2501–2530. | US-HEXA20 | [10.1002/nme.1056](https://doi.org/10.1002/nme.1056) |
| Zhou PL, Cen S, Huang JB, Li CF, Zhang Q (2017). IJNME 109(8):1130–1158. | US-ATFH8 | [10.1002/nme.5318](https://doi.org/10.1002/nme.5318) |
| Huang J, Cen S, Li Z, Li CF (2018). IJNME 116(12-13):759–783. | US-ATFHS8 | [10.1002/nme.5945](https://doi.org/10.1002/nme.5945) |
| Wu NC, Huang YQ, Chen HB (2023). CMAME 417:116396. | UHS11 | [10.1016/j.cma.2023.116396](https://doi.org/10.1016/j.cma.2023.116396) |
| He PQ, Sun Q, Liang K (2019). Comput Mech 63:755–781. Generalized modal element method. | US-MEM8S, S-MEM8S | [10.1007/s00466-018-1618-1](https://doi.org/10.1007/s00466-018-1618-1) |
| Pfefferkorn R, Betsch P (2021). IJNME 122(23):6924–6954. Petrov–Galerkin EAS, revisits MacNeal's theorem. | Q1U/E4, H1U/E12 | [10.1002/nme.6817](https://doi.org/10.1002/nme.6817) |

## 6. Nonlinear extensions (geometric, material, hyperelastic)

| Reference | Element | DOI |
|---|---|---|
| Ooi ET, Rajendran S, Yeo JH (2007). Eng Comput 24(4):407–431. | US-QUAD8/US-HEXA20 (TL) | [10.1108/02644400710748715](https://doi.org/10.1108/02644400710748715) |
| Li Z, Cen S, Wu CJ, Shang Y, Li CF (2018). IJNME 114(9):931–954. | US-ATFQ4 (UL) | [10.1002/nme.5771](https://doi.org/10.1002/nme.5771) |
| Li Z, Huang J, Cen S, Li CF (2019). IJNME 120(5):580–606. | US-ATFHS8 (CR) | [10.1002/nme.6149](https://doi.org/10.1002/nme.6149) |
| Li Z, Cen S, Huang J, Li CF (2020). IJNME 121(16):3702–3721. | US-ATFQ4 hyperelastic | [10.1002/nme.6378](https://doi.org/10.1002/nme.6378) |
| Chen X, Chen W, Xu J, Tu W, Rajendran S, Lu J (2021). IJNME 122(20):5540–5557. | US-QUAD8 elasto-plastic | [10.1002/nme.6763](https://doi.org/10.1002/nme.6763) |
| Pfefferkorn R, Betsch P (2023). IJNME 124(6):1307–1343. | PG-EAS H1U/h12 | [10.1002/nme.7166](https://doi.org/10.1002/nme.7166) |
| Wang YF, Cen S, Li CF, Zhang Q (2023). Int J Mech Sci 260:108647. | US-ATFH8 soft tissue | [10.1016/j.ijmecsci.2023.108647](https://doi.org/10.1016/j.ijmecsci.2023.108647) |
| Wu NC, Huang YQ, Chen HB (2026). FEAD 260:104610. | UHS11 hyperelastic | [10.1016/j.finel.2026.104610](https://doi.org/10.1016/j.finel.2026.104610) |

## 7. Dynamics and time integration

These papers matter for the dynamic, modal and frequency-domain analyses (`call_Dynamic.m`, `call_Modal.m`, `call_FrequencyDomain.m`).

| Reference | Key content | DOI |
|---|---|---|
| Wang J, Hillman MC (2022). CMAME 393:114738. | With non-symmetric $\mathbf M$/$\mathbf K$, energy is not conserved, and Newmark can be **unconditionally unstable** (negative generalized eigenvalues). Proposes an energy-conserving scheme. | [10.1016/j.cma.2022.114738](https://doi.org/10.1016/j.cma.2022.114738) |
| Zähringer F, Betsch P (2026). Preprint, KIT. | UQ8\* in elastodynamics. Standard Newmark (A1) is unstable, Wang–Hillman (A2) can drift, and a new mixed velocity–strain formulation (A3) conserves energy and stays distortion-insensitive. | preprint |
| Kumar S (2011). FEAD 47(4):373–377. | First transient (impact) application of PM20 with Newmark; stability is not examined. | [10.1016/j.finel.2010.12.016](https://doi.org/10.1016/j.finel.2010.12.016) |

## 8. Plates, beams and shells

| Reference | Element | DOI |
|---|---|---|
| Ju YL, Huang YQ, Chen HB (2025). FEAD 251:104408. | UQP4-RPIM | [10.1016/j.finel.2025.104408](https://doi.org/10.1016/j.finel.2025.104408) |
| Shang Y, Liu YD, Liu SX (2022). Eng Comput 38(Suppl 2):S1065–S1079. | US-TP4 | [10.1007/s00366-020-01254-x](https://doi.org/10.1007/s00366-020-01254-x) |

## 9. Couple-stress and flexoelectric continua

| Reference | Element | DOI |
|---|---|---|
| Shang Y, Qian ZH, Cen S, Li CF (2019). IJNME 119(9):807–825. | US-Q4-CS | [10.1002/nme.6073](https://doi.org/10.1002/nme.6073) |
| Shang Y, Li CF, Jia KY (2020). IJNME 121(12):2683–2700. | 8-node hex (MCST) | [10.1002/nme.6325](https://doi.org/10.1002/nme.6325) |
| Shang Y, Wu H (2022). Eng Comput 39(3):1097–1117. | US-Q8, US-Q8L | [10.1108/EC-03-2021-0145](https://doi.org/10.1108/EC-03-2021-0145) |
| Wu HP, Shang Y, Cen S, Li CF (2023). EABE 147:302–319. | Penalty C0 Q8/H20 (CCST) | [10.1016/j.enganabound.2022.12.008](https://doi.org/10.1016/j.enganabound.2022.12.008) |
| Deng Z, Liu SX, Cen S, Sun M, Shang Y (2025). IJNME 126(11):e70045. | CR flexoelectric membrane | [10.1002/nme.70045](https://doi.org/10.1002/nme.70045) |
| Shang Y, Li XJ, Pan MD, Cen S, Chen BB, Li CF (2025). IJNME 126(22):e70179. | Trefftz Kirchhoff plate (CCST) | [10.1002/nme.70179](https://doi.org/10.1002/nme.70179) |

## 10. Mesh distortion: theory and assessment

| Reference | Key content | DOI |
|---|---|---|
| Lee NS, Bathe KJ (1993). IJNME 36(20):3553–3576. | Classic study: angular distortion destroys the quadratic completeness of Serendipity elements. | [10.1002/nme.1620362009](https://doi.org/10.1002/nme.1620362009) |
| MacNeal RH (1992). IJNME 35(8):1589–1601. | "MacNeal's theorem": a symmetric element must either fail the patch test or lock. | [10.1002/nme.1620350804](https://doi.org/10.1002/nme.1620350804) |
| Rajendran S, Ooi ET, Yeo JH (2007). CNME 23(2):157–168. | Strong-form patch tests: only PM-QUAD8 passes all of them. | [10.1002/cnm.893](https://doi.org/10.1002/cnm.893) |
| Prathap G, Senthilkumar V, Manju S (2006). Sadhana 31(5):505–514. | Best-fit paradigm explanation of the distortion immunity of PM elements. | [10.1007/BF02715909](https://doi.org/10.1007/BF02715909) |
| Mukherjee S, Manju S (2011). SEM 38(3):261–281. | Variationally correct, distortion-immune 3-node bar (MBAR3). | [10.12989/sem.2011.38.3.261](https://doi.org/10.12989/sem.2011.38.3.261) |
| Manju S, Mukherjee S (2019). SEM 69(6):615–626. | Function-space analysis of the distorted 3-node Timoshenko metric beam. | [10.12989/sem.2019.69.6.615](https://doi.org/10.12989/sem.2019.69.6.615) |
| Engvall L, Evans JA (2020). CMAME 371:113305. | Mesh quality metrics for high-order (Bernstein–Bézier) elements. | [10.1016/j.cma.2020.113305](https://doi.org/10.1016/j.cma.2020.113305) |
| Quiriny A, Kučera V, Lambrechts J, Moës N, Remacle JF (2024). Manuscript. | Tempered FEM: a lower bound on $\det\mathbf J$ instead of a change of test/trial spaces. | preprint |
| Quiriny A, Lambrechts J, Moës N, Kučera V, Remacle JF (2026). SSRN preprint. | TFEM for sliver tetrahedra. | preprint |

---

## Synthesis: what the literature says, and what it means for this code

1. **Why UFEM.** Conventional (Bubnov–Galerkin) Serendipity elements lose their quadratic
   completeness under angular and curved-edge distortion (Lee & Bathe 1993). MacNeal's
   theorem shows that no symmetric element can both pass the patch test and be
   distortion-insensitive in bending. Petrov–Galerkin elements avoid this. They combine
   compatible (parametric) test functions with Cartesian-complete (metric) trial functions
   (Rajendran & Liew 2003; Rajendran 2010).
2. **Sufficient conditions for distortion immunity** (Rajendran 2010):
   - compatible test functions;
   - trial functions that are complete in Cartesian coordinates;
   - no $\det\mathbf J$ in the denominator of the integrand (adjugate form).

   This code meets all three for arbitrary polynomial degree $p$. The trial functions are
   tensor-product Lagrange polynomials in skew coordinates, which contain all polynomials of
   total degree $\le p$ in $(x, y)$, and $\mathbf B_L$ is built with $\operatorname{adj}\mathbf J$.
   Any displacement field of total degree $\le p$ is therefore reproduced exactly on
   arbitrarily distorted meshes (verified by `tests/integration/PatchTest.m`).
3. **Frame invariance.** Metric trial functions in global coordinates are frame-dependent,
   and their interpolation (Vandermonde) matrix can become singular (Ooi et al. 2008).
   Writing them in **skew coordinates** (Yuan, Huang & Pian 1993; Xie, Sze & Zhou 2016)
   fixes the frame dependence. That is the default here (`ShapeFunctions2D_Metric`,
   `method = 1`). The Vandermonde matrix can still be singular or ill-conditioned for
   degenerate elements, e.g. three nearly collinear corners.
4. **Variational inconsistency.** The PM formulation is not variationally correct (Prathap et
   al. 2007): the stiffness matrix is unsymmetric and Bubnov–Galerkin energy optimality is
   lost. For fields outside the trial space, the unsymmetric element can be slightly less
   accurate than the classical one on coarse meshes (Zähringer & Betsch 2024), but it keeps
   the optimal convergence rate on distorted meshes.
5. **Dynamics.** Unsymmetric $\mathbf K$ (and $\mathbf M$) break the equivalence between
   energy and momentum conservation. Standard Newmark can become unstable, and generalized
   eigenvalues can be complex (Wang & Hillman 2022; Ooi 2006; Zähringer & Betsch 2026).
   Interpret dynamic, modal and frequency-domain results of the unsymmetric formulation with
   care. With GLL nodes and nodal quadrature, the unsymmetric mass matrix is nevertheless
   diagonal, because both the parametric and the metric functions are interpolatory at the
   quadrature points.
6. **Open questions.**
   - Element-level stability, i.e. the conditioning of the metric Vandermonde matrix and the
     discrete inf-sup constant of the unsymmetric system, on strongly distorted high-order
     elements.
   - Energy-consistent time integration for unsymmetric spectral elements.
   - A convergence theory for these variationally inconsistent schemes.
   - Extension to triangles (`ShapeFunctions2D_Metric_Tria` is an unfinished start), 3D and
     axisymmetry (Section 1).
