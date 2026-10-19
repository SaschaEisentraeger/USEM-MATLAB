#!/usr/bin/env python3
"""Thick cylinder under internal pressure (Rajendran 2010, Section 5.9).

Geometry and mesh are generated with gmsh (Python API):
  * quarter annulus a <= r <= b, 0 <= theta <= pi/2 (built-in kernel,
    exact circular arcs),
  * structured (transfinite) mesh of 2N x 2N quadrilaterals, i.e. N x N
    "super-elements" of 2 x 2 elements each (Rajendran 2010, Fig. 19),
  * 8-node serendipity elements (Mesh.SecondOrderIncomplete = 1).

Mesh distortion (mesh sequence B of Rajendran 2010): the mid-side nodes of
the four internal edges of every super-element are rotated counter-
clockwise about the centre node of the super-element by the angle alpha.
The boundary of the domain is not changed.

For every (N, alpha) an ABAQUS input file is written that can be used with
code/Start_Sim.m directly (corner nodes numbered first, as required by
pSpectralElement_NodesTopology). The file contains
  * node and element sets Set-XSYM (x = 0), Set-YSYM (y = 0),
    Set-OUTER (r = b) and Pressure (r = a),
  * symmetry boundary conditions (Set-XSYM: u_x = 0, Set-YSYM: u_y = 0),
  * the internal pressure p on the inner arc (*Dsload, P).
In addition the mesh edges are written as polylines (pgfplots tables) for
the report.

Usage: python3 gen_cylinder_meshes.py [outdir]
"""
import math
import os
import sys

import gmsh
import numpy as np

A_IN, B_OUT = 1.0, 3.0           # inner and outer radius
E_MOD, NU = 1000.0, 0.25         # material (plane stress)
PRESSURE = 1000.0                # internal pressure
N_LIST = [1, 2, 4, 8]            # super-elements per direction
ALPHA_LIST = [0, 15, 30, 45, 60, 90]   # distortion angles (degree)
# corner pairs of the edges of an 8-node element (ABAQUS CPS8 ordering)
EDGES = [(0, 1), (1, 2), (2, 3), (3, 0)]


def gmsh_mesh(N):
    """Undistorted 8-node mesh of the quarter annulus with 2N x 2N elements.

    Returns corner-first renumbered nodes (n x 2) and connectivity (m x 8).
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("thick_cylinder")
    g = gmsh.model.geo
    c = g.addPoint(0, 0, 0)
    p1 = g.addPoint(A_IN, 0, 0)
    p2 = g.addPoint(B_OUT, 0, 0)
    p3 = g.addPoint(0, B_OUT, 0)
    p4 = g.addPoint(0, A_IN, 0)
    l1 = g.addLine(p1, p2)              # y = 0
    c2 = g.addCircleArc(p2, c, p3)      # r = b
    l3 = g.addLine(p3, p4)              # x = 0
    c4 = g.addCircleArc(p4, c, p1)      # r = a
    loop = g.addCurveLoop([l1, c2, l3, c4])
    s = g.addPlaneSurface([loop])
    for cv in (l1, c2, l3, c4):
        g.mesh.setTransfiniteCurve(cv, 2 * N + 1)
    g.mesh.setTransfiniteSurface(s, "Left", [p1, p2, p3, p4])
    g.mesh.setRecombine(2, s)
    g.synchronize()
    gmsh.model.addPhysicalGroup(2, [s], name="CYLINDER")
    gmsh.model.addPhysicalGroup(1, [l1], name="YSYM")
    gmsh.model.addPhysicalGroup(1, [c2], name="OUTER")
    gmsh.model.addPhysicalGroup(1, [l3], name="XSYM")
    gmsh.model.addPhysicalGroup(1, [c4], name="INNER")
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(2)

    tags, coord, _ = gmsh.model.mesh.getNodes()
    xyz = dict(zip(tags, np.reshape(coord, (-1, 3))[:, :2]))
    etypes, _, enodes = gmsh.model.mesh.getElements(2, s)
    assert len(etypes) == 1
    name, _, order, nnode, _, _ = gmsh.model.mesh.getElementProperties(etypes[0])
    assert nnode == 8, name
    conn = np.reshape(enodes[0], (-1, 8))   # gmsh QUA_8 = ABAQUS CPS8 order
    gmsh.finalize()

    # Corner nodes first (pSpectralElement_NodesTopology), then mid-side nodes
    corner = list(dict.fromkeys(conn[:, :4].ravel()))
    mid = list(dict.fromkeys(conn[:, 4:].ravel()))
    new = {t: i for i, t in enumerate(corner + mid)}
    nodes = np.array([xyz[t] for t in corner + mid])
    topo = np.vectorize(new.get)(conn)

    # Snap boundary nodes onto the exact boundary (gmsh positions carry
    # round-off of O(1e-9))
    r = np.hypot(nodes[:, 0], nodes[:, 1])
    th = np.arctan2(nodes[:, 1], nodes[:, 0])
    for R in (A_IN, B_OUT):
        on = np.abs(r - R) < 1e-6
        nodes[on] = R * np.c_[np.cos(th[on]), np.sin(th[on])]
    nodes[np.abs(nodes[:, 0]) < 1e-6, 0] = 0.0
    nodes[np.abs(nodes[:, 1]) < 1e-6, 1] = 0.0
    # Element orientation must be counter-clockwise
    for e in range(topo.shape[0]):
        P = nodes[topo[e, :4]]
        area = 0.5 * np.sum(P[:, 0] * np.roll(P[:, 1], -1) -
                            np.roll(P[:, 0], -1) * P[:, 1])
        assert area > 0
    return nodes, topo


def distort(nodes, topo, N, alpha_deg):
    """Rotate the internal mid-side nodes of each super-element about its
    centre node by alpha (counter-clockwise)."""
    nodes = nodes.copy()
    if alpha_deg == 0:
        return nodes
    a = math.radians(alpha_deg)
    R = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    # Structured indices of the corner nodes
    dr = (B_OUT - A_IN) / (2 * N)
    dth = 0.5 * math.pi / (2 * N)
    ncorner = int(topo[:, :4].max()) + 1
    P = nodes[:ncorner]
    ir = np.rint((np.hypot(P[:, 0], P[:, 1]) - A_IN) / dr).astype(int)
    it = np.rint(np.arctan2(P[:, 1], P[:, 0]) / dth).astype(int)
    index = {(i, j): n for n, (i, j) in enumerate(zip(ir, it))}
    # Mid-side node of each edge (pair of corner nodes)
    midnode = {}
    for e in range(topo.shape[0]):
        for k, (i0, i1) in enumerate(EDGES):
            key = frozenset((topo[e, i0], topo[e, i1]))
            midnode[key] = topo[e, 4 + k]
    moved = 0
    for I in range(1, 2 * N, 2):
        for J in range(1, 2 * N, 2):
            cnode = index[(I, J)]
            C = nodes[cnode]
            for (di, dj) in ((1, 0), (0, 1), (-1, 0), (0, -1)):
                m = midnode[frozenset((cnode, index[(I + di, J + dj)]))]
                nodes[m] = C + R @ (nodes[m] - C)
                moved += 1
    assert moved == 4 * N * N
    return nodes


def boundary_sets(nodes, topo):
    x, y = nodes[:, 0], nodes[:, 1]
    r = np.hypot(x, y)
    tol = 1e-10
    sets = {
        "XSYM": np.flatnonzero(np.abs(x) < tol),
        "YSYM": np.flatnonzero(np.abs(y) < tol),
        "INNER": np.flatnonzero(np.abs(r - A_IN) < tol),
        "OUTER": np.flatnonzero(np.abs(r - B_OUT) < tol),
    }
    # Elements with an edge (three CPS8 nodes) on each boundary part
    elsets = {}
    for key, ids in sets.items():
        on = set(ids.tolist())
        elsets[key] = np.array([e for e in range(topo.shape[0])
                                if sum(int(n) in on for n in topo[e]) >= 3])
    return sets, elsets


def write_list(f, ids, per_line=16):
    ids = [int(i) + 1 for i in ids]
    for k in range(0, len(ids), per_line):
        f.write(" " + ", ".join(str(i) for i in ids[k:k + per_line]) + "\n")


def write_inp(fname, nodes, topo, N, alpha):
    sets, elsets = boundary_sets(nodes, topo)
    name = os.path.splitext(os.path.basename(fname))[0]
    with open(fname, "w") as f:
        f.write("*Heading\n")
        f.write(f"** Job name: {name} Model name: Model-1\n")
        f.write("** Thick cylinder under internal pressure (Rajendran 2010, Sec. 5.9),\n")
        f.write(f"** quarter model a = {A_IN:g}, b = {B_OUT:g}, p = {PRESSURE:g}, "
                f"E = {E_MOD:g}, nu = {NU:g}, plane stress.\n")
        f.write(f"** {2*N} x {2*N} CPS8 elements ({N} x {N} super-elements), "
                f"distortion angle alpha = {alpha} deg.\n")
        f.write("** Geometry and mesh: gmsh (studies/thick_cylinder/mesh/gen_cylinder_meshes.py)\n")
        f.write("*Preprint, echo=NO, model=NO, history=NO, contact=NO\n")
        f.write("**\n** PARTS\n**\n*Part, name=Cylinder\n*Node\n")
        for i, (x, y) in enumerate(nodes):
            f.write(f"{i+1:7d}, {x:.17g}, {y:.17g}\n")
        f.write("*Element, type=CPS8\n")
        for e, row in enumerate(topo):
            f.write(f"{e+1:7d}, " + ", ".join(str(int(n) + 1) for n in row) + "\n")
        f.write(f"*Nset, nset=Set-all, generate\n 1, {nodes.shape[0]}, 1\n")
        f.write(f"*Elset, elset=Set-all, generate\n 1, {topo.shape[0]}, 1\n")
        f.write("** Section: Section-1\n*Solid Section, elset=Set-all, material=Material-1\n,\n")
        f.write("*End Part\n**\n**\n** ASSEMBLY\n**\n*Assembly, name=Assembly\n**\n")
        f.write("*Instance, name=Cylinder-1, part=Cylinder\n*End Instance\n**\n")
        # Read_boundary_condition_ABAQUS expects set names containing "Set"
        # (and an element set of the same name, from which
        # pSpectralElement_BCs adds the high-order edge nodes)
        for key, name in (("XSYM", "Set-XSYM"), ("YSYM", "Set-YSYM"),
                          ("OUTER", "Set-OUTER"), ("INNER", "Pressure")):
            f.write(f"*Nset, nset={name}, instance=Cylinder-1\n")
            write_list(f, sets[key])
            f.write(f"*Elset, elset={name}, instance=Cylinder-1\n")
            write_list(f, elsets[key])
        f.write("*End Assembly\n**\n** MATERIALS\n**\n*Material, name=Material-1\n")
        f.write("*Density\n 1.,\n*Elastic\n")
        f.write(f" {E_MOD:g}, {NU:g}\n")
        f.write("** ----------------------------------------------------------------\n")
        f.write("**\n** STEP: Step-1\n**\n*Step, name=Step-1, nlgeom=NO\n*Static\n1., 1., 1e-05, 1.\n")
        f.write("**\n** BOUNDARY CONDITIONS\n**\n")
        f.write("** Name: BC-1 Type: Displacement/Rotation\n*Boundary\nSet-XSYM, 1, 1\n")
        f.write("** Name: BC-2 Type: Displacement/Rotation\n*Boundary\nSet-YSYM, 2, 2\n")
        f.write("**\n** LOADS\n**\n** Name: Load-1   Type: Pressure\n*Dsload\n")
        f.write(f"Pressure, P, {PRESSURE:g}\n")
        f.write("**\n** OUTPUT REQUESTS\n**\n*Restart, write, frequency=0\n")
        f.write("*Output, field, variable=PRESELECT\n*Output, history, variable=PRESELECT\n")
        f.write("*End Step\n")


def edge_polylines(nodes, topo, nseg=24):
    """Element edges as polylines (parabolas through the three edge nodes;
    exact arcs on r = a and r = b, as used in the analysis)."""
    s = np.linspace(-1.0, 1.0, nseg + 1)
    L = np.c_[0.5 * s * (s - 1), 1 - s**2, 0.5 * s * (s + 1)]
    seen, lines = set(), []
    for e in range(topo.shape[0]):
        for k, (i0, i1) in enumerate(EDGES):
            n0, n1, nm = topo[e, i0], topo[e, i1], topo[e, 4 + k]
            key = frozenset((n0, n1))
            if key in seen:
                continue
            seen.add(key)
            P = L @ nodes[[n0, nm, n1]]
            r = np.hypot(*nodes[[n0, nm, n1]].T)
            for R in (A_IN, B_OUT):
                if np.all(np.abs(r - R) < 1e-10):
                    t0, t1 = (math.atan2(*nodes[n][::-1]) for n in (n0, n1))
                    t = 0.5 * (t0 + t1) + 0.5 * s * (t1 - t0)
                    P = R * np.c_[np.cos(t), np.sin(t)]
            lines.append(P)
    return lines


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "inp")
    os.makedirs(outdir, exist_ok=True)
    plotdir = os.path.join(outdir, "..", "plotdata")
    os.makedirs(plotdir, exist_ok=True)
    for N in N_LIST:
        nodes0, topo = gmsh_mesh(N)
        for alpha in ALPHA_LIST:
            nodes = distort(nodes0, topo, N, alpha)
            fname = os.path.join(outdir, f"ThickCylinder_N{N}_a{alpha:03d}.inp")
            write_inp(fname, nodes, topo, N, alpha)
            if N <= 2:
                with open(os.path.join(plotdir, f"mesh_N{N}_a{alpha:03d}.dat"), "w") as f:
                    f.write("x y\n")
                    for P in edge_polylines(nodes, topo):
                        for x, y in P:
                            f.write(f"{x:.6f} {y:.6f}\n")
                        f.write("\n")
            print(f"N = {N}, alpha = {alpha:3d}: {nodes.shape[0]} nodes, "
                  f"{topo.shape[0]} elements -> {os.path.basename(fname)}")


if __name__ == "__main__":
    main()
