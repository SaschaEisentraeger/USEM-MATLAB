#!/usr/bin/env python3
"""Square panel of Rajendran (2010), Sections 5.5-5.7.

Geometry and mesh are generated with gmsh (Python API):
  * panel [0, 2] x [-1, 1] (built-in kernel),
  * structured (transfinite) mesh of 2N x 2N quadrilaterals, i.e. N x N
    super-elements of 2 x 2 elements (Rajendran 2010, Figs. 9 and 10),
  * 8-node elements (Mesh.SecondOrderIncomplete = 1); the 9th node of the
    9-node elements is added in the analysis (centroid of the corners).

Mesh distortion (Rajendran 2010, Fig. 9): the mid-side nodes of the four
internal edges of every super-element (nodes 17-20 of Mesh a) are rotated
counter-clockwise about the centre node of the super-element (node 21) by
theta. Meshes a-f: theta = 0, 90, 135, 175, 225, 270 deg; convergence
sequences A, B, C: theta = 0, 90, 135 deg with N = 1, 2, 4.

For every (N, theta) an ABAQUS input file is written (corner nodes first,
as required by pSpectralElement_NodesTopology) with the material, the node
sets Set-ROLLER ((0,-1), (0,1): u_x = 0) and Set-PIN ((0,0): u_x = u_y = 0)
and the corresponding boundary conditions. The loads (equilibrated nodal
forces on x = 0 and x = 2) are applied in the analysis.

Usage: python3 gen_panel_meshes.py [outdir]
"""
import math
import os
import sys

import gmsh
import numpy as np

LX, Y0, Y1 = 2.0, -1.0, 1.0
E_MOD, NU = 2000.0, 0.25     # E = 2000 reproduces the exact values of the paper
CASES = [(1, t) for t in (0, 90, 135, 175, 225, 270)] + \
        [(N, t) for N in (2, 4) for t in (0, 90, 135)]
EDGES = [(0, 1), (1, 2), (2, 3), (3, 0)]


def gmsh_mesh(N):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("panel")
    g = gmsh.model.geo
    p = [g.addPoint(0, Y0, 0), g.addPoint(LX, Y0, 0), g.addPoint(LX, Y1, 0), g.addPoint(0, Y1, 0)]
    lines = [g.addLine(p[i], p[(i + 1) % 4]) for i in range(4)]
    s = g.addPlaneSurface([g.addCurveLoop(lines)])
    for cv in lines:
        g.mesh.setTransfiniteCurve(cv, 2 * N + 1)
    g.mesh.setTransfiniteSurface(s, "Left", p)
    g.mesh.setRecombine(2, s)
    g.synchronize()
    gmsh.model.addPhysicalGroup(2, [s], name="PANEL")
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 1)
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(2)
    tags, coord, _ = gmsh.model.mesh.getNodes()
    xyz = dict(zip(tags, np.reshape(coord, (-1, 3))[:, :2]))
    etypes, _, enodes = gmsh.model.mesh.getElements(2, s)
    conn = np.reshape(enodes[0], (-1, 8))
    gmsh.finalize()
    corner = list(dict.fromkeys(conn[:, :4].ravel()))
    mid = list(dict.fromkeys(conn[:, 4:].ravel()))
    new = {t: i for i, t in enumerate(corner + mid)}
    nodes = np.array([xyz[t] for t in corner + mid])
    nodes = np.round(nodes * 2 * N) / (2 * N)      # exact grid coordinates
    nodes = np.round(nodes, 15)
    topo = np.vectorize(new.get)(conn)
    for e in range(topo.shape[0]):
        P = nodes[topo[e, :4]]
        area = 0.5 * np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1])
        assert area > 0
    return nodes, topo


def distort(nodes, topo, N, theta):
    nodes = nodes.copy()
    if theta == 0:
        return nodes
    a = math.radians(theta)
    R = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    h = 1.0 / N
    nc = int(topo[:, :4].max()) + 1
    ix = np.rint(nodes[:nc, 0] / h).astype(int)
    iy = np.rint((nodes[:nc, 1] - Y0) / h).astype(int)
    index = {(i, j): n for n, (i, j) in enumerate(zip(ix, iy))}
    midnode = {}
    for e in range(topo.shape[0]):
        for k, (i0, i1) in enumerate(EDGES):
            midnode[frozenset((topo[e, i0], topo[e, i1]))] = topo[e, 4 + k]
    for I in range(1, 2 * N, 2):
        for J in range(1, 2 * N, 2):
            c = index[(I, J)]
            C = nodes[c].copy()
            for di, dj in ((1, 0), (0, 1), (-1, 0), (0, -1)):
                m = midnode[frozenset((c, index[(I + di, J + dj)]))]
                nodes[m] = C + R @ (nodes[m] - C)
    return nodes


def write_inp(fname, nodes, topo, N, theta):
    x, y = nodes[:, 0], nodes[:, 1]
    roller = [i for i in range(len(x)) if abs(x[i]) < 1e-12 and abs(abs(y[i]) - 1) < 1e-12]
    pin = [i for i in range(len(x)) if abs(x[i]) < 1e-12 and abs(y[i]) < 1e-12]
    assert len(roller) == 2 and len(pin) == 1
    name = os.path.splitext(os.path.basename(fname))[0]
    with open(fname, "w") as f:
        f.write("*Heading\n")
        f.write(f"** Job name: {name} Model name: Model-1\n")
        f.write("** Square panel [0,2] x [-1,1] of Rajendran (2010), Secs. 5.5-5.7 (nodal rotation\n")
        f.write(f"** patch test / convergence test), {2*N} x {2*N} CPS8 elements, theta = {theta} deg.\n")
        f.write("** E = 2000 reproduces the exact solutions given in the paper (which states E = 1500).\n")
        f.write("** Loads (equilibrated nodal forces on x = 0 and x = 2) are applied in the analysis.\n")
        f.write("** Geometry and mesh: gmsh (studies/rajendran_panel/mesh/gen_panel_meshes.py)\n")
        f.write("*Preprint, echo=NO, model=NO, history=NO, contact=NO\n")
        f.write("**\n** PARTS\n**\n*Part, name=Panel\n*Node\n")
        for i, (a, b) in enumerate(nodes):
            f.write(f"{i+1:7d}, {a:.17g}, {b:.17g}\n")
        f.write("*Element, type=CPS8\n")
        for e, row in enumerate(topo):
            f.write(f"{e+1:7d}, " + ", ".join(str(int(n) + 1) for n in row) + "\n")
        f.write(f"*Nset, nset=Set-all, generate\n 1, {nodes.shape[0]}, 1\n")
        f.write(f"*Elset, elset=Set-all, generate\n 1, {topo.shape[0]}, 1\n")
        f.write("** Section: Section-1\n*Solid Section, elset=Set-all, material=Material-1\n,\n")
        f.write("*End Part\n**\n**\n** ASSEMBLY\n**\n*Assembly, name=Assembly\n**\n")
        f.write("*Instance, name=Panel-1, part=Panel\n*End Instance\n**\n")
        f.write("*Nset, nset=Set-ROLLER, instance=Panel-1\n " + ", ".join(str(i + 1) for i in roller) + "\n")
        f.write("*Nset, nset=Set-PIN, instance=Panel-1\n " + ", ".join(str(i + 1) for i in pin) + "\n")
        f.write("*End Assembly\n**\n** MATERIALS\n**\n*Material, name=Material-1\n")
        f.write(f"*Density\n 1.,\n*Elastic\n {E_MOD:g}, {NU:g}\n")
        f.write("** ----------------------------------------------------------------\n")
        f.write("**\n** STEP: Step-1\n**\n*Step, name=Step-1, nlgeom=NO\n*Static\n1., 1., 1e-05, 1.\n")
        f.write("**\n** BOUNDARY CONDITIONS\n**\n")
        f.write("** Name: BC-1 Type: Displacement/Rotation\n*Boundary\nSet-ROLLER, 1, 1\n")
        f.write("** Name: BC-2 Type: Displacement/Rotation\n*Boundary\nSet-PIN, 1, 2\n")
        f.write("**\n** OUTPUT REQUESTS\n**\n*Restart, write, frequency=0\n")
        f.write("*Output, field, variable=PRESELECT\n*Output, history, variable=PRESELECT\n*End Step\n")


def edge_polylines(nodes, topo, nseg=24):
    s = np.linspace(-1.0, 1.0, nseg + 1)
    L = np.c_[0.5 * s * (s - 1), 1 - s**2, 0.5 * s * (s + 1)]
    seen, lines = set(), []
    for e in range(topo.shape[0]):
        for k, (i0, i1) in enumerate(EDGES):
            key = frozenset((topo[e, i0], topo[e, i1]))
            if key not in seen:
                seen.add(key)
                lines.append(L @ nodes[[topo[e, i0], topo[e, 4 + k], topo[e, i1]]])
    return lines


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "inp")
    plotdir = os.path.join(outdir, "..", "plotdata")
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(plotdir, exist_ok=True)
    base = {}
    for N, theta in CASES:
        if N not in base:
            base[N] = gmsh_mesh(N)
        nodes0, topo = base[N]
        nodes = distort(nodes0, topo, N, theta)
        fname = os.path.join(outdir, f"Panel_N{N}_t{theta:03d}.inp")
        write_inp(fname, nodes, topo, N, theta)
        with open(os.path.join(plotdir, f"mesh_N{N}_t{theta:03d}.dat"), "w") as f:
            f.write("x y\n")
            for P in edge_polylines(nodes, topo):
                for a, b in P:
                    f.write(f"{a:.6f} {b:.6f}\n")
                f.write("\n")
        print(f"N = {N}, theta = {theta:3d}: {nodes.shape[0]} nodes, {topo.shape[0]} elements")


if __name__ == "__main__":
    main()
