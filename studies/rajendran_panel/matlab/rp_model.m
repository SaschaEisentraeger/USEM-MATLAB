function model = rp_model(inpFile, pdeg)
%RP_MODEL Spectral element model of the panel from a gmsh/ABAQUS file.
%   Edge nodes on the parabolic 8-node edges, interior nodes from the
%   bilinear map of the four corner nodes. Built with
%   pSpectralElement_NodesTopology on the 4-node corner mesh, the edge nodes
%   are then moved onto the 8-node edges.
[nodes, topo8] = Read_mesh_ABAQUS(inpFile);
model = usemBuildModel(nodes, topo8(:, 1:4), pdeg);
[xi, eta] = cyl_local_coords(model.p);
nb = 4 + 4*(pdeg - 1);                      % corner and edge nodes
N8 = ShapeFunctions2D_8nodes(xi(1:nb).', eta(1:nb).');
for e = 1:model.num_Elem
    ids = model.Topology_Node(e, 1:nb);
    model.xCoor(ids) = N8.'*nodes(topo8(e, :), 1);
    model.yCoor(ids) = N8.'*nodes(topo8(e, :), 2);
end
x = model.xCoor;  y = model.yCoor;  tol = 1e-12;
roller = find(abs(x) < tol & abs(abs(y) - 1) < tol);
pin = find(abs(x) < tol & abs(y) < tol);
model.A = find(abs(x - 2) < tol & abs(y) < tol);
model.fixed = sort([2*roller - 1; 2*pin - 1; 2*pin]);
end
