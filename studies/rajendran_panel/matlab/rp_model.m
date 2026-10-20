function model = rp_model(inpFile, pdeg, center)
%RP_MODEL Spectral element model of the panel from a gmsh/ABAQUS file.
%   center = 'serendipity' : code base (usemReadAbaqus ->
%                            pSpectralElement_NodesTopology): all nodes at
%                            the images of the reference nodes under the
%                            8-node (serendipity) map, i.e. the exact
%                            isoparametric representation of the 8-node mesh.
%                            Stops (coincidentNodes) if an interior node
%                            falls onto an edge node (e.g. Mesh b, even p).
%          | 'bilinear'    : a different geometry: edge nodes on the
%                            parabolic 8-node edges as above, interior nodes
%                            from the bilinear map of the four corner nodes;
%                            for p = 2 the 9th node lies at the centroid of
%                            the corners, as in PM-QUAD9* of Rajendran 2010,
%                            Sec. 5.7. Built with
%                            pSpectralElement_NodesTopology on the 4-node
%                            corner mesh, the edge nodes are then moved onto
%                            the 8-node edges.
if nargin < 3, center = 'serendipity'; end
p = struct('pxi', pdeg, 'peta', pdeg, 'ndist', 2);
if strcmp(center, 'serendipity')
    tmp = tempname;  mkdir(tmp);  old = cd(tmp);
    try
        evalc('S = usemReadAbaqus(inpFile, p);');
    catch err
        cd(old);  rmdir(tmp, 's');  rethrow(err);
    end
    cd(old);  rmdir(tmp, 's');
    model = S;
else
    tmp = tempname;  mkdir(tmp);  old = cd(tmp);
    try
        [~, name, ext] = fileparts(inpFile);
        copyfile(inpFile, [name ext]);
        Adjust_Set_names_input_file_ABAQUS([name ext]);
        evalc('[nodes, topo8] = Read_mesh_ABAQUS([name ''_adjusted'' ext]);');
    catch err
        cd(old);  rmdir(tmp, 's');  rethrow(err);
    end
    cd(old);  rmdir(tmp, 's');
    model = usemBuildModel(nodes, topo8(:, 1:4), pdeg);
    [xi, eta] = cyl_local_coords(p);
    nb = 4 + 4*(pdeg - 1);                      % corner and edge nodes
    N8 = ShapeFunctions2D_8nodes(xi(1:nb).', eta(1:nb).');
    for e = 1:model.num_Elem
        ids = model.Topology_Node(e, 1:nb);
        model.xCoor(ids) = N8.'*nodes(topo8(e, :), 1);
        model.yCoor(ids) = N8.'*nodes(topo8(e, :), 2);
    end
end
x = model.xCoor;  y = model.yCoor;  tol = 1e-12;
model.roller = find(abs(x) < tol & abs(abs(y) - 1) < tol);
model.pin = find(abs(x) < tol & abs(y) < tol);
model.A = find(abs(x - 2) < tol & abs(y) < tol);
% Constrained DOFs from the ABAQUS file must agree with the geometric sets
fixed = sort([2*model.roller - 1; 2*model.pin - 1; 2*model.pin]);
model.fixed = fixed;
end
