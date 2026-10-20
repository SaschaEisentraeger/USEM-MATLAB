function model = rp_model(inpFile, pdeg, center)
%RP_MODEL Spectral element model of the panel from a gmsh/ABAQUS file.
%   The 8-node mesh is lifted to degree pdeg with the code base
%   (usemReadAbaqus -> pSpectralElement_NodesTopology).
%   center = 'serendipity' (code base: interior nodes from the 8-node map)
%          | 'bilinear'    (interior nodes from the bilinear map of the four
%                           corner nodes; for p = 2 the 9th node lies at the
%                           centroid of the corners, as in PM-QUAD9* of
%                           Rajendran 2010, Sec. 5.7)
%   Edge nodes are those of the code base in both cases. With 'serendipity',
%   interior nodes of strongly distorted elements can coincide with edge
%   nodes (e.g. p = 2, Mesh b), which makes the metric interpolation
%   singular.
if nargin < 3, center = 'serendipity'; end
p = struct('pxi', pdeg, 'peta', pdeg, 'ndist', 2);
tmp = tempname;  mkdir(tmp);  old = cd(tmp);
try
    evalc('S = usemReadAbaqus(inpFile, p);');
catch err
    cd(old);  rmdir(tmp, 's');  rethrow(err);
end
cd(old);  rmdir(tmp, 's');
model = S;
if strcmp(center, 'bilinear')
    [xi, eta] = cyl_local_coords(p);
    nb = 4 + 4*(pdeg - 1);                      % corner and edge nodes
    N4 = ShapeFunctions2D_4nodes(xi(nb+1:end).', eta(nb+1:end).');
    for e = 1:S.num_Elem
        c = S.Topology_Node(e, 1:4);
        m = S.Topology_Node(e, nb+1:end);
        model.xCoor(m) = N4.'*S.xCoor(c(:));
        model.yCoor(m) = N4.'*S.yCoor(c(:));
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
