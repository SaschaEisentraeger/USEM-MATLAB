function model = cyl_model(inpFile, pdeg, ndist)
%CYL_MODEL Spectral element model of the thick cylinder from a gmsh mesh.
%
%   model = cyl_model(inpFile, pdeg, ndist) reads the ABAQUS input file
%   (8-node elements) and lifts it to spectral elements of degree pdeg with
%   the routines of the code base (usemReadAbaqus ->
%   pSpectralElement_NodesTopology). The lifted edges on the circular arcs
%   r = a and r = b are parabolas through the three CPS8 nodes; they are
%   replaced by the exact arcs (interpolated at the GLL points of the edge)
%   by a Gordon-Hall (linear blending) correction. Element edges in the
%   interior of the domain are not modified, so the mesh stays conforming.
%
%   Additional fields: nodes8, topo8 (input mesh), arcEdges [elem, edge, R],
%   set.xsym / ysym / inner / outer (node numbers), minDetJ, negDetJ (number
%   of elements with det J <= 0 somewhere), area (signed, should be 2*pi).

if nargin < 3, ndist = 2; end
A_IN = 1;  B_OUT = 3;

p = struct('pxi', pdeg, 'peta', pdeg, 'ndist', ndist);
tmp = tempname;  mkdir(tmp);  old = cd(tmp);
try
    evalc('S = usemReadAbaqus(inpFile, p);');
catch err
    cd(old);  rmdir(tmp, 's');  rethrow(err);
end
cd(old);  rmdir(tmp, 's');

model = S;
model.nodes8 = S.nodes;
model.topo8 = S.topology;
x = S.xCoor;  y = S.yCoor;
[xiN, etaN] = cyl_local_coords(p);
pairs = [1 2; 2 3; 4 3; 1 4];         % edge k runs from pairs(k,1) to pairs(k,2)
model.arcEdges = zeros(0, 3);
for e = 1:S.num_Elem
    c = S.topology(e, :);
    for k = 1:4
        Xs = S.nodes(c(pairs(k,1)), 1:2);
        Xe = S.nodes(c(pairs(k,2)), 1:2);
        Xm = S.nodes(c(4+k), 1:2);
        R = [];
        for Rc = [A_IN, B_OUT]
            if abs(norm(Xs) - Rc) < 1e-10 && abs(norm(Xe) - Rc) < 1e-10
                R = Rc;
            end
        end
        if isempty(R), continue, end
        model.arcEdges(end+1, :) = [e, k, R];
        if mod(k, 2) == 1, s = xiN; else, s = etaN; end
        switch k
            case 1, w = (1 - etaN)/2;
            case 2, w = (1 + xiN)/2;
            case 3, w = (1 + etaN)/2;
            case 4, w = (1 - xiN)/2;
        end
        Xq = Xs.*(s.*(s-1)/2) + Xm.*(1 - s.^2) + Xe.*(s.*(s+1)/2);
        ts = atan2(Xs(2), Xs(1));  te = atan2(Xe(2), Xe(1));
        th = (ts + te)/2 + s*(te - ts)/2;
        d = R*[cos(th), sin(th)] - Xq;
        ids = S.Topology_Node(e, :).';
        x(ids) = x(ids) + w.*d(:,1);
        y(ids) = y(ids) + w.*d(:,2);
    end
end
model.xCoor = x;  model.yCoor = y;

% Boundary node sets
r = hypot(x, y);  tol = 1e-10;
model.set.xsym  = find(abs(x) < tol);
model.set.ysym  = find(abs(y) < tol);
model.set.inner = find(abs(r - A_IN) < tol);
model.set.outer = find(abs(r - B_OUT) < tol);

% Jacobian check on a dense grid (plus the area as a geometry check)
g = linspace(-1, 1, 41);
[G1, G2] = meshgrid(g, g);
[~, dNxi, dNeta] = ShapeFunctions2D_SEM(G1(:), G2(:), p);
[q, wq] = lgwt(2*pdeg + 4, -1, 1);
[Q1, Q2] = meshgrid(q, q);  W = wq*wq.';
[~, dQxi, dQeta] = ShapeFunctions2D_SEM(Q1(:), Q2(:), p);
model.minDetJ = inf;  model.negDetJ = 0;  model.area = 0;
for e = 1:S.num_Elem
    ids = S.Topology_Node(e, :).';
    [~, ~, ~, ~, dJ] = cyl_jacobian(dNxi, dNeta, x(ids), y(ids));
    model.minDetJ = min(model.minDetJ, min(dJ));
    model.negDetJ = model.negDetJ + any(dJ <= 0);
    [~, ~, ~, ~, dJ] = cyl_jacobian(dQxi, dQeta, x(ids), y(ids));
    model.area = model.area + sum(dJ(:).*W(:));
end
end
