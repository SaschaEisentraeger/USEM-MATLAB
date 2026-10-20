function dg = cyl_diagnostics(model, C, ng)
%CYL_DIAGNOSTICS Properties of the trial functions of a mesh.
%   For the two trial spaces of the unsymmetric element (ShapeFunctions2D_
%   Metric, method 1: skew coordinates, suffix PM; method 2: hybrid space,
%   suffix PH) and the conventional element (PP):
%   dg.lebesgue / lebesguePH : max over the elements of the Lebesgue
%                 constant max_x sum_i |M_i(x)|, sampled at the images of an
%                 ng x ng grid of the reference element
%   dg.lebesgueGLL : Lebesgue constant of the parametric (GLL) interpolation
%                 (the value of both trial spaces on an affine element)
%   dg.condV / condVPH : max over the elements of cond(V) of the generalized
%                 Vandermonde matrix that ShapeFunctions2D_Metric solves
%                 with (basis of the code base)
%   dg.condKePM / condKePH / condKePP : max over the elements of
%                 sigma_max/sigma_(2n-3) of the element stiffness matrix
%                 (the three rigid-body modes excluded)
%   A trial space that the code base rejects on an element (singular
%   interpolation or det J = 0) gives NaN for its quantities.
if nargin < 3, ng = 16; end
p = model.p;
g = linspace(-1, 1, ng);
[G1, G2] = meshgrid(g, g);
N = ShapeFunctions2D_SEM(G1(:), G2(:), p);
dg.lebesgueGLL = max(sum(abs(N), 1));
names = {'', 'PH'};
for k = 1:2
    dg.(['lebesgue' names{k}]) = 0;
    dg.(['condV' names{k}]) = 0;
end
dg.condKePM = 0;  dg.condKePH = 0;  dg.condKePP = 0;
ok = [true true];
for e = 1:model.num_Elem
    ids = model.Topology_Node(e, :).';
    x = model.xCoor(ids);  y = model.yCoor(ids);
    for k = 1:2
        if ~ok(k), continue, end
        try
            [M, ~, ~, V] = ShapeFunctions2D_Metric(N.'*x, N.'*y, x, y, p, k, [G1(:), G2(:)]);
            q = p;  q.metricMethod = k;
            s = svd(cyl_Ke_PM(q, C, x, y, 2*p.pxi));
        catch err
            if ~startsWith(err.identifier, 'ShapeFunctions2D_Metric:'), rethrow(err); end
            ok(k) = false;
            continue
        end
        dg.(['lebesgue' names{k}]) = max(dg.(['lebesgue' names{k}]), max(sum(abs(M), 1)));
        dg.(['condV' names{k}]) = max(dg.(['condV' names{k}]), cond(V));
        f = {'condKePM', 'condKePH'};
        dg.(f{k}) = max(dg.(f{k}), s(1)/s(end-3));
    end
    s = svd(cyl_Ke_PP(p, C, x, y, 2*p.pxi));
    dg.condKePP = max(dg.condKePP, s(1)/s(end-3));
end
if ~ok(1), dg.lebesgue = NaN;  dg.condV = NaN;  dg.condKePM = NaN; end
if ~ok(2), dg.lebesguePH = NaN;  dg.condVPH = NaN;  dg.condKePH = NaN; end
end
