function dg = cyl_diagnostics(model, C, ng)
%CYL_DIAGNOSTICS Properties of the metric trial functions of a mesh.
%   dg.lebesgue : max over the elements of the Lebesgue constant of the
%                 metric interpolation, max_x sum_i |M_i(x)|, sampled at the
%                 images of an ng x ng grid of the reference element
%                 (for comparison: the GLL tensor Lebesgue constant is the
%                 value on an affine element)
%   dg.lebesgueGLL : Lebesgue constant of the parametric (GLL) interpolation
%   dg.condV    : max over the elements of cond(V) of the generalized
%                 Vandermonde matrix in ShapeFunctions2D_Metric
%   dg.condKePM, dg.condKePP : max over the elements of
%                 sigma_max/sigma_(2n-3) of the element stiffness matrix
%                 (the three rigid-body modes excluded), unsymmetric (PM)
%                 and conventional (PP) element
if nargin < 3, ng = 31; end
p = model.p;
g = linspace(-1, 1, ng);
[G1, G2] = meshgrid(g, g);
N = ShapeFunctions2D_SEM(G1(:), G2(:), p);
dg.lebesgueGLL = max(sum(abs(N), 1));
dg.lebesgue = 0;  dg.condV = 0;  dg.condKePM = 0;  dg.condKePP = 0;
for e = 1:model.num_Elem
    ids = model.Topology_Node(e, :).';
    x = model.xCoor(ids);  y = model.yCoor(ids);
    M = ShapeFunctions2D_Metric(N.'*x, N.'*y, x, y, p);
    dg.lebesgue = max(dg.lebesgue, max(sum(abs(M), 1)));
    [xs, ys] = CorotationalSkewCoordinates(x, y, x, y, p);
    V = zeros(numel(x));  c = 0;
    for i = 0:p.peta
        for j = 0:p.pxi
            c = c + 1;  V(c, :) = (xs.^j).*(ys.^i);
        end
    end
    dg.condV = max(dg.condV, cond(V));
    s = svd(cyl_Ke_PM(p, C, x, y, 2*p.pxi));
    dg.condKePM = max(dg.condKePM, s(1)/s(end-3));
    s = svd(cyl_Ke_PP(p, C, x, y, 2*p.pxi));
    dg.condKePP = max(dg.condKePP, s(1)/s(end-3));
end
end
