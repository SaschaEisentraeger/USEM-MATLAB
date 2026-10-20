function Qd = cyl_quad(p, nq)
%CYL_QUAD Gauss rule and shape functions on the reference element (cached).
%   Qd.q, Qd.wq : 1D Gauss points/weights (nq), Qd.xi, Qd.eta, Qd.W : tensor
%   rule (column vectors), Qd.N, Qd.dNxi, Qd.dNeta : ShapeFunctions2D_SEM
%   at the tensor points; Qd.edge{k}.N, .dNxi, .dNeta at the Gauss points
%   of edge k (1: eta = -1, 2: xi = 1, 3: eta = 1, 4: xi = -1). Evaluating ShapeFunctions2D_SEM is expensive for
%   high p, so the results are cached per (p, ndist, nq).
persistent cache
if isempty(cache), cache = containers.Map(); end
key = sprintf('%d_%d_%d_%d', p.pxi, p.peta, p.ndist, nq);
if isKey(cache, key)
    Qd = cache(key);
    return
end
[q, wq] = lgwt(nq, -1, 1);
[Q1, Q2] = meshgrid(q, q);  W = wq*wq.';
Qd.q = q;  Qd.wq = wq;
Qd.xi = Q1(:);  Qd.eta = Q2(:);  Qd.W = W(:);
[Qd.N, Qd.dNxi, Qd.dNeta] = ShapeFunctions2D_SEM(Qd.xi, Qd.eta, p);
o = ones(nq, 1);
pts = {[q, -o], [o, q], [q, o], [-o, q]};
for k = 1:4
    [Ne, dNxe, dNee] = ShapeFunctions2D_SEM(pts{k}(:,1), pts{k}(:,2), p);
    Qd.edge{k} = struct('N', Ne, 'dNxi', dNxe, 'dNeta', dNee);
end
cache(key) = Qd;
end
