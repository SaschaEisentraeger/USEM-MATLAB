function F = cyl_loads(model, prob, nq)
%CYL_LOADS Consistent load vector (parametric test functions).
%   Body force  : F_i = int N_i f det(J) dxi deta  (signed det J)
%   Traction    : F_i = int N_i (sigma_exact n_h) ds on the arcs r = a, b
%                 (n_h: outward normal of the discrete, degree-p boundary,
%                 so that the problem posed on the discrete domain has the
%                 exact field as solution). On the straight edges only
%                 Dirichlet conditions are prescribed.
p = model.p;
if nargin < 3 || isempty(nq), nq = 3*p.pxi + 2; end
T = model.Topology_Node;
F = zeros(model.num_DoF, 1);
Qd = cyl_quad(p, nq);
wq = Qd.wq;  W = Qd.W;  N = Qd.N;  dNxi = Qd.dNxi;  dNeta = Qd.dNeta;
for e = 1:size(T, 1)
    ids = T(e, :).';
    x = model.xCoor(ids);  y = model.yCoor(ids);
    [~, ~, ~, ~, detJ] = cyl_jacobian(dNxi, dNeta, x, y);
    xg = (N.'*x);  yg = (N.'*y);
    w = W(:).*detJ(:);
    F(2*ids-1) = F(2*ids-1) + N*(prob.fx(xg, yg).*w);
    F(2*ids)   = F(2*ids)   + N*(prob.fy(xg, yg).*w);
end
% Tractions on the circular arcs
for k = 1:size(model.arcEdges, 1)
    e = model.arcEdges(k, 1);  edge = model.arcEdges(k, 2);  R = model.arcEdges(k, 3);
    Ne = Qd.edge{edge}.N;  dNxi_e = Qd.edge{edge}.dNxi;  dNeta_e = Qd.edge{edge}.dNeta;
    ids = T(e, :).';
    x = model.xCoor(ids);  y = model.yCoor(ids);
    [J11, J12, J21, J22] = cyl_jacobian(dNxi_e, dNeta_e, x, y);
    % Tangent of the discrete (degree-p) edge and its normal n_h
    if mod(edge, 2) == 1, tx = J11; ty = J12; else, tx = J21; ty = J22; end
    tx = tx(:);  ty = ty(:);
    ds = hypot(tx, ty);
    xg = Ne.'*x;  yg = Ne.'*y;
    nx = ty./ds;  ny = -tx./ds;
    % Orientation: outward normal -e_r on r = a, +e_r on r = b
    sgn = sign(nx.*xg + ny.*yg).*(1 - 2*(R < 2));
    nx = sgn.*nx;  ny = sgn.*ny;
    txv = prob.sxx(xg, yg).*nx + prob.sxy(xg, yg).*ny;
    tyv = prob.sxy(xg, yg).*nx + prob.syy(xg, yg).*ny;
    w = wq.*ds;
    F(2*ids-1) = F(2*ids-1) + Ne*(txv.*w);
    F(2*ids)   = F(2*ids)   + Ne*(tyv.*w);
end
end
