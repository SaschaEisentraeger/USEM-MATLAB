function [K, M] = ms_assemble_KM(model, C, form, family)
%MS_ASSEMBLE_KM Dense global stiffness and consistent mass matrix for the
%classical (PP), Q_p(s,t) (PM) and hybrid (PH) spectral element.
%   [K, M] = ms_assemble_KM(model, C, form, family) integrates the element
%   matrices with the calibrated Gauss rule of the mesh family
%   (ms_rule/ms_quadrature_calibration). Without family, the former
%   default rule with 2p x 2p points is used.
P = model.p;
if nargin < 4 || isempty(family)
    nq = 2*P.pxi;
else
    nq = ms_rule(family, form, P.pxi);
end
T = Topology_DoF_2D(model.Topology_Node);
n = model.num_DoF;
K = zeros(n); M = zeros(n);
for e = 1:model.num_Elem
    nd = model.Topology_Node(e, :);
    x = model.xCoor(nd); y = model.yCoor(nd);
    K(T(e,:), T(e,:)) = K(T(e,:), T(e,:)) + ms_Ke(P, C, x, y, nq, form);
    M(T(e,:), T(e,:)) = M(T(e,:), T(e,:)) + ms_Me(P, x, y, nq, form);
end
end
