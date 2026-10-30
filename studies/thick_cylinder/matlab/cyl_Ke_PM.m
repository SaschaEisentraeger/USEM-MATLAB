function K = cyl_Ke_PM(p, C, x, y, nq)
%CYL_KE_PM Unsymmetric spectral element stiffness (vectorized Ke_2D).
%   Same formulation as code/subroutines/ElementMatrices/Ke_2D.m:
%       K = int B_L*' C B_R dxi deta,
%   B_L* from adj(J)*grad_xi(N) (parametric test functions, |J| cancelled),
%   B_R from the metric trial functions (ShapeFunctions2D_Metric; trial
%   space selected by p.metricMethod: 1 skew coordinates (default), 2
%   hybrid), nq x nq Gauss points. As in Ke_2D, the local coordinates of
%   the integration points are passed to ShapeFunctions2D_Metric. The loop
%   over the integration points is vectorized; the result agrees with
%   Ke_2D(p,C,x,y,nq,nq) to round-off.
nn = numel(x);  ndof = 2*nn;
Qd = cyl_quad(p, nq);
W = Qd.W.';  N = Qd.N;  dNxi = Qd.dNxi;  dNeta = Qd.dNeta;
[J11, J12, J21, J22] = cyl_jacobian(dNxi, dNeta, x, y);
J11 = J11(:).';  J12 = J12(:).';  J21 = J21(:).';  J22 = J22(:).';
% adj(J)*[dN/dxi; dN/deta]  (rows: nodes, columns: integration points)
dNx = dNxi.*J22 - dNeta.*J12;
dNy = -dNxi.*J21 + dNeta.*J11;
xg = N.'*x;  yg = N.'*y;
[~, dMx, dMy] = ShapeFunctions2D_Metric(xg, yg, x, y, p, [], [Qd.xi, Qd.eta]);
% K = sum_i BL_i' C BR_i w_i, assembled blockwise for (u_x, u_y)
c = C;
A = @(a, b) (a.*W)*b.';
Kxx = c(1,1)*A(dNx, dMx) + c(1,3)*A(dNx, dMy) + c(3,1)*A(dNy, dMx) + c(3,3)*A(dNy, dMy);
Kxy = c(1,2)*A(dNx, dMy) + c(1,3)*A(dNx, dMx) + c(3,2)*A(dNy, dMy) + c(3,3)*A(dNy, dMx);
Kyx = c(2,1)*A(dNy, dMx) + c(2,3)*A(dNy, dMy) + c(3,1)*A(dNx, dMx) + c(3,3)*A(dNx, dMy);
Kyy = c(2,2)*A(dNy, dMy) + c(2,3)*A(dNy, dMx) + c(3,2)*A(dNx, dMy) + c(3,3)*A(dNx, dMx);
K = zeros(ndof);
K(1:2:end, 1:2:end) = Kxx;  K(1:2:end, 2:2:end) = Kxy;
K(2:2:end, 1:2:end) = Kyx;  K(2:2:end, 2:2:end) = Kyy;
end
