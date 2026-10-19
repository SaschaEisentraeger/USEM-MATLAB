function K = cyl_Ke_PP(p, C, x, y, nq)
%CYL_KE_PP Conventional (Bubnov-Galerkin) spectral element stiffness.
%   Parametric shape functions as test and trial functions (PP), Gauss
%   quadrature with nq x nq points. Same formulation as
%   tests/helpers/usemKeSymmetric (which uses nq = 2p), vectorized over the
%   integration points.
ndof = 2*numel(x);
Qd = cyl_quad(p, nq);
W = Qd.W.';  dNxi = Qd.dNxi;  dNeta = Qd.dNeta;
[J11, J12, J21, J22, detJ] = cyl_jacobian(dNxi, dNeta, x, y);
J11 = J11(:).';  J12 = J12(:).';  J21 = J21(:).';  J22 = J22(:).';  detJ = detJ(:).';
dNx = ( dNxi.*J22 - dNeta.*J12)./detJ;
dNy = (-dNxi.*J21 + dNeta.*J11)./detJ;
Wd = W.*detJ;
A = @(a, b) (a.*Wd)*b.';
c = C;
Kxx = c(1,1)*A(dNx, dNx) + c(1,3)*A(dNx, dNy) + c(3,1)*A(dNy, dNx) + c(3,3)*A(dNy, dNy);
Kxy = c(1,2)*A(dNx, dNy) + c(1,3)*A(dNx, dNx) + c(3,2)*A(dNy, dNy) + c(3,3)*A(dNy, dNx);
Kyx = c(2,1)*A(dNy, dNx) + c(2,3)*A(dNy, dNy) + c(3,1)*A(dNx, dNx) + c(3,3)*A(dNx, dNy);
Kyy = c(2,2)*A(dNy, dNy) + c(2,3)*A(dNy, dNx) + c(3,2)*A(dNx, dNy) + c(3,3)*A(dNx, dNx);
K = zeros(ndof);
K(1:2:end, 1:2:end) = Kxx;  K(1:2:end, 2:2:end) = Kxy;
K(2:2:end, 1:2:end) = Kyx;  K(2:2:end, 2:2:end) = Kyy;
end
