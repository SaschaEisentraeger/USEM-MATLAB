function K = rp_Ke(p, C, x, y, nq, form, method)
%RP_KE Element stiffness matrices of Rajendran (2010), Sec. 5, degree p.
%   form = 'PM' : parametric test / metric trial functions, |J| cancelled
%                 (PM-QUAD9*; Ke_2D of the code base for method = 1)
%          'MM' : metric functions on both sides, int B_R' D B_R det J
%          'PP' : parametric functions on both sides (conventional)
%   method (metric functions, ShapeFunctions2D_Metric):
%          0 = Cartesian coordinates relative to the centroid
%              (Rajendran 2010, Eq. (20)), 1 = corotational skew
%              coordinates (code-base default)
if strcmp(form, 'PP')
    K = cyl_Ke_PP(p, C, x, y, nq);
    return
end
nn = numel(x);
Qd = cyl_quad(p, nq);
W = Qd.W.';
[J11, J12, J21, J22, detJ] = cyl_jacobian(Qd.dNxi, Qd.dNeta, x, y);
J11 = J11(:).';  J12 = J12(:).';  J21 = J21(:).';  J22 = J22(:).';  detJ = detJ(:).';
xg = Qd.N.'*x;  yg = Qd.N.'*y;
[~, dMx, dMy] = ShapeFunctions2D_Metric(xg, yg, x, y, p, method);
if strcmp(form, 'PM')
    aX = Qd.dNxi.*J22 - Qd.dNeta.*J12;      % adj(J) grad_xi N
    aY = -Qd.dNxi.*J21 + Qd.dNeta.*J11;
    Wt = W;
else                                       % MM
    aX = dMx;  aY = dMy;
    Wt = W.*detJ;
end
A = @(a, b) (a.*Wt)*b.';
c = C;
Kxx = c(1,1)*A(aX, dMx) + c(1,3)*A(aX, dMy) + c(3,1)*A(aY, dMx) + c(3,3)*A(aY, dMy);
Kxy = c(1,2)*A(aX, dMy) + c(1,3)*A(aX, dMx) + c(3,2)*A(aY, dMy) + c(3,3)*A(aY, dMx);
Kyx = c(2,1)*A(aY, dMx) + c(2,3)*A(aY, dMy) + c(3,1)*A(aX, dMx) + c(3,3)*A(aX, dMy);
Kyy = c(2,2)*A(aY, dMy) + c(2,3)*A(aY, dMx) + c(3,2)*A(aX, dMy) + c(3,3)*A(aX, dMx);
K = zeros(2*nn);
K(1:2:end, 1:2:end) = Kxx;  K(1:2:end, 2:2:end) = Kxy;
K(2:2:end, 1:2:end) = Kyx;  K(2:2:end, 2:2:end) = Kyy;
end
