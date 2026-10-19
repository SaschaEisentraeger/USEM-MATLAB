function K = usemKeSymmetric(p, C, x, y)
%USEMKESYMMETRIC Classical (Bubnov-Galerkin) spectral element stiffness.
%
%   K = usemKeSymmetric(p, C, x, y) uses the parametric shape functions
%   both as test and as trial functions. It is NOT part of the code base;
%   the tests use it as a reference to show that a mesh is distorted
%   enough to spoil the conventional element, while the unsymmetric
%   element (Ke_2D) stays exact. Same quadrature as Ke_2D (2p points).

ndof = 2*numel(x);
K = zeros(ndof);
[xi, wxi]   = lgwt(2*p.pxi,  -1, 1);
[eta, weta] = lgwt(2*p.peta, -1, 1);
[Xi, Eta] = meshgrid(xi, eta);
[Wxi, Weta] = meshgrid(wxi, weta);
IP = [Xi(:), Eta(:)];
IW = Wxi(:).*Weta(:);

[~, dNxi, dNeta] = ShapeFunctions2D_SEM(IP(:,1), IP(:,2), p);
for i = 1:numel(IW)
    J = Jacobi_Matrix(p, IP(i,:), x, y);
    dN = J\[dNxi(:,i)'; dNeta(:,i)'];
    B = zeros(3, ndof);
    B(1, 1:2:end) = dN(1,:);
    B(2, 2:2:end) = dN(2,:);
    B(3, 1:2:end) = dN(2,:);
    B(3, 2:2:end) = dN(1,:);
    K = K + B'*C*B*det(J)*IW(i);
end

end
