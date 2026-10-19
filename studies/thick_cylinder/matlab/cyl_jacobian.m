function [J11, J12, J21, J22, detJ] = cyl_jacobian(dNxi, dNeta, x, y)
%CYL_JACOBIAN Jacobian entries at several points (same convention as
%   Jacobi_Matrix: J = [dx/dxi, dy/dxi; dx/deta, dy/deta]).
J11 = (dNxi.'*x).';   J12 = (dNxi.'*y).';
J21 = (dNeta.'*x).';  J22 = (dNeta.'*y).';
detJ = J11.*J22 - J12.*J21;
end
