function C = OrthotropicStiffness(D)
%ORTHOTROPICSTIFFNESS Symmetric 6x6 stiffness matrix from ABAQUS orthotropic constants.

% Symmetric 6x6 matrix from the upper-triangle entries of the ABAQUS
% *Elastic, type=ORTHOTROPIC data line. ABAQUS orders the shear terms
% 12, 13, 23; this code uses the Voigt order 12, 23, 13 (rows 4, 5, 6),
% hence D1313 -> (6,6) and D2323 -> (5,5).
row = [1 1 2 1 2 3 4 6 5];
col = [1 2 2 3 3 3 4 6 5];
U = sparse(row,col,D);
C = full(U + U' - diag(diag(U)));
end
