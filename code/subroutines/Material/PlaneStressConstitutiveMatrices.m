%% Elasticity matrix - Plane stress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 29.03.2011 (University of Magdeburg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% C:                Three-dimensional elasticity matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% C_sigma:          Two-dimensional elasticity matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C_sigma] = PlaneStressConstitutiveMatrices(C)
 
% 2D constitutive matrices
C_sigma = zeros(3,3);

% Reduce 3D to 2D
C_sigma(1,1) = C(1,1)-C(1,3)*C(1,3)/C(3,3);
C_sigma(1,2) = C(1,2)-C(1,3)*C(2,3)/C(3,3);
C_sigma(2,1) = C(1,2)-C(1,3)*C(2,3)/C(3,3);
C_sigma(2,2) = C(2,2)-C(2,3)*C(2,3)/C(3,3);
C_sigma(3,3) = C(4,4);

end