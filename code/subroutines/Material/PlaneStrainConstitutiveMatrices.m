%% Elasticity matrix - Plane strain
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
function [C_epsilon] = PlaneStrainConstitutiveMatrices(C)

% D constitutive matrices
C_epsilon = zeros(3,3);

% Reduce 3D to 2D
C_epsilon(1,1) = C(1,1);
C_epsilon(1,2) = C(1,2);
C_epsilon(2,1) = C(1,2);
C_epsilon(2,2) = C(2,2);
C_epsilon(3,3) = C(4,4);

end