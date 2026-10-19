%% Elasticity matrix (3D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 11.02.2010 (University of Magdeburg)
%
% Three-dimensional elasticity matrix for a linear elastic material
% determined by Young's modulus and Poisson's ratio.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% E:        Young's modulus
% nu:       Poisson's ration
% 
% Shear modulus
% G = E/2*(1+nu)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% C:    Elasticity matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [C] = HookeMatrix3d(E,nu)

C = zeros(6,6);
C(1,1) = (1-nu)*    	E/((1+nu)*(1-2*nu));
C(1,2) = nu*            E/((1+nu)*(1-2*nu));
C(1,3) = nu*            E/((1+nu)*(1-2*nu));
C(2,1) = nu*            E/((1+nu)*(1-2*nu));
C(2,2) = (1-nu)*     	E/((1+nu)*(1-2*nu));
C(2,3) = nu*            E/((1+nu)*(1-2*nu));
C(3,1) = nu*            E/((1+nu)*(1-2*nu));
C(3,2) = nu*            E/((1+nu)*(1-2*nu));
C(3,3) = (1-nu)*        E/((1+nu)*(1-2*nu));
C(4,4) = (1-2*nu)/2*    E/((1+nu)*(1-2*nu));
C(5,5) = (1-2*nu)/2*    E/((1+nu)*(1-2*nu));
C(6,6) = (1-2*nu)/2*    E/((1+nu)*(1-2*nu));
    
end