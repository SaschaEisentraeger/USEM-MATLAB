%% Elasticity matrix for an orthotropic material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 17.05.2010 (University of Magdeburg)
%
% Compute the elasticity matrix for an orthotropic material.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% E1, E2, E3:       Young's moduli in the three directions of orthotropy
% G12, G13, G23:    Shear moduli in three orthogonal planes
% nu12, nu13, nu23: Poisson's ratios for the orthogonal planes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% C:                Elasticity matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C] = HookeMatrix3d_ortho(E,G,nu)
                               
% Young's moduli
E1 = E(1);
E2 = E(2);
E3 = E(3);

% Shear moduli
G12 = G(1);
G13 = G(2);
G23 = G(3);

% Poisson's ratios
nu12 = nu(1);
nu13 = nu(2);
nu23 = nu(3);
if numel(nu) == 3
    nu21 = E2/E1*nu12;
    nu31 = E3/E1*nu13;
    nu32 = E3/E2*nu23;
elseif numel(nu) == 6
	nu21 = nu(4);
    nu31 = nu(5);
    nu32 = nu(6);
end

% Auxiliary variable
D = 1 - nu12*nu21 - nu13*nu31 - nu23*nu32 - 2*nu12*nu23*nu31;

% Elasticity matrix
C = zeros(6,6);
C(1,1) = (1 - nu23*nu32)/D*E1;
C(1,2) = (nu13*nu32 + nu12)/D*E2;
C(1,3) = (nu12*nu23 + nu13)/D*E3;
C(2,2) = (1 - nu13*nu31)/D*E2;
C(2,3) = (nu21*nu13 + nu23)/D*E3;
C(3,3) = (1 - nu12*nu21)/D*E3;
C(4,4) = G12;
C(5,5) = G23;
C(6,6) = G13;

% Symmetry
C(2,1) = C(1,2);
C(3,1) = C(1,3);
C(3,2) = C(2,3);
    
end