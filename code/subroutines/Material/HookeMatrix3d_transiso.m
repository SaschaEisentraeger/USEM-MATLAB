%% Elasticity matrix for a transversely isotropic material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 17.05.2010 (University of Magdeburg)
%
% Compute the elasticity matrix for a transversely isotropic material. The
% direction of the fibers is aligned with the global x-direction.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% EL:       Young's modulus in the direction of the fibers (longitudinal)
% ET:       Young's modulus in orthogonal to the fiber direction
% GLT:      Shear modulus in the direction of the fibers
% nuLT:     Poisson's ratio for loading in fiber direction
% nuTT:     Poisson's ratio for loading orthogonal to the fiber direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% C:        Elasticity matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C] = HookeMatrix3d_transiso(EL,ET,GLT,nuLT,nuTT)

C = zeros(6,6);
C(1,1) = EL*(1-nuTT)/(1-nuTT-2*(ET/EL)*nuLT^2);
C(1,2) = ET*(nuLT)/(1-nuTT-2*(ET/EL)*nuLT^2);
C(1,3) = C(1,2);
C(2,1) = C(1,2);
C(2,2) = ET*(1-(ET/EL)*nuLT^2)/((1-nuTT-2*(ET/EL)*nuLT^2)*(1+nuTT));
C(2,3) = ET*(nuTT+(ET/EL)*nuLT^2)/((1-nuTT-2*(ET/EL)*nuLT^2)*(1+nuTT));
C(3,1) = C(1,2);
C(3,2) = C(2,3);
C(3,3) = C(2,2);
C(4,4) = GLT;
C(5,5) = 1/2*(C(2,2)-C(2,3));
C(6,6) = C(4,4);

end