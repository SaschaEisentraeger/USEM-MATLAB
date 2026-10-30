%% Elasticity matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 11.02.2010 (University of Magdeburg)
%
% Compute the plane stress elasticity matrix based on the Young's modulus (E) and
% Poisson's ratio (nu) as input parameters. Note that a linear elastic
% material behavior is assumed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% E:        Young's modulus
% nu:       Poisson's ratio
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% C:    	Elasticity (constitutive) matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [C] = HookeMatrix2d(E,nu)

C = zeros(3,3);
C(1,1) =                E/(1-nu^2);
C(1,2) = nu*            E/(1-nu^2);
C(2,1) = nu*            E/(1-nu^2);
C(2,2) =                E/(1-nu^2);
C(3,3) = (1-nu)/2*      E/(1-nu^2);

end
