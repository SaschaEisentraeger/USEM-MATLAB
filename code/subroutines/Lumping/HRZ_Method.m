%% Diagonalization of the mass matrix by the HRZ algorithm
% 10.04.2017
% University of Magdeburg (OvGU)
% Sascha Duczek
%
% The HRZ algorithm is named after Hinton, Rock and Zienkiewicz who first
% proposed this method (often also referred to as diagonal scaling).
%
% Hinton, Rock, Zienkiewicz 1976 : "Note on Mass Lumping and Related 
% Processes in the Finite Element Method"
%
% Hughes "The Finite Element Method: Linear Static and Dynamic Finite 
% Element Analysis" p. 436
%-------------------------------------------------------------------------%
% Input
% M:        Global mass matrix
%-------------------------------------------------------------------------%
% Output
% Mdiag:	Lumped mass matrix
%-------------------------------------------------------------------------%
function [Mdiag] = HRZ_Method(M)

% Sum of the diagonal components of the consistent mass matrix
beta = sum(diag(M));

% Total mass of the finite element
m = sum(sum(M));

% Scale the components on the main diagonal of the consistent mass matrix 
% by the ratio m/beta
m_vec = m/beta*diag(M);

% Diagonal mass matrix;
Mdiag = diag(m_vec);

end