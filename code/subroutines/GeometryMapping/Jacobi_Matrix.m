%% Jacobi-Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 11.02.2010 (University of Magdeburg)
%
% Computation of the Jacobi matrix --> mapping between the local and global
% domain.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% p:      	Polynomial degree structure
% p.pxi:  	Polynomial degree in xi
% p.peta: 	Polynomial degree in eta
% GaussP:   Coordinates of the integration points
% xCoor:    x-coordinates of the nodes
% yCoor:    y-coordinates of the nodes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% J:        Jacobi-Matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J] = Jacobi_Matrix(p,GaussP,xCoor,yCoor)

% Derivatives of the shape functions
[~,dNxi,dNeta] = ShapeFunctions2D_SEM(GaussP(1),GaussP(2),p);

% Compute the components of the Jacobi-Matrix
J_11 = dot(dNxi,xCoor);
J_12 = dot(dNxi,yCoor);
J_21 = dot(dNeta,xCoor);
J_22 = dot(dNeta,yCoor);

% 2D Jacobi-matrix
J = [J_11, J_12; J_21, J_22];

end