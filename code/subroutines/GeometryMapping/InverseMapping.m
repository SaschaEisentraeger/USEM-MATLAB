%% Inverse mapping function: Global to local
%
% The global coordinates of a point (x,y) are given as input data. A
% nonlinear system of equations is solved to obtain the corresponding local 
% coordinates of the point (xi,eta).
% In this function we compute the global coordinates of the point from the
% iteratively determined local coordinates and compare them with the
% original global coordinates. From this residuum we know how accurate the
% local coordinates are and if the iteration process has to be continued.
% Therefore, the residuum is passed to the function fsolve that checks the
% stop criterion. See also the function: Global2LocalCoordinatesQuad.m
%
%-------------------------------------------------------------------------%
%% Input
%
% Coor:         Local coordinates of the point (as computed by fsolve)
% x:            x-coordinates of the element nodes (geometry)
% y:            y-coordinates of the element nodes (geometry)
% P:           	Global coordinates of the point
% p:          	Polynomial degree structure
% p.pxi:       	Polynomial degree in xi
% p.peta:      	Polynomial degree in eta
% p.ndist:    	Nodal distribution
%
%-------------------------------------------------------------------------%
%% Output
%
% F:            Residuum of the nonlinear system of equations
%
%-------------------------------------------------------------------------%
function F = InverseMapping(Coor,x,y,P,p)

% Matrix of shape functions
N = ShapeFunctions2D_SEM(Coor(:,1),Coor(:,2),p)';

% Residuum
F = [N*x-P(:,1), N*y-P(:,2)];

end