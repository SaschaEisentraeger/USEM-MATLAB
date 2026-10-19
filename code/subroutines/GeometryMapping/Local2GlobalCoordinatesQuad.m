% Mapping from the local to the global domain - Quadrilateral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The global coordinates of a point are computed from the local coordinates
% associated to the reference/master element; i.e. mapping from the local
% to the global domain.
%
% Sascha Duczek (25.08.2016)
%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% p:            Polynomial degree structure
% p.pxi:      	Polynomial degree in xi
% p.peta:      	Polynomial degree in eta
% p.ndist:   	Nodal distribution
% P:            Local coordinates of the point (that is mapped)
% NodalCoor:    Global coordinates of the element nodes (geometry
%               description)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% P_Global: 	Global coordinates of the point (that is mapped)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P_Global] = Local2GlobalCoordinatesQuad(p,P,NodalCoor)

% Map all points
% Shape functions
N = ShapeFunctions2D_SEM(P(:,1),P(:,2),p);
	
% Global coordinates
x = N'*NodalCoor(:,1);
y = N'*NodalCoor(:,2); 

% Assign results
P_Global = [x, y];

end