% Inverse mapping from the global to the local domain - Quadrilateral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the local coordinates from given global coordinates. The inverse
% mapping results in a nonlinear system of equations that has to be solved
% iteratively.
%
% Sascha Duczek (25.08.2016)
%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input variables
%
% P_Global:   	Global coordinates of the point
% P_Element:    Global coordinates of the nodes describing the geometry of
% 				a finite element
% p:            Polynomial degree structure
% p.pxi:       	Polynomial degree in xi
% p.peta:      	Polynomial degree in eta
% p.ndist:      Nodal distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output variables
%
% P_Local: 		Local coordinates of the mapped point
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P_Local] = Global2LocalCoordinatesQuad(P_Global,P_Element,p)

% Determine local coordinates
x = P_Element(:,1);
y = P_Element(:,2);
options = optimset('FunValCheck','on','Diagnostics','off', ...
                   'Display','off','TolX',1e-15,'TolFun',1e-15);

% Iterative solution of nonlinear system of equations for the local
% coordinates
[CoorLocal,~] = fsolve(@InverseMapping,-ones(size(P_Global)),options,...
                       x,y,P_Global,p);

% Assign results
xi  = CoorLocal(:,1);
eta = CoorLocal(:,2);

% Round the results to 12 significant digits
xi  = round(xi,12) ;
eta = round(eta,12);

P_Local = [xi eta];

end