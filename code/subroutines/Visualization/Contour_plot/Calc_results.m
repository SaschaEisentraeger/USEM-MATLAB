%% Compute the element results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 27.11.2014 (University of Magdeburg)
%
% The results within a finite element are computed/interpolated based on
% the nodal displacement values and the shape functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% p:              	Polynomial degree structure
% p.pxi:          	Polynomial degree in xi
% p.peta:           Polynomial degree in eta
% p.ndist:          Nodal distribution
% xElem:            x-coordinates of all nodes of a finite element
% yElem:          	y-coordinates of all nodes of a finite element
% uElem:            Displacement vector  of a finite element
% ndiv:             Number of sub-divisions (in each direction) for the
%                   visualization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% x,y:              (Global) coordinates of the evaluation points
% ux,uy:           	(Result) displacement vectors in x- and y-directions
% EP:               (Local) coordinates of the evaluation points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

function [x,y,ux,uy,EP] = Calc_results(p,xElem,yElem,uElem,ndiv)

% Interval boundary
max_coor = 1;

% Sub-division of the element for visualization purposes
xi   = linspace(-max_coor,max_coor,ndiv(1));
eta  = linspace(-max_coor,max_coor,ndiv(2));

% Save evaluation points in a vector
[Xi,Eta] = meshgrid(xi,eta);
EP = [reshape(Xi,[],1),reshape(Eta,[],1)];

% Separate the displacement vector in x-, y-parts
Ux = uElem(1:2:end);
Uy = uElem(2:2:end);

%------------------------Parametric shape functions-----------------------%
% % Matrix of shape functions
% % Parametric shape functions (xi,eta)
% N = ShapeFunctions2D_SEM(EP(:,1),EP(:,2),p);
%--------------------------Metric shape functions-------------------------%
% Matrix of shape functions
% Metric shape functions (x,y)
IPg = Local2GlobalCoordinatesQuad(p,EP,[xElem,yElem]);
N = ShapeFunctions2D_Metric(IPg(:,1),IPg(:,2),xElem,yElem,p);
%-------------------------------------------------------------------------%

% Compute the displacements and global coordinates
ux = N'*Ux;
uy = N'*Uy;
x  = N'*xElem;
y  = N'*yElem;

end