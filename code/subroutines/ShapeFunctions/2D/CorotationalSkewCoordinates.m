%% Corotational skew coordinates
%-------------------------------------------------------------------------%
%
% Sascha Eisenträger, Otto von Guericke University, 17.01.2025
%
% Compute non-dimensional corotational skew coordinates.
%
% [1] Xie, Q. / Sze, K. Y. / Zhou, Y. X.:
% Modified and Trefftz Unsymmetric Finite Element Models, 2014,
% International Journal of Mechanics and Materials in Design, Vol. 12,
% No. 1, Springer Science and Business Media LLC, p. 53-70
%
%-------------------------------------------------------------------------%
%% Input
%
% (xQ,yQ):      Global coordinates of the query point(s)
% (xN,yN):      Global coordinates of the nodes
% p:            Polynomial degree template
% p.pxi:     	Polynomial degree in x
% p.peta:     	Polynomial degree in y
%-------------------------------------------------------------------------%
%% Output
%
% (xiQ,etaQ):   Skew corotational coordinates of the query point(s)
%-------------------------------------------------------------------------%
function [xiQ,etaQ,J] = CorotationalSkewCoordinates(xQ,yQ,xN,yN,p)

% Centroid of the element
xiEV = 0;
etaEV = 0;

% Global coordinates of the element centroid
N = ShapeFunctions2D_SEM(xiEV,etaEV,p);
xEV = N'*xN;
yEV = N'*yN;

% Jacobian matrix
[J] = Jacobi_Matrix(p,[xiEV,etaEV],xN,yN);
Jinv = [J(2,2), -J(1,2); -J(2,1), J(1,1)]/det(J);
JinvT = Jinv';

% Non-dimensional skew coordinates (query points)
CoorQ = JinvT*([xQ,yQ]' - [xEV,yEV]'); CoorQ = CoorQ';
xiQ = CoorQ(:,1);
etaQ = CoorQ(:,2);

end