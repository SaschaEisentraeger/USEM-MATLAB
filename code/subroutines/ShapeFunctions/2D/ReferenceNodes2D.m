%% Local coordinates of the nodes of a spectral element
%-------------------------------------------------------------------------%
%
% Sascha Eisenträger, Otto von Guericke University, 20.10.2026
%
% Local coordinates (xi,eta) of the element nodes in the order used by
% ShapeFunctions2D_SEM: corners, interior nodes of edges 1-4, interior
% nodes row by row (xi fastest).
%
%-------------------------------------------------------------------------%
%% Input
%
% p:            Polynomial degree template
% p.pxi:     	Polynomial degree in xi
% p.peta:     	Polynomial degree in eta
% p.ndist:      Nodal distribution
%-------------------------------------------------------------------------%
%% Output
%
% xi, eta:      Local coordinates of the nodes (column vectors)
%-------------------------------------------------------------------------%
function [xi,eta] = ReferenceNodes2D(p)

x = nodal_distribution(p.pxi,p.ndist);
e = nodal_distribution(p.peta,p.ndist);
x = x(:);
e = e(:);

% Corners
xi  = [-1; 1; 1; -1];
eta = [-1; -1; 1; 1];

% Edges 1 (eta = -1), 2 (xi = 1), 3 (eta = 1), 4 (xi = -1)
xi  = [xi;  x(2:end-1);        ones(p.peta-1,1); x(2:end-1);       -ones(p.peta-1,1)];
eta = [eta; -ones(p.pxi-1,1);  e(2:end-1);       ones(p.pxi-1,1);  e(2:end-1)];

% Interior nodes
[XI,ETA] = meshgrid(x(2:end-1),e(2:end-1));
XI  = XI';
ETA = ETA';
xi  = [xi;  XI(:)];
eta = [eta; ETA(:)];

end
