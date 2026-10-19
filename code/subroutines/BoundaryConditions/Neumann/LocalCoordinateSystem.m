%% Local coordinate system (1D reference element)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 02.09.2016 (University of Magdeburg)
%
% Compute the local coordinate system at a point P on the edge of a
% two-dimensional spectral element.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% xi:       (Local) nodal coordinate of the origin of the coordinate system
% x:        x-coordinates of the element edge
% y:        y-coordinates of the element edge
% ndist:    Nodal distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% V1:       Tangent vector
% V2:       Normal vector
% V2norm:  	Unit normal vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V1,V2,V2norm] = LocalCoordinateSystem(xi,x,y,ndist)

%-------------------------------------------------------------------------%
% Test uniqueness of coordinates --> Determine whether the element is
% parallel to one of the global coordinate planes. Adjust the nodal
% coordinates if necessary (round-off errors in the computation of the
% high order nodes lead to small force components that should not exist).
xU = unique(x);
yU = unique(y);
% Difference between the unique coordinates
diff_xU = diff(xU)./xU(2:end);
diff_yU = diff(yU)./yU(2:end);
% If the difference is larger than a prescribed error threshold the nodal
% coordinates are adjusted
err_thr = 1e-12;
idx = find(diff_xU < err_thr);
idy = find(diff_yU < err_thr);
% Adjust nodal coordinates
% x-coordinates
if ~isempty(idx)
    id = setdiff(1:numel(xU),idx+1);
    xU_New = xU(id);
    % Number of nodes
    nN = numel(xU_New);
    % Change nodal coordinates
    for i = 1:nN
        x(abs(x-xU_New(i)) < err_thr) = xU_New(i);
    end
end
% y-coordinates
if ~isempty(idy)
    id = setdiff(1:numel(yU),idy+1);
    yU_New = yU(id);
    % Number of nodes
    nN = numel(yU_New);
    % Change nodal coordinates
    for i = 1:nN
        y(abs(y-yU_New(i)) < err_thr) = yU_New(i);
    end
end
%-------------------------------------------------------------------------%

% Direction vectors
V1 = zeros(2,1);
V2 = zeros(2,1);

% Number of nodes
nNodes = length(x);

% Polynomial degree
p = nNodes - 1;

% Derivative of the one-dimensional shape functions
[~,dN] = ShapeFunctions1D_SEM(p,xi,ndist);

% Tangent vector: tangential to the element edge at the considered point
V1(1) = dot(dN,x);
V1(2) = dot(dN,y);

% Round-off error
V1(abs(V1) < round(max(V1)*eps*1e15)/1e12) = 0;

% Vector normal to V1
V2(1) =  V1(2);
V2(2) = -V1(1);
V2norm = V2/norm(V2);

end