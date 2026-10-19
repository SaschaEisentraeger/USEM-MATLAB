%% Direction of the normal vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 02.09.2016 (University of Magdeburg)
%
% Determine the correct direction of the normal vector (inside/outside).
% This should normally follow from the definition of the element topology
% but is not always done consistently.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% n:                Unit normal vector
% xCoor:            x-coordinates of the nodes (of an element)
% yCoor:            y-coordinates of the nodes (of an element)
% EdgeNumber:       (Local) edge number
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% FlagDirection:    Variable indicating if the vector points inside (-1) or
%                   outside (+1) the element
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FlagDirection] = NormalVectorInsideOutside(n,xCoor,yCoor,...
                           EdgeNumber)

% Centroid of a quadrilateral element
xi  = 0;
eta = 0;
% Matrix of shape functions
N = ShapeFunctions2D_4nodes(xi,eta);
xc = dot(N,xCoor);
yc = dot(N,yCoor);

% Centroid of the edge (local coordinates)
if EdgeNumber == 1
    xi  =  0;
    eta = -1;
elseif EdgeNumber == 2
    xi  =  1;
    eta =  0;
elseif EdgeNumber == 3
    xi  =  0;
    eta =  1;
elseif EdgeNumber == 4
    xi  = -1;
    eta =  0;
end
% Matrix of shape functions
N = ShapeFunctions2D_4nodes(xi,eta);
% Centroid of the edge (global coordinates)
xfc = dot(N,xCoor);
yfc = dot(N,yCoor);

% Direction vector of the edge (only exact/constant if the edge is
% straight)
EdgeVector = [xfc - xc;...
              yfc - yc];

% Angle between edge direction vector and normal vector
Value = dot(n,EdgeVector);

% Sign indicating if the orientation of the normal vector has to be
% adjusted
if Value > 0
	FlagDirection = -1;
else
	FlagDirection =  1;
end

end