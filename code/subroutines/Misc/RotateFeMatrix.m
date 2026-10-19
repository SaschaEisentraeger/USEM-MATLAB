%% Rotate element matrices (2D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 16.11.2016 (University of Magdeburg)
%
% Rotate the element matrices of a two-dimensional finite element.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% alpha:        Rotation angle
% Mat:          Matrix (that should be rotated)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% MatNew:    	Rotated matrix
% T:            Rotation matrix
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [MatNew,T] = RotateFeMatrix(alpha,Mat)

% Two-dimensional rotation matrix
T = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];

% Number of nodes
NumNodes = size(Mat,1)/2;

% Assemble the rotation matrix for the input Mat
Trot = zeros(size(Mat));
for i = 1:NumNodes
    Idx = 2*(i-1) + [1 2];
    Trot(Idx,Idx) = T;
end

% Compute the rotated matrix
MatNew = Trot'*Mat*Trot;

end