%% Rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 17.05.2010 (University of Magdeburg)
%
% Rotation of the elasticity matrix about one of the global axis.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% C:        Elasticity matrix in the original coordinate system
% phi:      Rotation angle
% axis:		Rotation axis; 'x', 'y', or 'z'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% C:        Elasticity matrix in the rotated coordinate system
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C] = RotationHooke(C,phi,axis)

if strcmp(axis,'x')
%-------------------------------------------------------------------------%
% Rotation around the x-axis
if phi ~= 0
    T = zeros(6,6);
    T(1,1) = 1;
    T(2,2) = cos(phi)^2;
    T(2,3) = sin(phi)^2;
    T(2,5) = 2*cos(phi)*sin(phi);
    T(3,2) = sin(phi)^2;
    T(3,3) = cos(phi)^2;
    T(3,5) = -2*cos(phi)*sin(phi);
    T(4,4) = cos(phi);
    T(4,6) = sin(phi);
    T(5,2) = -cos(phi)*sin(phi);
    T(5,3) = cos(phi)*sin(phi);
    T(5,5) = cos(phi)^2-sin(phi)^2;
    T(6,4) = -sin(phi);
    T(6,6) = cos(phi);

    % Engineering-Tensor notation interchange matrix (gammy_ij = 2*epsilon_ij)
    R = zeros(6,6);
    R(1,1) = 1;
    R(2,2) = 1;
    R(3,3) = 1;
    R(4,4) = 2;
    R(5,5) = 2;
    R(6,6) = 2;
    
    C = T\C*R*T/R;
end
%-------------------------------------------------------------------------%
elseif strcmp(axis,'y')
% Rotation around the y-axis
if phi ~= 0
    T = zeros(6,6);
    T(1,1) = cos(phi)^2;
    T(1,3) = sin(phi)^2;
    T(1,6) = -2*cos(phi)*sin(phi);
    T(2,2) = 1;
    T(3,1) = sin(phi)^2;
    T(3,3) = cos(phi)^2;
    T(3,6) = 2*cos(phi)*sin(phi);
    T(4,4) = cos(phi);
    T(4,5) = -sin(phi);
    T(5,4) = sin(phi);
    T(5,5) = cos(phi);
    T(6,1) = cos(phi)*sin(phi);
    T(6,3) = -cos(phi)*sin(phi);
    T(6,6) = cos(phi)^2-sin(phi)^2;
    
    % Engineering-Tensor notation interchange matrix (gammy_ij = 2*epsilon_ij)
    R = zeros(6,6);
    R(1,1) = 1;
    R(2,2) = 1;
    R(3,3) = 1;
    R(4,4) = 2;
    R(5,5) = 2;
    R(6,6) = 2;
    
    C = T\C*R*T/R;
end
%-------------------------------------------------------------------------%
elseif strcmp(axis,'z')
% Rotation around the z-axis
if phi ~= 0
    T = zeros(6,6);
    T(1,1) = cos(phi)^2;
    T(1,2) = sin(phi)^2;
    T(1,4) = 2*sin(phi)*cos(phi);  
    T(2,1) = sin(phi)^2;
    T(2,2) = cos(phi)^2;
    T(2,4) = -2*sin(phi)*cos(phi);
    T(3,3) = 1;
    T(4,1) = -sin(phi)*cos(phi);
    T(4,2) = sin(phi)*cos(phi);
    T(4,4) = cos(phi)^2-sin(phi)^2;
    T(5,5) = cos(phi);
    T(5,6) = -sin(phi);
    T(6,5) = sin(phi);
    T(6,6) = cos(phi);

    % Engineering-Tensor notation interchange matrix (gammy_ij = 2*epsilon_ij)
    R = zeros(6,6);
    R(1,1) = 1;
    R(2,2) = 1;
    R(3,3) = 1;
    R(4,4) = 2;
    R(5,5) = 2;
    R(6,6) = 2;
    
    C = T\C*R*T/R;
end
%-------------------------------------------------------------------------%
end

C(abs(C) < 1e-12*max(max(abs(C)))) = 0;

end