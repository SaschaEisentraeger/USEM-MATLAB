%% Loading function (spatial variation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek, 28/06/2019, UNSW Sydney
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% Param:        Parameter of the loading function
% Points:       Nodal coordinates of the evaluation points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% LoadVal:   	Load value at the evaluation points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LoadVal] = SpatialLoadingFunction(Param,Points)

% Coordinates of the evaluation points
x = Points(:,1);
y = Points(:,2);

% Amplitude of the function
Amp = Param(1);

%% Examples
% User-defined function for the spatial variation of the load
%-------------------------------------------------------------------------%
% LoadVal = Amp*sin(pi/2*x).*sin(pi*y);
%-------------------------------------------------------------------------%
% Lee & Bathe: Effects of Element Distortions on the Performance of
% Isoparametric Elements, International Journal for Numerical Methods in
% Engineering, 1993
% Quadratic patch test:
%
% LoadVal = -Amp*120*y;
%
%-------------------------------------------------------------------------%
% Lee & Bathe: Effects of Element Distortions on the Performance of
% Isoparametric Elements, International Journal for Numerical Methods in
% Engineering, 1993
% Cubic patch test:
%
% Load identifier
Ident = Param(2);

if Ident == 1
	% Pressure load
    LoadVal = Amp*120*y;
elseif Ident == 2
    % Traction load
    LoadVal = Amp*(6 - 6*y.^2);
end
%
%-------------------------------------------------------------------------%

end