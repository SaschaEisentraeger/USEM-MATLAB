%% Displacement function (spatial variation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek, 28/06/2019, UNSW Sydney
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% Param:     	Parameter of the displacement function
% Points:       Nodal coordinates of the evaluation points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% DispVal:   	Displacement value at the evaluation points
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DispVal] = SpatialDisplacementFunction(Param,Points)

% Coordinates of the evaluation points
x = Points(:,1);
y = Points(:,2);

% Amplitude of the function
Amp = Param(1);

% Direction identifier: (1) ux; (2) uy
Ident = Param(2);

%% Examples
% User-defined function for the spatial variation of the displacement field
%-------------------------------------------------------------------------%
% DispVal = Amp*sin(pi*x).*sin(pi/2*y);
%-------------------------------------------------------------------------%
% Higher order patch test, version B: The displacement field is prescribed
% at the boundaries of the structure by an arbitrary (complete) polynomial
% function. Note that the function has to fulfill the equilibrium
% conditions div(sigma) = 0 (we assume no body load).
p = 2;
if p == 1
    % Linear polynomial (Note that each linear polynomial displacement
    % field fulfills the equilibrium conditions as second order
    % derivatives are involved).
    a = 1:3;
    b = 1:3;

    if Ident == 1
        % Displacement in x-direction
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y);
    elseif Ident == 2
        % Displacement in y-direction
        DispVal = Amp*(b(1) + b(2)*x + b(3)*y);
    else
        % Displacement in both directions
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y);
    end
elseif p == 2
    % Quadratic polynomial (Note that the coefficients have to be chosen
    % according to the equilibrium conditions! Suitable coefficients are
    % provided in: HigherOrderDisplacementField.m)
    a = 1:6;
    b = 1:6;
    nu = 0.3; % Poisson's ratio
    a(5) = (8*(nu - 4))/(nu + 1);
    b(5) = (4*(3*nu - 7))/(nu + 1);

    if Ident == 1
        % Displacement in x-direction
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2);
    elseif Ident == 2
        % Displacement in y-direction
        DispVal = Amp*(b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*x.*y + ...
                       b(6)*y.^2);
    else
        % Displacement in both directions
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2);
    end
elseif p == 3
    % Cubic polynomial (Note that the coefficients have to be chosen
    % according to the equilibrium conditions! Suitable coefficients are
    % provided in: HigherOrderDisplacementField.m)
    a = 1:10;
    b = 1:10;
    nu = 0.3; % Poisson's ratio
    a(5) =  (8*(nu - 4))/(nu + 1);
    b(5) =  (4*(3*nu - 7))/(nu + 1);
    a(8) =  (3*(7*nu - 13))/(nu + 3);
    b(8) = -(3*(17*nu - 3))/(2*nu);
    a(9) = -(3*(17*nu + 3))/(2*nu);
    b(9) =  (6*(5*nu - 2))/(nu + 3);

    if Ident == 1
        % Displacement in x-direction
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3);
    elseif Ident == 2
        % Displacement in y-direction
        DispVal = Amp*(b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*x.*y + ...
                       b(6)*y.^2 + b(7)*x.^3 + b(8)*x.^2.*y + ...
                       b(9)*x.*y.^2 + b(10)*y.^3);
    else
        % Displacement in both directions
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3);
    end
elseif p == 4
    % Quartic polynomial (Note that the coefficients have to be chosen
    % according to the equilibrium conditions! Suitable coefficients are
    % provided in: HigherOrderDisplacementField.m)
    a = 1:15;
    b = 1:15;
    nu = 0.3; % Poisson's ratio
    a(5)  =  (8*(nu - 4))/(nu + 1);
    b(5)  =  (4*(3*nu - 7))/(nu + 1);
    a(8)  =  (3*(7*nu - 13))/(nu + 3);
    b(8)  = -(3*(17*nu - 3))/(2*nu);
    a(9)  = -(3*(17*nu + 3))/(2*nu);
    b(9)  =  (6*(5*nu - 2))/(nu + 3);
    a(12) =  (4*(11*nu + 15))/(nu + 1);
    b(12) = -(4*(13*nu + 9))/(nu + 1);
    a(13) = -78;
    b(13) = -78;
    a(14) = -(4*(13*nu + 17))/(nu + 1);
    b(14) =  (4*(15*nu + 11))/(nu + 1);

    if Ident == 1
        % Displacement in x-direction
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3 + a(11)*x.^4 + ...
                       a(12)*x.^3.*y + a(13)*x.^2.*y.^2 + ...
                       a(14)*x.*y.^3 + a(15)*y.^4);
    elseif Ident == 2
        % Displacement in y-direction
        DispVal = Amp*(b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*x.*y + ...
                       b(6)*y.^2 + b(7)*x.^3 + b(8)*x.^2.*y + ...
                       b(9)*x.*y.^2 + b(10)*y.^3 + b(11)*x.^4 + ...
                       b(12)*x.^3.*y + b(13)*x.^2.*y.^2 + ...
                       b(14)*x.*y.^3 + b(15)*y.^4);
    else
        % Displacement in both directions
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3 + a(11)*x.^4 + ...
                       a(12)*x.^3.*y + a(13)*x.^2.*y.^2 + ...
                       a(14)*x.*y.^3 + a(15)*y.^4);
    end
elseif p == 5
    % Quintic polynomial (Note that the coefficients have to be chosen
    % according to the equilibrium conditions! Suitable coefficients are
    % provided in: HigherOrderDisplacementField.m)
    a = 1:21;
    b = 1:21;
    nu = 0.3; % Poisson's ratio
    a(5)  =  (8*(nu - 4))/(nu + 1);
    b(5)  =  (4*(3*nu - 7))/(nu + 1);
    a(8)  =  (3*(7*nu - 13))/(nu + 3);
    b(8)  = -(3*(17*nu - 3))/(2*nu);
    a(9)  = -(3*(17*nu + 3))/(2*nu);
    b(9)  =  (6*(5*nu - 2))/(nu + 3);
    a(12) =  (4*(11*nu + 15))/(nu + 1);
    b(12) = -(4*(13*nu + 9))/(nu + 1);
    a(13) = -78;
    b(13) = -78;
    a(14) = -(4*(13*nu + 17))/(nu + 1);
    b(14) =  (4*(15*nu + 11))/(nu + 1);
    a(17) =  (5*(16*nu + 37))/(nu + 2);
    b(17) = -(5*(11*nu + 53))/(3*nu + 1);
    a(18) = -(10*(11*nu - 21))/(3*nu + 1);
    b(18) = -(5*(37*nu + 69))/(nu + 2);
    a(19) = -(5*(37*nu + 79))/(nu + 2);
    b(19) = -(20*(13*nu - 8))/(3*nu + 1);
    a(20) = -(10*(13*nu + 29))/(3*nu + 1);
    b(20) =  (5*(21*nu + 37))/(nu + 2);

    if Ident == 1
        % Displacement in x-direction
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3 + a(11)*x.^4 + ...
                       a(12)*x.^3.*y + a(13)*x.^2.*y.^2 + ...
                       a(14)*x.*y.^3 + a(15)*y.^4 + a(16)*x.^5 + ...
                       a(17)*x.^4.*y + a(18)*x.^3.*y.^2 + ...
                       a(19)*x.^2.*y.^3 + a(20)*x.*y.^4 + a(21)*y.^5);
    elseif Ident == 2
        % Displacement in y-direction
        DispVal = Amp*(b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*x.*y + ...
                       b(6)*y.^2 + b(7)*x.^3 + b(8)*x.^2.*y + ...
                       b(9)*x.*y.^2 + b(10)*y.^3 + b(11)*x.^4 + ...
                       b(12)*x.^3.*y + b(13)*x.^2.*y.^2 + ...
                       b(14)*x.*y.^3 + b(15)*y.^4 + b(16)*x.^5 + ...
                       b(17)*x.^4.*y + b(18)*x.^3.*y.^2 + ...
                       b(19)*x.^2.*y.^3 + b(20)*x.*y.^4 + b(21)*y.^5);
    else
        % Displacement in both directions
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3 + a(11)*x.^4 + ...
                       a(12)*x.^3.*y + a(13)*x.^2.*y.^2 + ...
                       a(14)*x.*y.^3 + a(15)*y.^4 + a(16)*x.^5 + ...
                       a(17)*x.^4.*y + a(18)*x.^3.*y.^2 + ...
                       a(19)*x.^2.*y.^3 + a(20)*x.*y.^4 + a(21)*y.^5);
    end
elseif p == 6
    % 6th order polynomial (Note that the coefficients have to be chosen
    % according to the equilibrium conditions! Suitable coefficients are
    % provided in: HigherOrderDisplacementField.m)
    a = 1:28;
    b = 1:28;
    nu = 0.3; % Poisson's ratio
    a(5)  =  (8*(nu - 4))/(nu + 1);
    b(5)  =  (4*(3*nu - 7))/(nu + 1);
    a(8)  =  (3*(7*nu - 13))/(nu + 3);
    b(8)  = -(3*(17*nu - 3))/(2*nu);
    a(9)  = -(3*(17*nu + 3))/(2*nu);
    b(9)  =  (6*(5*nu - 2))/(nu + 3);
    a(12) =  (4*(11*nu + 15))/(nu + 1);
    b(12) = -(4*(13*nu + 9))/(nu + 1);
    a(13) = -78;
    b(13) = -78;
    a(14) = -(4*(13*nu + 17))/(nu + 1);
    b(14) =  (4*(15*nu + 11))/(nu + 1);
    a(17) =  (5*(16*nu + 37))/(nu + 2);
    b(17) = -(5*(11*nu + 53))/(3*nu + 1);
    a(18) = -(10*(11*nu - 21))/(3*nu + 1);
    b(18) = -(5*(37*nu + 69))/(nu + 2);
    a(19) = -(5*(37*nu + 79))/(nu + 2);
    b(19) = -(20*(13*nu - 8))/(3*nu + 1);
    a(20) = -(10*(13*nu + 29))/(3*nu + 1);
    b(20) =  (5*(21*nu + 37))/(nu + 2);
	a(23) =  (4*(33*nu - 17))/(nu + 1);
    b(23) = -(8*(4*nu + 29))/(nu + 1);
    a(24) = -80;
    b(24) = -80;
    a(25) = -(80*(4*nu - 21))/(3*(nu + 1));
    b(25) = -(40*(17*nu - 33))/(3*(nu + 1));
    a(26) = -170;
    b(26) = -170;
    a(27) = -(4*(17*nu + 67))/(nu + 1);
    b(27) =  (8*(21*nu - 4))/(nu + 1);

    if Ident == 1
        % Displacement in x-direction
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3 + a(11)*x.^4 + ...
                       a(12)*x.^3.*y + a(13)*x.^2.*y.^2 + ...
                       a(14)*x.*y.^3 + a(15)*y.^4 + a(16)*x.^5 + ...
                       a(17)*x.^4.*y + a(18)*x.^3.*y.^2 + ...
                       a(19)*x.^2.*y.^3 + a(20)*x.*y.^4 + a(21)*y.^5 + ...
                       a(22)*x.^6 + a(23)*x.^5.*y + a(24)*x.^4.*y.^2 + ...
                       a(25)*x.^3.*y.^3 + a(26)*x.^2.*y.^4 + ...
                       a(27)*x.*y.^5 + a(28)*y.^6);
    elseif Ident == 2
        % Displacement in y-direction
        DispVal = Amp*(b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*x.*y + ...
                       b(6)*y.^2 + b(7)*x.^3 + b(8)*x.^2.*y + ...
                       b(9)*x.*y.^2 + b(10)*y.^3 + b(11)*x.^4 + ...
                       b(12)*x.^3.*y + b(13)*x.^2.*y.^2 + ...
                       b(14)*x.*y.^3 + b(15)*y.^4 + b(16)*x.^5 + ...
                       b(17)*x.^4.*y + b(18)*x.^3.*y.^2 + ...
                       b(19)*x.^2.*y.^3 + b(20)*x.*y.^4 + b(21)*y.^5 + ...
                       b(22)*x.^6 + b(23)*x.^5.*y + b(24)*x.^4.*y.^2 + ...
                       b(25)*x.^3.*y.^3 + b(26)*x.^2.*y.^4 + ...
                       b(27)*x.*y.^5 + b(28)*y.^6);
    else
        % Displacement in both directions
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3 + a(11)*x.^4 + ...
                       a(12)*x.^3.*y + a(13)*x.^2.*y.^2 + ...
                       a(14)*x.*y.^3 + a(15)*y.^4 + a(16)*x.^5 + ...
                       a(17)*x.^4.*y + a(18)*x.^3.*y.^2 + ...
                       a(19)*x.^2.*y.^3 + a(20)*x.*y.^4 + a(21)*y.^5 + ...
                       a(22)*x.^6 + a(23)*x.^5.*y + a(24)*x.^4.*y.^2 + ...
                       a(25)*x.^3.*y.^3 + a(26)*x.^2.*y.^4 + ...
                       a(27)*x.*y.^5 + a(28)*y.^6);
    end
elseif p == 7
    % 7th order polynomial (Note that the coefficients have to be chosen
    % according to the equilibrium conditions! Suitable coefficients are
    % provided in: HigherOrderDisplacementField.m)
    a = 1:36;
    b = 1:36;
    nu = 0.3; % Poisson's ratio
    a(5)  =  (8*(nu - 4))/(nu + 1);
    b(5)  =  (4*(3*nu - 7))/(nu + 1);
    a(8)  =  (3*(7*nu - 13))/(nu + 3);
    b(8)  = -(3*(17*nu - 3))/(2*nu);
    a(9)  = -(3*(17*nu + 3))/(2*nu);
    b(9)  =  (6*(5*nu - 2))/(nu + 3);
    a(12) =  (4*(11*nu + 15))/(nu + 1);
    b(12) = -(4*(13*nu + 9))/(nu + 1);
    a(13) = -78;
    b(13) = -78;
    a(14) = -(4*(13*nu + 17))/(nu + 1);
    b(14) =  (4*(15*nu + 11))/(nu + 1);
    a(17) =  (5*(16*nu + 37))/(nu + 2);
    b(17) = -(5*(11*nu + 53))/(3*nu + 1);
    a(18) = -(10*(11*nu - 21))/(3*nu + 1);
    b(18) = -(5*(37*nu + 69))/(nu + 2);
    a(19) = -(5*(37*nu + 79))/(nu + 2);
    b(19) = -(20*(13*nu - 8))/(3*nu + 1);
    a(20) = -(10*(13*nu + 29))/(3*nu + 1);
    b(20) =  (5*(21*nu + 37))/(nu + 2);
	a(23) =  (4*(33*nu - 17))/(nu + 1);
    b(23) = -(8*(4*nu + 29))/(nu + 1);
    a(24) = -80;
    b(24) = -80;
    a(25) = -(80*(4*nu - 21))/(3*(nu + 1));
    b(25) = -(40*(17*nu - 33))/(3*(nu + 1));
    a(26) = -170;
    b(26) = -170;
    a(27) = -(4*(17*nu + 67))/(nu + 1);
    b(27) =  (8*(21*nu - 4))/(nu + 1);
    a(30) =  (21*(29*nu + 5))/(3*nu + 5);
    b(30) = -(21*(41*nu + 17))/(2*(2*nu + 1));
    a(31) = -(21*(123*nu + 65))/(2*(2*nu + 1));
    b(31) = -(42*(11*nu + 40))/(3*nu + 5);
    a(32) = -(70*(11*nu - 25))/(3*nu + 5);
    b(32) =  (35*(65*nu + 29))/(2*nu + 1);
    a(33) =  (35*(65*nu + 36))/(2*nu + 1);
    b(33) = -(35*(43*nu - 15))/(3*nu + 5);
    a(34) = -(21*(43*nu + 115))/(3*nu + 5);
    b(34) = -(21*(137*nu + 65))/(2*(2*nu + 1));
    a(35) = -(7*(137*nu + 79))/(2*(2*nu + 1));
    b(35) =  (14*(54*nu + 25))/(3*nu + 5);

    if Ident == 1
        % Displacement in x-direction
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3 + a(11)*x.^4 + ...
                       a(12)*x.^3.*y + a(13)*x.^2.*y.^2 + ...
                       a(14)*x.*y.^3 + a(15)*y.^4 + a(16)*x.^5 + ...
                       a(17)*x.^4.*y + a(18)*x.^3.*y.^2 + ...
                       a(19)*x.^2.*y.^3 + a(20)*x.*y.^4 + a(21)*y.^5 + ...
                       a(22)*x.^6 + a(23)*x.^5.*y + a(24)*x.^4.*y.^2 + ...
                       a(25)*x.^3.*y.^3 + a(26)*x.^2.*y.^4 + ...
                       a(27)*x.*y.^5 + a(28)*y.^6 + a(29)*x.^7 + ...
                       a(30)*x.^6.*y + a(31)*x.^5.*y.^2 + ...
                       a(32)*x.^4.*y.^3 + a(33)*x.^3.*y.^4 + ...
                       a(34)*x.^2.*y.^5 + a(35)*x.*y.^6 + a(36)*y.^7);
    elseif Ident == 2
        % Displacement in y-direction
        DispVal = Amp*(b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*x.*y + ...
                       b(6)*y.^2 + b(7)*x.^3 + b(8)*x.^2.*y + ...
                       b(9)*x.*y.^2 + b(10)*y.^3 + b(11)*x.^4 + ...
                       b(12)*x.^3.*y + b(13)*x.^2.*y.^2 + ...
                       b(14)*x.*y.^3 + b(15)*y.^4 + b(16)*x.^5 + ...
                       b(17)*x.^4.*y + b(18)*x.^3.*y.^2 + ...
                       b(19)*x.^2.*y.^3 + b(20)*x.*y.^4 + b(21)*y.^5 + ...
                       b(22)*x.^6 + b(23)*x.^5.*y + b(24)*x.^4.*y.^2 + ...
                       b(25)*x.^3.*y.^3 + b(26)*x.^2.*y.^4 + ...
                       b(27)*x.*y.^5 + b(28)*y.^6 + b(29)*x.^7 + ...
                       b(30)*x.^6.*y + b(31)*x.^5.*y.^2 + ...
                       b(32)*x.^4.*y.^3 + b(33)*x.^3.*y.^4 + ...
                       b(34)*x.^2.*y.^5 + b(35)*x.*y.^6 + b(36)*y.^7);
    else
        % Displacement in both directions
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3 + a(11)*x.^4 + ...
                       a(12)*x.^3.*y + a(13)*x.^2.*y.^2 + ...
                       a(14)*x.*y.^3 + a(15)*y.^4 + a(16)*x.^5 + ...
                       a(17)*x.^4.*y + a(18)*x.^3.*y.^2 + ...
                       a(19)*x.^2.*y.^3 + a(20)*x.*y.^4 + a(21)*y.^5 + ...
                       a(22)*x.^6 + a(23)*x.^5.*y + a(24)*x.^4.*y.^2 + ...
                       a(25)*x.^3.*y.^3 + a(26)*x.^2.*y.^4 + ...
                       a(27)*x.*y.^5 + a(28)*y.^6 + a(29)*x.^7 + ...
                       a(30)*x.^6.*y + a(31)*x.^5.*y.^2 + ...
                       a(32)*x.^4.*y.^3 + a(33)*x.^3.*y.^4 + ...
                       a(34)*x.^2.*y.^5 + a(35)*x.*y.^6 + a(36)*y.^7);
    end
elseif p == 8
    % 8th order polynomial (Note that the coefficients have to be chosen
    % according to the equilibrium conditions! Suitable coefficients are
    % provided in: HigherOrderDisplacementField.m)
    a = 1:45;
    b = 1:45;
    nu = 0.3; % Poisson's ratio
    a(5)  =  (8*(nu - 4))/(nu + 1);
    b(5)  =  (4*(3*nu - 7))/(nu + 1);
    a(8)  =  (3*(7*nu - 13))/(nu + 3);
    b(8)  = -(3*(17*nu - 3))/(2*nu);
    a(9)  = -(3*(17*nu + 3))/(2*nu);
    b(9)  =  (6*(5*nu - 2))/(nu + 3);
    a(12) =  (4*(11*nu + 15))/(nu + 1);
    b(12) = -(4*(13*nu + 9))/(nu + 1);
    a(13) = -78;
    b(13) = -78;
    a(14) = -(4*(13*nu + 17))/(nu + 1);
    b(14) =  (4*(15*nu + 11))/(nu + 1);
	a(17) =  (5*(16*nu + 37))/(nu + 2);
    b(17) = -(5*(11*nu + 53))/(3*nu + 1);
    a(18) = -(10*(11*nu - 21))/(3*nu + 1);
    b(18) = -(5*(37*nu + 69))/(nu + 2);
    a(19) = -(5*(37*nu + 79))/(nu + 2);
    b(19) = -(20*(13*nu - 8))/(3*nu + 1);
    a(20) = -(10*(13*nu + 29))/(3*nu + 1);
    b(20) =  (5*(21*nu + 37))/(nu + 2);
	a(23) =  (4*(33*nu - 17))/(nu + 1);
    b(23) = -(8*(4*nu + 29))/(nu + 1);
    a(24) = -80;
    b(24) = -80;
    a(25) = -(80*(4*nu - 21))/(3*(nu + 1));
    b(25) = -(40*(17*nu - 33))/(3*(nu + 1));
    a(26) = -170;
    b(26) = -170;
    a(27) = -(4*(17*nu + 67))/(nu + 1);
    b(27) =  (8*(21*nu - 4))/(nu + 1);
    a(30) =  (21*(29*nu + 5))/(3*nu + 5);
    b(30) = -(21*(41*nu + 17))/(2*(2*nu + 1));
    a(31) = -(21*(123*nu + 65))/(2*(2*nu + 1));
    b(31) = -(42*(11*nu + 40))/(3*nu + 5);
    a(32) = -(70*(11*nu - 25))/(3*nu + 5);
    b(32) =  (35*(65*nu + 29))/(2*nu + 1);
    a(33) =  (35*(65*nu + 36))/(2*nu + 1);
    b(33) = -(35*(43*nu - 15))/(3*nu + 5);
    a(34) = -(21*(43*nu + 115))/(3*nu + 5);
    b(34) = -(21*(137*nu + 65))/(2*(2*nu + 1));
    a(35) = -(7*(137*nu + 79))/(2*(2*nu + 1));
    b(35) =  (14*(54*nu + 25))/(3*nu + 5);
    a(38) =  (8*(37*nu + 41))/(nu + 1);
    b(38) = -(8*(39*nu + 35))/(nu + 1);
    a(39) = -1092;
    b(39) = -1092;
    a(40) = -(56*(39*nu + 43))/(nu + 1);
    b(40) =  (56*(41*nu + 37))/(nu + 1);
    a(41) =  2870;
    b(41) =  2870;
    a(42) =  (56*(41*nu + 45))/(nu + 1);
    b(42) = -(56*(43*nu + 39))/(nu + 1);
    a(43) = -1204;
    b(43) = -1204;
    a(44) = -(8*(43*nu + 47))/(nu + 1);
    b(44) =  (8*(45*nu + 41))/(nu + 1);

    if Ident == 1
        % Displacement in x-direction
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3 + a(11)*x.^4 + ...
                       a(12)*x.^3.*y + a(13)*x.^2.*y.^2 + ...
                       a(14)*x.*y.^3 + a(15)*y.^4 + a(16)*x.^5 + ...
                       a(17)*x.^4.*y + a(18)*x.^3.*y.^2 + ...
                       a(19)*x.^2.*y.^3 + a(20)*x.*y.^4 + a(21)*y.^5 + ...
                       a(22)*x.^6 + a(23)*x.^5.*y + a(24)*x.^4.*y.^2 + ...
                       a(25)*x.^3.*y.^3 + a(26)*x.^2.*y.^4 + ...
                       a(27)*x.*y.^5 + a(28)*y.^6 + a(29)*x.^7 + ...
                       a(30)*x.^6.*y + a(31)*x.^5.*y.^2 + ...
                       a(32)*x.^4.*y.^3 + a(33)*x.^3.*y.^4 + ...
                       a(34)*x.^2.*y.^5 + a(35)*x.*y.^6 + a(36)*y.^7 + ...
                       a(37)*x.^8 + a(38)*x.^7.*y + a(39)*x.^6.*y.^2 + ...
                       a(40)*x.^5.*y.^3 + a(41)*x.^4.*y.^4 + ...
                       a(42)*x.^3.*y.^5 + a(43)*x.^2.*y.^6 + ...
                       a(44)*x.*y.^7 + a(45)*y.^8);
    elseif Ident == 2
        % Displacement in y-direction
        DispVal = Amp*(b(1) + b(2)*x + b(3)*y + b(4)*x.^2 + b(5)*x.*y + ...
                       b(6)*y.^2 + b(7)*x.^3 + b(8)*x.^2.*y + ...
                       b(9)*x.*y.^2 + b(10)*y.^3 + b(11)*x.^4 + ...
                       b(12)*x.^3.*y + b(13)*x.^2.*y.^2 + ...
                       b(14)*x.*y.^3 + b(15)*y.^4 + b(16)*x.^5 + ...
                       b(17)*x.^4.*y + b(18)*x.^3.*y.^2 + ...
                       b(19)*x.^2.*y.^3 + b(20)*x.*y.^4 + b(21)*y.^5 + ...
                       b(22)*x.^6 + b(23)*x.^5.*y + b(24)*x.^4.*y.^2 + ...
                       b(25)*x.^3.*y.^3 + b(26)*x.^2.*y.^4 + ...
                       b(27)*x.*y.^5 + b(28)*y.^6 + b(29)*x.^7 + ...
                       b(30)*x.^6.*y + b(31)*x.^5.*y.^2 + ...
                       b(32)*x.^4.*y.^3 + b(33)*x.^3.*y.^4 + ...
                       b(34)*x.^2.*y.^5 + b(35)*x.*y.^6 + b(36)*y.^7 + ...
                       b(37)*x.^8 + b(38)*x.^7.*y + b(39)*x.^6.*y.^2 + ...
                       b(40)*x.^5.*y.^3 + b(41)*x.^4.*y.^4 + ...
                       b(42)*x.^3.*y.^5 + b(43)*x.^2.*y.^6 + ...
                       b(44)*x.*y.^7 + b(45)*y.^8);
    else
        % Displacement in both directions
        DispVal = Amp*(a(1) + a(2)*x + a(3)*y + a(4)*x.^2 + a(5)*x.*y + ...
                       a(6)*y.^2 + a(7)*x.^3 + a(8)*x.^2.*y + ...
                       a(9)*x.*y.^2 + a(10)*y.^3 + a(11)*x.^4 + ...
                       a(12)*x.^3.*y + a(13)*x.^2.*y.^2 + ...
                       a(14)*x.*y.^3 + a(15)*y.^4 + a(16)*x.^5 + ...
                       a(17)*x.^4.*y + a(18)*x.^3.*y.^2 + ...
                       a(19)*x.^2.*y.^3 + a(20)*x.*y.^4 + a(21)*y.^5 + ...
                       a(22)*x.^6 + a(23)*x.^5.*y + a(24)*x.^4.*y.^2 + ...
                       a(25)*x.^3.*y.^3 + a(26)*x.^2.*y.^4 + ...
                       a(27)*x.*y.^5 + a(28)*y.^6 + a(29)*x.^7 + ...
                       a(30)*x.^6.*y + a(31)*x.^5.*y.^2 + ...
                       a(32)*x.^4.*y.^3 + a(33)*x.^3.*y.^4 + ...
                       a(34)*x.^2.*y.^5 + a(35)*x.*y.^6 + a(36)*y.^7 + ...
                       a(37)*x.^8 + a(38)*x.^7.*y + a(39)*x.^6.*y.^2 + ...
                       a(40)*x.^5.*y.^3 + a(41)*x.^4.*y.^4 + ...
                       a(42)*x.^3.*y.^5 + a(43)*x.^2.*y.^6 + ...
                       a(44)*x.*y.^7 + a(45)*y.^8);
    end
end
%-------------------------------------------------------------------------%

end