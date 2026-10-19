%% Generic function for the use in manufactured solutions
% Malaya, N. / Estacio-Hiroms, K. C. / Stogner, R. H. / Schulz, K. W. /
% Bauman, P. T. / Carey, G. F.:
% MASA: A Library for Verification Using Manufactured and Analytical
% Solutions, Engineering with Computers, Vol. 29, No. 4, p. 487-496, 2012 
% Springer Science and Business Media LLC 
%-------------------------------------------------------------------------%
% u0, ux, uy, ut, aux, auy, aut, L are real constants
% Instead of sine-functions, it is also possible to use cosine-functions or
% polynomials.
clearvars; clc; close all; format long

% Sum = 1; Product = 2;
method = 1;

% Parameters
u0 = 1;
ux = 1;
uy = 1;
ut = 1;
aux = 1;
auy = 1;
aut = 1;
L = 1;

x = 0:0.01:L;
y = 0:0.01:L;
t = 0:0.01:L;

if method == 1
	% Simple sum of different terms
	u = @(x,y,t) u0 + ux*sin(aux*pi*x/L) + uy*sin(auy*pi*y/L) + ut*sin(aut*pi*t/L);
elseif method == 2
	% Product terms
	u = @(x,y,t) u0 + ux*uy*ut*sin(aux*pi*x/L)*sin(auy*pi*y/L)*sin(aut*pi*t/L);
else
    u = [];
end