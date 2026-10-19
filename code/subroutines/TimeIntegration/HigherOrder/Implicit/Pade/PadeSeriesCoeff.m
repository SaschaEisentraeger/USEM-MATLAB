%% Pade series approximant
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by Sascha Eisentraeger (12/2020)
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%
% Calculate the coefficients of the numerator and denominator polynomials
% for a Pade series.
%
%-------------------------------------------------------------------------%
%% Input
%
% func:         Symbolic function that is approximated by a Pade series
% L:            Order of the numerator polynomial
% M:            Order of the denominator polynomial
% flag_norm:    If true the maximum coefficient in the expansion is scaled
%               to unity
%
%-------------------------------------------------------------------------%
%% Output
%
% CoeffN:   Coefficients of the numerator polynomial (x^L,x^(L-1),...,x^0)
% CoeffD:   Coefficients of the denominator polynomial (x^L,x^(L-1),...,x^0)
% R:        Roots of the denominator polynomial
%
%-------------------------------------------------------------------------%
% % Example
% L = 2;
% M = 2;
% % Initialize symbolic variable x
% syms x
% % Exponential function
% func = exp(x);
% % Pade series
% [CoeffN,CoeffD,R] = PadeSeriesCoeff(func,L,M,false);
%-------------------------------------------------------------------------%
function [CoeffN,CoeffD,R] = PadeSeriesCoeff(func,L,M,flag_norm)

% Check the number of input arguments
if nargin < 3
	error('The minimum is three input parameters!')
elseif nargin == 3
	flag_norm = false;
end

% Pade series of order L/M
PS = pade(func,'Order',[L,M]);

% Extract Numerator and Denominator
[N,D] = numden(PS);

% Numerator: Coefficients (sorted form x^0 to x^L)
C1 = coeffs(N,'All');
% Denominator: Coefficients (sorted form x^0 to x^M)
C2 = coeffs(D,'All');

% Numerical values --> sort from x^n to x^0 (use: polyval, polyder, polyint
% to calculate the value, derivative, or integral of a polynomial function)
CoeffN = flip(eval(C1));
CoeffD = flip(eval(C2));

if flag_norm
	% Find identical coefficients in the numerator and denominator
	CoeffDN = intersect(CoeffN,CoeffD);
	% Largest identical coefficient
	maxCoeff = max(abs(CoeffDN));

	% Scale the coefficients: max(Coeff) = 1
	if ~isempty(CoeffDN)
		CoeffN = CoeffN/maxCoeff;
		CoeffD = CoeffD/maxCoeff;
	end
end

% Compute the roots of the denominator polynomial
R = roots(flip(CoeffD));

end