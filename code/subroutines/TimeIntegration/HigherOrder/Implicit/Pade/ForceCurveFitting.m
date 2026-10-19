%% Force vector --> polynomial curve fitting
syms A x

% Exponential function
func = exp(x);
% Order of the numerator
L = 4;
% Order of the denominator
M = L;
% Order of the polynomial approximating the force
N = L;
% Pade series
[CoeffN,CoeffD,R] = PadeSeriesCoeff(func,L,M,false);

% Denominator and numerator polynomials for the Pade-series of the matrix
% exponential
P = A.^(0:L).*CoeffN;
P = sum(sign(P(1))*P);
Q = A.^(0:M).*CoeffD;
Q = sum(sign(Q(1))*Q);

% Coefficients
C = sym('C', [1 N+1]);
for i = 0:N
    if i == 0
        C(i+1) = A\(P-Q);
    else
        C(i+1) = A\(i*C(i) + (-1/2)^(i)*(P-(-1)^(i)*Q));
    end
end

% Result
disp(expand(C))

%% Approximated force terms
% Lagrange interpolation based on GLL-points
syms s f0 f1 f2 f3
xi = 2*s-1;

% Linear
% Quadratic approximation
N1 = 1/2*(1-xi);
N2 = 1/2*(1+xi);
% simplify(expand(N1))
% simplify(expand(N2))
dN1 = diff(N1);
dN2 = diff(N2);

Fmn1 = subs(N1,0.5)*f0 + subs(N2,0.5)*f1;
Fmn1 = 1/factorial(0)*Fmn1;
dFmn1 = subs(dN1,0.5)*f0 + subs(dN2,0.5)*f1;
dFmn1 = 1/factorial(1)*dFmn1;
fmn1_L = [Fmn1;dFmn1];

% Quadratic approximation
N1 = 1/2*xi*(xi-1);
N2 = 1 - xi^2;
N3 = 1/2*xi*(xi+1);
% simplify(expand(N1))
% simplify(expand(N2))
% simplify(expand(N3))
dN1 = diff(N1);
dN2 = diff(N2);
dN3 = diff(N3);
ddN1 = diff(dN1);
ddN2 = diff(dN2);
ddN3 = diff(dN3);

Fmn2 = subs(N1,0.5)*f0 + subs(N2,0.5)*f1 + subs(N3,0.5)*f2;
Fmn2 = 1/factorial(0)*Fmn2;
dFmn2 = subs(dN1,0.5)*f0 + subs(dN2,0.5)*f1 + subs(dN3,0.5)*f2;
dFmn2 = 1/factorial(1)*dFmn2;
ddFmn2 = subs(ddN1,0.5)*f0 + subs(ddN2,0.5)*f1 + subs(ddN3,0.5)*f2;
ddFmn2 = 1/factorial(2)*ddFmn2;
fmn2_L = [Fmn2;dFmn2;ddFmn2];

% Cubic approximation
N1 = -1/8*(5*xi^2-1)*(xi-1);
N2 = 5/8*(xi^2-1)*(sqrt(5)*xi-1);
N3 = -5/8*(xi^2-1)*(sqrt(5)*xi+1);
N4 = 1/8*(5*xi^2-1)*(xi+1);
% simplify(expand(N1))
% simplify(expand(N2))
% simplify(expand(N3))
% simplify(expand(N4))
dN1 = diff(N1);
dN2 = diff(N2);
dN3 = diff(N3);
dN4 = diff(N4);
ddN1 = diff(dN1);
ddN2 = diff(dN2);
ddN3 = diff(dN3);
ddN4 = diff(dN4);
dddN1 = diff(ddN1);
dddN2 = diff(ddN2);
dddN3 = diff(ddN3);
dddN4 = diff(ddN4);

Fmn3 = subs(N1,0.5)*f0 + subs(N2,0.5)*f1 + subs(N3,0.5)*f2 + subs(N4,0.5)*f3;
Fmn3 = 1/factorial(0)*Fmn3;
dFmn3 = subs(dN1,0.5)*f0 + subs(dN2,0.5)*f1 + subs(dN3,0.5)*f2 + subs(dN4,0.5)*f3;
dFmn3 = 1/factorial(1)*dFmn3;
ddFmn3 = subs(ddN1,0.5)*f0 + subs(ddN2,0.5)*f1 + subs(ddN3,0.5)*f2 + subs(ddN4,0.5)*f3;
ddFmn3 = 1/factorial(2)*ddFmn3;
dddFmn3 = subs(dddN1,0.5)*f0 + subs(dddN2,0.5)*f1 + subs(dddN3,0.5)*f2 + subs(dddN4,0.5)*f3;
dddFmn3 = 1/factorial(3)*dddFmn3;
fmn3_L = [Fmn3;dFmn3;ddFmn3;dddFmn3];

%% Implementation
syms f0 f1 f2 f3
% Linear
xi = lglnodes(1);
xi = flip(xi);
% Mapping to the interval [0,1]
s  = 1/2*(xi + 1);
a1 = (s-0.5).^[ 0 1];
fmn1 = (a1'*a1)\(a1')*[f0; f1];

% Quadratic
% GLL-points in the interval [-1,1]
xi = lglnodes(2);
xi = flip(xi);
% Mapping to the interval [0,1]
s  = 1/2*(xi + 1);
a2 = (s-0.5).^[ 0 1 2];
fmn2 = (a2'*a2)\(a2')*[f0; f1; f2];

% Cubic
% GLL-points in the interval [-1,1]
xi = lglnodes(3);
xi = flip(xi);
% Mapping to the interval [0,1]
s  = 1/2*(xi + 1);
a3 = (s-0.5).^[ 0 1 2 3];
fmn3 = (a3'*a3)\(a3')*[f0; f1; f2; f3];