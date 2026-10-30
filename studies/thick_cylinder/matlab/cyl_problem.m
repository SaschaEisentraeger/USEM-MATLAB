function prob = cyl_problem(kind, E, nu)
%CYL_PROBLEM Exact solution, strains, stresses and body force.
%
%   prob = cyl_problem(kind, E, nu) for plane stress with Young's modulus E
%   and Poisson's ratio nu. The fields are derived symbolically:
%       epsilon = sym(grad u),  sigma = 2 mu epsilon + lambda tr(epsilon) I,
%       f = -div(sigma).
%
%   kind = 'lame'   : thick cylinder under internal pressure (Lame solution,
%                     a = 1, b = 3, p = 1000), Rajendran (2010), Eq. (25)
%          'poly'   : manufactured polynomial field of total degree 5
%          'smooth' : manufactured smooth non-polynomial (analytic) field
%
%   prob.ux, uy, exx, eyy, gxy, sxx, syy, sxy, fx, fy : function handles
%   @(x,y) acting elementwise (gxy is the engineering shear strain).

syms x y
X = x/3;  Y = y/3;          % non-dimensional coordinates on [0,1]^2
switch kind
    case 'lame'
        a = sym(1);  b = sym(3);  p = sym(1000);
        r2 = x^2 + y^2;
        % u_r/r of Eq. (25) of Rajendran (2010)
        ur_r = a^2*p/(sym(E)*(b^2 - a^2))*((1 - sym(nu)) + b^2/r2*(1 + sym(nu)));
        u = [ur_r*x, ur_r*y];
    case 'poly'
        u = [ (1 + X + 2*Y + X*Y - 3*X^2*Y + X^3 + 2*X^4*Y - X^2*Y^3 + X^5 - 2*Y^5)/10, ...
             (1 - X + Y^2 + 2*X^2*Y - X^3 + X*Y^3 - 3*X^3*Y^2 + Y^5 + 2*X*Y^4 - X^5)/10 ];
    case 'smooth'
        u = [ (sin(3*X)*cos(2*Y) + exp(X*Y))/10, ...
              (cos(2*X)*sin(3*Y) + Y*exp(-X))/10 ];
    otherwise
        error('cyl_problem:kind', 'Unknown problem "%s".', kind);
end

% Plane stress Lame parameters
Es = sym(E);  nus = sym(nu);
lambda = Es*nus/(1 - nus^2);
mu = Es/(2*(1 + nus));

grad = @(u) [diff(u(1), x), diff(u(2), x); diff(u(1), y), diff(u(2), y)];
epsilon = @(u) 1/2*(grad(u) + (grad(u)).');
sigma = @(u) 2*mu*epsilon(u) + lambda*trace(epsilon(u))*eye(2);
divergence = @(s) [diff(s(1,1), x) + diff(s(2,1), y); diff(s(1,2), x) + diff(s(2,2), y)];
fFun = @(u) -divergence(sigma(u));

epsS = epsilon(u);
sigS = sigma(u);
fS = fFun(u);

prob.kind = kind;
prob.ux  = toHandle(u(1), x, y);
prob.uy  = toHandle(u(2), x, y);
prob.exx = toHandle(epsS(1,1), x, y);
prob.eyy = toHandle(epsS(2,2), x, y);
prob.gxy = toHandle(2*epsS(1,2), x, y);
prob.sxx = toHandle(sigS(1,1), x, y);
prob.syy = toHandle(sigS(2,2), x, y);
prob.sxy = toHandle(sigS(1,2), x, y);
prob.fx  = toHandle(fS(1), x, y);
prob.fy  = toHandle(fS(2), x, y);
prob.symbolic.u = u;
prob.symbolic.f = fS;

end

function h = toHandle(expr, x, y)
% Elementwise function handle; constants are expanded to the input size.
g = matlabFunction(simplify(expr), 'Vars', [x, y]);
h = @(xx, yy) g(xx, yy) + zeros(size(xx));
end
