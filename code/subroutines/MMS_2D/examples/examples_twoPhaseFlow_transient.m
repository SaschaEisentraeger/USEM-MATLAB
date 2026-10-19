%% example: transient manufactured solution of 2-phase flow
% This is the transient version of the first example in 
% "A stabilized finite element method for modeling dispersed multiphase 
% flows using orthogonal subgrid scales", Gravenkamp et al.
% see Section 6.1.1

%% problem-dependent functions

mu1 = 0.1;                                                                  % viscosities
mu2 = 0.2;
rho1 = 1;                                                                   % densities
rho2 = 1;
dragCoef = 0.5;                                                             % drag coefficient
A1 = 0.3;                                                                   % volume fractions
A2 = 0.7;

f  = @(x) x.^2.*(1-x).^2;                                                   % function f(x) used in spatial dependency 
df = @(x) 2*x.*(1-x).^2 - 2*x.^2.*(1-x);                                    % df/dx
ddf = @(x) 2*(1-x).^2 - 8*x.*(1-x) + 2*x.^2;                                % d²f/dx²
dddf = @(x) 24*x-12;                                                        % d³f/dx³
g  = @(t) cos(pi*t).*exp(-t);                                               % time-dependent part of solution (denoted f_t in paper)
dg = @(t) -pi*sin(pi*t).*exp(-t) - cos(pi*t).*exp(-t);                      % derivative wrt t

UgradU = @(x,y) 100^2*[...                                                  % u*grad(u), phase 1
    f(x).*df(y).*df(x).*df(y)-df(x).*f(y).*f(x).*ddf(y);...
    -f(x).*df(y).*ddf(x).*f(y)+df(x).*f(y).*df(x).*df(y)];
UgradU2 = @(x,y) 4*UgradU(x,y);                                             % phase 2

p  = @(x,y)  100* x.^2;                                                     % pressure
gradP = @(x,y) 100*[2*x;0];                                                 % pressure gradient
ux = @(x,y)  100*x.^2.*(1-x).^2.*(2*y.*(1-y).^2 + y.^2.*2.*(1-y)*(-1));     % velocity, x-component, phase 1
uy = @(x,y) -100*y.^2.*(1-y).^2.*(2*x.*(1-x).^2 + x.^2.*2.*(1-x)*(-1));     % velocity, y-component, phase 1

u1 = @(x,y) [ ...                                                           % velocity, phase 1
    +100*x.^2.*(1-x).^2.*(2*y.*(1-y).^2 + y.^2.*2.*(1-y)*(-1)); ...
    -100*y.^2.*(1-y).^2.*(2*x.*(1-x).^2 + x.^2.*2.*(1-x)*(-1))];
u2 = @(x,y) 2*[ ...                                                         % velocity, phase 2
    +100*x.^2.*(1-x).^2.*(2*y.*(1-y).^2 + y.^2.*2.*(1-y)*(-1)); ...
    -100*y.^2.*(1-y).^2.*(2*x.*(1-x).^2 + x.^2.*2.*(1-x)*(-1))] ;

U=@(x,y,t) [                                                                % collect all dofs
    ux(x,y)*g(t); uy(x,y)*g(t);...
    2*ux(x,y)*g(t); 2*uy(x,y)*g(t); ...
    A1; A2; p(x,y) ];


%% evaluate PDE
% terms need to be treated separately because of time-dependent function

epsi = @(u) 1/2*(grad(u)+(grad(u)).');                                      % strain rate
f1Fun = @(u1,u2) -2*mu1*div(A1*epsi(u1)) -dragCoef*(u2-u1);                 % part of momentum equation multiplied by g(t)
f1 = manufacturedSolution({u1,u2},f1Fun);                                   % corresponding rhs
f2Fun = @(u1,u2) -2*mu2*div(A2*epsi(u2)) -dragCoef*(u1-u2);                 % same for phase 2
f2 = manufacturedSolution({u1,u2},f2Fun);

f1 = @(x,y) [f1(x,y)/A1;f2(x,y)/A2;0;0;0];                                  % collect loads to be multiplied by g(t)
f2 = @(x,y) [rho1*u1(x,y);rho2*u2(x,y);0;0;0];                              % term to be multiplied by dg(t)/dt
f3 = @(x,y) [gradP(x,y);gradP(x,y);0;0;0];                                  % time-independent term
f4 = @(x,y) [rho1*UgradU(x,y);rho2*UgradU2(x,y);0;0;0];                     % terms to be multiplied by g(t)^2

ft1 = g;                                                                    % time function of terms g(t)
ft2 = dg;                                                                   % time function for mass term dg(t)/dt
ft3 = @(t) 1;                                                               % constant in time
ft4 = @(t) g(t).^2;                                                         % g(t)^2

f = @(x,y,t) f1(x,y).*ft1(t) + f2(x,y).*ft2(t)...
    + f3(x,y).*ft3(t)+f4(x,y).*ft4(t);                                      % assemble time-dependent load

%% plotting
domain = [[0,1];[0,1]];                                                     % x- and y-limits for plotting
varNames = {'u_1^x','u_1^y','u_2^x','u_2^y','\alpha_1','\alpha_2','p'};     % names of each dof for plots
plotManufactured(U,f,domain,1,100,varNames);                                % call plot function