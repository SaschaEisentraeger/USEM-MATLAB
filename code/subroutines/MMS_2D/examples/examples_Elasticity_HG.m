
%% example: Structural mechanics -- Statics
%% Material properties
E = 70e9;
nu = 0.3;
antype = 1;
% Lame parameters
lambda = E*nu/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));

%% Constitutive matrix
[C] = HookeMatrix3d(E,nu);
if antype == 1
    C2d = PlaneStressConstitutiveMatrices(C);
else
    C2d = PlaneStrainConstitutiveMatrices(C);
end
%% problem-dependent functions
u = @(x,y) [sin(pi*x)*sin(pi*y); sin(pi*x)*sin(pi*y)];

%% Strain definition
epsilon = @(u) 1/2*(grad(u)+(grad(u)).');                                   % strain

%% Stress definition
sigma = @(u) 2*mu*epsilon(u) + lambda*trace(epsilon(u))*eye(2);

%% evaluate PDE
fFun = @(u) -div(sigma(u));                                                 % define PDE (left-hand side)
f = manufacturedSolution({u},fFun);                                         % compute body load

%% plotting
domain = [[0,1];[0,1]];                                                     % x- and y-limits for plotting
varNames = {'u_x','u_y'};                                                   % names of each dof for plots
U = @(x,y) [ux(x,y); uy(x,y)];
plotManufactured(U,f,domain,0,100,varNames);                                % call plot function


