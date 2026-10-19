%% example: Structural mechanics -- Statics
clearvars; clc; close all; dbstop if error; format long

version = 2;

%% Material properties
E = 1;%70e9;
nu = 0;%0.3;
rho = 1;%2700;
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
ux = @(x,y) sin(pi*x).*sin(pi*y);                                           % unknown field (amplitude)
uy = @(x,y) sin(pi*x).*sin(pi*y);
u  = FunctionDisplacement(ux,uy,true);

g  = @(t) cos(pi*t);                                                        % time signal
[dg,ddg]  = FunctionTimeSignal(g,true);

uxt = @(x,y,t) ux(x,y).*g(t);                                               % unknown field (time-dependent)
uyt = @(x,y,t) uy(x,y).*g(t);
ut  = FunctionDisplacement(uxt,uyt,true);

%% Strain definition
if version == 1
    gradUx = grad(ux);
    gradUy = grad(uy);
    epsilon = 1/2*([gradUx, gradUy] + [gradUx, gradUy]');
else
    epsilon = @(u) 1/2*(grad(u)+(grad(u)).');
end

%% Stress definition
if version == 1
    sigma = 2*mu*epsilon + lambda*trace(epsilon)*eye(2);
else
    sigma = @(u) 2*mu*epsilon(u) + lambda*trace(epsilon(u))*eye(2);
end

%% evaluate PDE
if version == 1
    f1Fun = @(ux,uy) -div(sigma);                                           % define PDE (left-hand side)
    f1 = manufacturedSolution({ux,uy},f1Fun,true,'BodyLoadFuncDiv.m');
else
    f1Fun = @(u) -div(sigma(u));                                            % define PDE (left-hand side)
    f1 = manufacturedSolution({u},f1Fun,true,'BodyLoadFuncDiv.m');
end
f2Fun = @(u) rho*u';
f2 = manufacturedSolution({u},f2Fun,true,'BodyLoadFuncRho.m');
f = @(x,y,t) f1(x,y).*g(t) + f2(x,y)*ddg(t);                                % assemble time-dependent load
f = manufacturedSolutionTransient(f,true);

%% plotting
domain = [[0,1];[0,1]];                                                     % x- and y-limits for plotting
varNames = {'u_x','u_y'};                                                   % names of each dof for plots
plotManufactured(ut,f,domain,0.5,100,varNames);                             % call plot function