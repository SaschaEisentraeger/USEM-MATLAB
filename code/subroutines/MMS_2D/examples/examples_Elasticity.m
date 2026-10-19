%% example: Structural mechanics -- Statics
clearvars; clc; close all; dbstop if error; format long

version = 2;

%% Material properties
E = 70e9;
nu = 0.3;
rho = 2700;
antype = 1;                                                                 % (1): plane stress / (2) plane strain
flag_plot = 0;

if antype == 1
    % Plane stress -- Lame parameters
    lambda = E*nu/(1-nu^2);
    mu = E/(2*(1+nu));
elseif antype == 2
    % Plane strain -- Lame parameters
    lambda = E*nu/((1+nu)*(1-2*nu));
    mu = E/(2*(1+nu));
end

%% problem-dependent functions
%----------------------Define displacement functions----------------------%

syms x y
% MMS example
ampx = 1;%1e-5;
ampy = 1;%1e-5;
ux_str = 'sin(5*pi*x).*cos(8*pi*y)*tanh(2*pi*x*y)';
uy_str = 'cos(7*pi*x).*sin(3*pi*y)*tanh(2*pi*x*y)';

% % Lee and Bathe (1993): Linear bending problem
% ampx = 1;
% ampy = 1;
% L = 10;
% c = 2;
% ux_str = '(-92/(c*L)*(y+1).^3 + 120/(c*L)*x.^2.*(y+1) -  60/L*x.^2 - 240/c*x.*(y+1) + 138/L*(y+1).^2 + 120*x    - 46*c/L*(y+1))/E';
% uy_str = '(-40/(c*L)*x.^3 -  36/(c*L)*x.*(y+1).^2 + 120/c*x.^2 +  36/L*x.*(y+1) +  36/c*(y+1).^2 + 46*c/L*x - 36*(y+1))/E';
% ux_str = replace(ux_str, 'c', num2str(c));
% ux_str = replace(ux_str, 'L', num2str(L));
% ux_str = replace(ux_str, 'E', num2str(E));
% uy_str = replace(uy_str, 'c', num2str(c));
% uy_str = replace(uy_str, 'L', num2str(L));
% uy_str = replace(uy_str, 'E', num2str(E));

% % Test for cubic elements
% ampx = 1;
% ampy = 1;
% ux_str = '1 + x + y + x^2 + x*y + y^2 + x^3 + x^2*y + x*y^2 + y^3';
% uy_str = '1 + x + y + x^2 + x*y + y^2 + x^3 + x^2*y + x*y^2 + y^3';
%------------------------------Add amplitude------------------------------%
ux_sym = ampx*str2sym(ux_str);
uy_sym = ampy*str2sym(uy_str);
%-------------------------------------------------------------------------%
% Create function handle
ux = matlabFunction(ux_sym,'Vars',[x,y]);
uy = matlabFunction(uy_sym,'Vars',[x,y]);
% Write m-file for Dirichlet boundary conditions 
u  = FunctionDisplacement(ux,uy,true);

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
    fFun = @(ux,uy) -div(sigma);                                            % define PDE (left-hand side)
    % Write m-file for body load (source term)
    f = manufacturedSolution({ux,uy},fFun,true);                            % compute body load
else
    fFun = @(u) -div(sigma(u));                                             % define PDE (left-hand side)
    % Write m-file for body load (source term)
    f = manufacturedSolution({u},fFun,true);                                % compute body load
end

%% plotting
if flag_plot == 1
    domain = [[0,1];[0,1]];                                                % x- and y-limits for plotting
    varNames = {'u_x','u_y'};                                              % names of each dof for plots
    plotManufactured(u,f,domain,0,100,varNames);                           % call plot function
end

%% Compute strain energy over the domain [0,1] x [0,1]
% Elasticity matrix (3D)
C3D = HookeMatrix3d(E,nu);

% Elasticity matrix (2D)
if antype == 1
    % Plane stress -- Lame parameters
    [C] = PlaneStressConstitutiveMatrices(C3D);
elseif antype == 2
    % Plane strain -- Lame parameters
    [C] = PlaneStrainConstitutiveMatrices(C3D);
end

% Convert function handle into a symbolic expression
epsilon_sym = sym(epsilon(u));
epsilon_sym_vec = [epsilon_sym(1,1); epsilon_sym(2,2); 2*epsilon_sym(1,2)];

% % Compute strain energy symbolically
% inner_int = int(transpose(epsilon_sym_vec)*C*epsilon_sym_vec, x, 0, 1);
% outer_int = 0.5 * int(inner_int, y, 0, 1);
% Energy = double(outer_int);

% Compute strain energy numerically
integrand = transpose(epsilon_sym_vec)*C*epsilon_sym_vec;
f = matlabFunction(integrand, 'Vars', [x y]);
Energy = 0.5 * integral2(f, 0, 1, 0, 1, 'AbsTol', 1e-15, 'RelTol', 1e-15);

% % Compute strain energy numerically
% integrand = transpose(epsilon_sym_vec)*C*epsilon_sym_vec;
% f = matlabFunction(integrand, 'Vars', [x y]);
% Energy = 0.5 * integral2(f, 0, L, -c/2, c/2, 'AbsTol', 1e-15, 'RelTol', 1e-15);

%% --------------------------------------------------------------------- %%
% syms x y ampx ampy mu lambda
% % MMS example
% ux = ampx*sin(5*pi*x).*cos(8*pi*y)*tanh(2*pi*x*y);
% uy = ampy*cos(7*pi*x).*sin(3*pi*y)*tanh(2*pi*x*y);
% 
% % Displacement derivatives
% ux_dx = derX(ux,1);
% uy_dy = derY(uy,1);
% ux_dy = derY(ux,1);
% uy_dx = derX(uy,1);
% 
% % Strain
% eps_xx = ux_dx;
% eps_yy = uy_dy;
% gamma_xy = 1/2*(ux_dy + uy_dx);
% 
% % Stress
% sigma_xx = 2*mu *eps_xx + lambda*(eps_xx+eps_yy);
% sigma_yy = 2*mu *eps_yy + lambda*(eps_xx+eps_yy);
% tau_xy   = 2*mu *gamma_xy;
% 
% % Source term
% f_x = diff(sigma_xx,x) + diff(tau_xy,y);
% f_y = diff(sigma_yy,y) + diff(tau_xy,x);