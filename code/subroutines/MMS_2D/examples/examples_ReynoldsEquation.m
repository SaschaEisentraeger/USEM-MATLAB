%% examples: Reynolds equation with cavitation
% Stationary manufactured solutions in 
% "Stabilized finite elements for the solution of the Reynolds equation 
% considering cavitation", Gravenkamp et al.
% see Sections 8.1, 8.2

%% choose example
exampleNo = 1;                                                              % choose 1 or 2 for examples in sections 8.1 or 8.2 

%% problem-dependent functions
switch exampleNo
    case 1                                                                  % section 8.1
        u = @(x,y) 1/3*(1-cos(2*x)).*sin(x).*0.5.*(1+cos(pi*y));            % unknown field
        g = @(x,y) 1/pi*atan(1/3*(1-cos(2*x))*sin(x)*0.5*(1+cos(pi*y))...
            /(1-0.98))+0.5;                                                 % switch function

    case 2
        u= @(x,y) 0.25*((1-exp(100*x/(2*pi)))/(1-exp(100))-...
            1+(cos(x/2)+1)/2).*(1+cos(pi*y));                               % unknown field
        g = @(x,y) 1/pi*atan(( 0.25*((1-exp(100*x/(2*pi)))/(1-exp(100))...
            -1+(cos(x/2)+1)/2).*(1+cos(pi*y)))/(1-0.98)) + 1/2;             % switch function
end

%% evaluate PDE
H = @(x,y) 1-0.5*cos(x-pi);                                                 % gap function
mu=1;                                                                       % viscosity
fFun = @(u,g,H) -1/12*div(H.^3/mu*grad(g*u))+derX(H,1)-derX((g-1)*H*u,1);   % define PDE
[f,fSym] = manufacturedSolution({u,g,H},fFun);                              % compute body load

%% plotting
domain = [[0,2*pi];[-1,1]];                                                 % x- and y-limits for plotting
plotManufactured(u,f,domain,0,100);                                         % call plot function