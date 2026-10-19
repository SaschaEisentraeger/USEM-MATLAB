%% example: Poisson equation

%% problem-dependent functions
u = @(x,y) sin(pi*x)*sin(pi*y);                                             % unknown field

%% evaluate PDE
fFun = @(u) lapl(u);                                                        % define PDE (left-hand side)
f = manufacturedSolution({u},fFun);                                         % compute body load

%% plotting
domain = [[0,1];[0,1]];                                                     % x- and y-limits for plotting
plotManufactured(u,f,domain,0,100);                                         % call plot function