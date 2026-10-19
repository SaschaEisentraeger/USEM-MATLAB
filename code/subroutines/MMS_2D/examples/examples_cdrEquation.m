%% examples: convection-diffusion-reaction equation
% similar solutions as used for the Reynolds equation, but scaled to the
% domain (0,1;0,1)

%% choose example
exampleNo = 1;                                                              % 1: smooth, 2: large gradients
parameterChoice = 0;                                                        % see choice below

%% problem-dependent functions
switch exampleNo
    case 1                                                                  % simple smooth solution
        u = @(x,y) 1/3*(1-cos(4*pi*x)).*sin(2*pi*x).*0.5.*(sin(pi*y));      % unknown field

    case 2                                                                  % solution involving large gradients
        u= @(x,y) 0.25*((1-exp(100*x))/(1-exp(100))-...
            1+(cos(x*pi)+1)/2).*(sin(pi*y));                                % unknown field
end

switch parameterChoice
    case 0                                                                  % same diffusivity, convection and reaction

        k = 1;                                                              % diffusivity
        a = [1, 0];                                                         % convection
        s = 1;                                                              % reaction

    case 1                                                                  % highly diffusive problem

        k = 1;                                                              % diffusivity
        a = [1e-4, 0];                                                      % convection
        s = 1e-4;                                                           % reaction

    case 2                                                                  % convection-dominated problem

        k = 1e-4;                                                           % diffusivity
        a = [1, 0];                                                         % convection
        s = 1e-4;                                                           % reaction

    case 3                                                                  % reaction-dominated problem 

        k = 1;                                                              % diffusivity
        a = [1e-4, 0];                                                      % convection
        s = 1;                                                              % reaction

end

%% evaluate PDE

fFun = @(u) -k*lapl(u) + a*grad(u) + s*u;                                   % define PDE
[f,fSym] = manufacturedSolution({u},fFun);                                  % compute body load

%% plotting
domain = [[0,1];[0,1]];                                                     % x- and y-limits for plotting
plotManufactured(u,f,domain,0,100);                                         % call plot function