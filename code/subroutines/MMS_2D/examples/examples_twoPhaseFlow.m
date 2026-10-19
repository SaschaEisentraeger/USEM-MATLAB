%% example: stationary manufactured solution of 2-phase flow
% Stationary manufactured solutions in 
% "A stabilized finite element method for modeling dispersed multiphase 
% flows using orthogonal subgrid scales", Gravenkamp et al.
% see Sections 6.1.1, 6.1.2

%% choose example
exampleNo = 1;                                                              % choose 1 or 2 for examples in sections 6.1.1 or 6.1.2 

%% problem-dependent functions
switch exampleNo
    case 1                                                                  % section 6.1.1

        MPF_dragModel = 'constant';                                         % constant drag term 
        dragCoef=0.5;                                                       % drag coefficient
        rho1 = 1;                                                           % densities
        rho2 = 1;
        mu1 = 0.1;                                                          % viscosities
        mu2 = 0.2;

        u1=@(x,y) 100*[...                                                  % velocities
            +x.^2.*(1-x).^2 .* (2*y.*(1-y).^2 - 2*y.^2.*(1-y)); ...
            -y.^2.*(1-y).^2 .* (2*x.*(1-x).^2 - 2*x.^2.*(1-x))  ]+1;
        u2=@(x,y) 100*2*[...
            +x.^2.*(1-x).^2 .* (2*y.*(1-y).^2 - 2*y.^2.*(1-y)); ...
            -y.^2.*(1-y).^2 .* (2*x.*(1-x).^2 - 2*x.^2.*(1-x))  ]+2;

        alpha1 = @(x,y) ones(size(x))*0.3;                                  % volume fractions
        alpha2 = @(x,y) ones(size(x))*0.7;
        p = @(x,y) 100*x.^2;                                                % pressure

    case 2                                                                  % section 6.1.2

        MPF_dragModel = 'Hiltunen';                                         % drag model taken from Hiltunen's paper
        Cd = 0.44;                                                          % coefficient in drag model
        d = 3/4*Cd;                                                         % particle diameter (such that g12 = alpha1 |u2-u1|)
        rho1 = 2;                                                           % densities
        rho2 = 1/3;
        mu1 = 1;                                                            % viscosities
        mu2 = 1;

        u1 = @(x,y)[x+1;   -(1.2*x+0.6+0.2)./(0.6*x+0.2).*y+4.5];           % velocities
        u2 = @(x,y)[x.^2+1; (1.8*x.^2+0.6-1.6*x)./(0.8-0.6*x).*y+1];

        alpha1=@(x,y) 0.2+0.6*x;                                            % volume fractions
        alpha2=@(x,y) 0.8-0.6*x;
        p=@(x,y) x+y+1;                                                     % pressure


end

U = @(x,y) [u1(x,y); u2(x,y); alpha1(x,y); alpha2(x,y); p(x,y)];            % collect all degrees of freedom

%% evaluate PDE

epsi=@(u) 1/2*(grad(u)+(grad(u)).');                                        % strain rate

% momentum equations
if strcmp(MPF_dragModel,'constant')                                         % constant drag term
    f1Fun=@(u1,u2,alpha1,alpha2,p) ...
        rho1*div(alpha1*u1*u1.') - 2*mu1*div(alpha1*epsi(u1)) + ...
        alpha1*grad(p) - dragCoef*(u2-u1);

    f2Fun=@(u1,u2,alpha1,alpha2,p) ...
        rho2*div(alpha2*u2*u2.') -2*mu2*div(alpha2*epsi(u2)) +...
        alpha2*grad(p) - dragCoef*(u1-u2);

else                                                                        % drag term based on Hiltunen
    f1Fun=@(u1,u2,alpha1,alpha2,p) ...
        rho1*div(alpha1*u1*u1.') - 2*mu1*div(alpha1*epsi(u1)) + ...
        alpha1*grad(p) - 3/4*Cd/d*alpha1*norm(u1-u2)*(u2-u1);

    f2Fun=@(u1,u2,alpha1,alpha2,p) ...
        rho2*div(alpha2*u2*u2.') -2*mu2*div(alpha2*epsi(u2)) +...
        alpha2*grad(p) - 3/4*Cd/d*alpha1*norm(u1-u2)*(u1-u2);

end
% continuity equations
f01Fun=@(u1,u2,alpha1,alpha2,p) ...
    div(alpha1*u1) ;
f02Fun=@(u1,u2,alpha1,alpha2,p) ...
    div(alpha2*u2) ;

% compute right-hand sides
[f1alpha1,f1Sym] = manufacturedSolution({u1,u2,alpha1,alpha2,p},f1Fun);             % momentum equations
[f2alpha2,f2Sym] = manufacturedSolution({u1,u2,alpha1,alpha2,p},f2Fun);
[f01,f01Sym]     = manufacturedSolution({u1,u2,alpha1,alpha2,p},f01Fun);            % continuity equations (should vanish)
[f02,f02Sym]     = manufacturedSolution({u1,u2,alpha1,alpha2,p},f02Fun);

f = @(x,y) [f1alpha1(x,y)./alpha1(x,y); f2alpha2(x,y)./alpha2(x,y);0;0;0 ]; % assemble in one vector, divide by alpha to get load

%% plotting
domain = [[0,1];[0,1]];                                                     % x- and y-limits for plotting
varNames = {'u_1^x','u_1^y','u_2^x','u_2^y','\alpha_1','\alpha_2','p'};     % names of each dof for plots
plotManufactured(U,f,domain,0,100,varNames);                                % call plot function