%% Explicit time integration: Zhang and Xing
% 04.09.2018
% University of Magdeburg
% Sascha Duczek
%-------------------------------------------------------------------------%
% Input
% K:            Stiffness matrix
% M:            Mass matrix
% C:            Damping matrix
% F:            Force vector
% Fdt:          Force vector (time-dependent)
% U0:			Initial displacement vector
% V0:			Initial velocity vector
% NumDoF:		Number of degrees of freedom
% DOFsave:		Degrees of freedom for which the results are saved
% ActiveDoFs:   Active degrees of freedom (without Dirichlet boundary
%               conditions)
% ExciteParam:  Parameters describing the excitation signal
% SignalType:   Type of the excitation signal
% dtmax:        Vector of different allowed maximum time increments
%               [dt_fex,dt_CDM,dt_imp]; the critical time increment of the
% 				central difference method (CDM) is estimated by 2/omega_max
% tsim:         Simulation time
%-------------------------------------------------------------------------%
% Output
% T:            Time vector
% U:            Displacements
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zhang, H. M. & Xing, Y. F.: Two Novel Explicit Time Integration Methods
% based on Displacement-Velocity Relations for Structural Dynamics,
% Computers & Structures, Elsevier BV, 2019, 221, 127-141
%
% Two novel explicit time integration methods are proposed based on
% displacement-velocity relations in this paper for structural dynamics.
% They avoid the factorization of damping and stiffness matrices, and are
% truly self-starting due to the exclusion of acceleration vectors. The
% first method employs the motion equation of expanded form and a linear
% relation of the displacement and velocity vectors. The recommended
% parameters, derived from linear analysis, enable the method to possess
% first-order accuracy mostly and second-order accuracy in the absence of
% numerical and physical damping. The second method adopts the idea of
% composite methods, and employs two motion equations of different expanded
% forms per step. Theoretical analysis indicates that this method can
% achieve a maximum stability limit of 4, and provides a single-parameter
% optimal scheme, controlled by the degree of numerical dissipation, for
% this method. The resulting scheme is second-order accurate with the
% stability limit ranged from 3.5708 to 4 for the undamped case, and some
% numerical experiments show that it has better numerical performance
% compared with some up-to-date explicit methods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = ZhangXingMethod1(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

%% Explicit time integration scheme
% Inverted mass matrix
if isdiag(M)
    flag_diag = true;
    Minv = 1./diag(M);
else
    flag_diag = false;
	% Matrix decomposition for solving linear systems
    dM = decomposition(M);
end

% Degree of numerical dissipation
rho = 1.0;

% Parameter
gamma = (3-rho)/(2+2*rho);

% Stability limit
tauCr = @(x) (-x + sqrt(x^2 + 2*gamma))/gamma;

if nnz(C) == 0
    zeta = 0;
elseif nnz(C) ~= 0
    % The damping ratio zeta is frequency-dependent for Rayleigh damping an
    % therefore, it is unclear how to set zeta; in the article the authors
    % used zeta = 0.1 but did not comment on the use in the FEM
    zeta = 0.1;
end
tau = tauCr(zeta);

% Determine the time increment used for the time integration:
% 2 is the stability limit for the CDM and tau is the corresponding limit
% for the novel explicit scheme.
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim,tau/2);

% Initial conditions (the method is self-starting: no initial
% acceleration needed)
x0  = U0;
xd0 = V0;

% Result variables
tvec = zeros(NumTimeStep+1,1);
uvec = zeros(NumTimeStep+1,numel(DOFsave));
% Assign the initial conditions to the results
tvec(1) = 0;
uvec(1,:) = x0(DOFsave);

pb = ProgressBar('Zhang-Xing method 1 (explicit)',NumTimeStep);

disp('Solve the equations of motion using the Zhang-Xing method 1 (explicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time integration scheme
for i = 1:NumTimeStep
    % Simulation time
    tE = i*dt;                  % Time at the end of the current interval
    tg = tE + (gamma-1)*dt;     % Time within the interval

    % Compute the force vector at the time t + gamma*dt
    [s,ds,dds] = SignalHistory(SignalType,ExciteParam,tg);
    Fg = LoadAtTime(F,Fdt,s,ds,dds);

    % Compute the displacement vector
	if flag_diag
		x1 = x0 + dt*xd0 + dt^2/2*Minv.*(Fg - K*(x0 + gamma*dt*xd0) - C*xd0);
	else
		x1 = x0 + dt*xd0 + dt^2/2*(dM\(Fg - K*(x0 + gamma*dt*xd0) - C*xd0));
	end

    % Compute the velocity vector
    xd1 = 2/dt*(x1 - x0) - xd0;

    % Save variables
    tvec(i+1)   = tE;
    uvec(i+1,:) = x1(DOFsave);

    % Update the variables
    x0 = x1;
    xd0 = xd1;
    
    pb.update(i);
end
pb.close();

% Results (all DOFs if every active DOF was saved)
U = ExpandToAllDoFs(NumDoF,ActiveDoFs,uvec);
T = tvec;

end