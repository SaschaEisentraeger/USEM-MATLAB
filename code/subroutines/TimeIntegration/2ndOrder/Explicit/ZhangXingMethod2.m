%% Explicit time integration: Kim and Lee
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
function [T,U] = ZhangXingMethod2(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

% (De-)activate the waitbar
flag_wb = true;

% Number of active degrees of freedom
NumDoFActive = numel(ActiveDoFs);

% Number of degrees of freedom that are saved
NumDoFSave = numel(DOFsave);

% Test whether the stiffness matrix is a triangular matrix
if (istriu(K) || istril(K)) && ~isscalar(K)
    K = K + K';
    M = M + M';
    C = C + C';
end

if isstruct(F_dt)
	Fdt.Neumann = F_dt.Neumann;
	Fdt.K       = F_dt.K;
	Fdt.C       = F_dt.C;
	Fdt.M       = F_dt.M;
else
	Fdt.Neumann = F_dt;
	Fdt.K       = zeros(NumDoFActive,1);
	Fdt.C       = zeros(NumDoFActive,1);
	Fdt.M       = zeros(NumDoFActive,1);
end

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
rho = 0.0;

% Parameter
if rho ~= 1
    gamma = (2-sqrt(2*rho+2))/(2-2*rho);
    alpha = (2*rho-2)/(rho+1-sqrt(2*rho+2));
    beta = (rho+3-2*sqrt(2*rho+2))*(rho+1-sqrt(2*rho+2))/...
           ((1-rho)^2*(-2+sqrt(2*rho+2)));
else
    gamma = 1/4;
    alpha = 4;
    beta = 1/8;
end

% Stability limit
if rho == 1
    tauCr = @(x) 4*sqrt(4*x^2+1)-8*x;
else
    tauCr = @(x) (2*(1-rho)*sqrt((1-rho)^2*x^2-(rho+3-2*sqrt(2*rho+2))*...
                 (rho+1-2*sqrt(2*rho+2)))-(1-rho)^2*x)/...
                 (-(rho+3-2*sqrt(2*rho+2))*(rho+1-2*sqrt(2*rho+2)));
end

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
dtmax(2) = dtmax(2)*tau/2;
dt = min(dtmax);

% Number of time steps
NumTimeStep = ceil(round(tsim/dt,3,'decimal'));

% Save the full displacement field for a few time steps
nUAll = 0;
if nUAll > 0
    nSave = ceil(NumTimeStep/nUAll);
    UAll = cell(nUAll,1);
end

if flag_wb
	% Less updates of the waitbar (avoid slowing the simulation down)
	nUpdate = ceil(NumTimeStep/1000);
end

% Adjusted time increment
dt = tsim/NumTimeStep;

fprintf('Time step size used: %6.3e\n\n',dt);

% Initial conditions
% Displacement
x0  = U0;
% Velocity
xd0 = V0;

% Result variables
tvec = zeros(NumTimeStep+1,1);
uvec = zeros(NumTimeStep+1,NumDoFSave);
% Assign the initial conditions to the results
tvec(1) = 0;
uvec(1,:) = x0(DOFsave);

if flag_wb
	% Waitbar
	h = waitbar(0,'Zhang-Xing method 2 (explicit): ');
end

disp('Solve the equations of motion using the Zhang-Xing method 2 (explicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time integration scheme
for i = 1:NumTimeStep
    % Simulation time
    tE = i*dt;                  % Time at the end of the current interval
    tg = tE + (gamma-1)*dt;     % Time within the interval
    th = tE - 0.5*dt;           % Middle of the time step

    % Compute the force vector at the time t + gamma*dt
    Fg = F                                                             + ...
         Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,tg)       + ...
         Fdt.K*ExcitationSignal(SignalType,ExciteParam,tg)             + ...
         Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,tg,1) + ...
         Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,tg,2) ;
    % Compute the force vector at the time t + dt/2
    Fh = F                                                             + ...
         Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,th)       + ...
         Fdt.K*ExcitationSignal(SignalType,ExciteParam,th)             + ...
         Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,th,1) + ...
         Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,th,2) ;

    % Compute auxiliary displacement vectors
	if flag_diag
		xa = alpha*x0 + (1+alpha/2)*dt*xd0 + ...
			 dt^2*Minv.*(Fg - K*(x0 + gamma*dt*xd0) - C*xd0);
	else
		xa = alpha*x0 + (1+alpha/2)*dt*xd0 + ...
				 dt^2*dM\(Fg - K*(x0 + gamma*dt*xd0) - C*xd0);
	end
    xb = (1-alpha*beta)*x0 - 1/2*(alpha*beta+2*beta-1)*dt*xd0 + beta*xa;

    % Compute auxiliary velocity vector
    xda = 1/(2*dt)*(xa-alpha*x0) + (1/2-alpha/4)*xd0;

    % Compute the velocity vector
	if flag_diag
		xd1 = xd0 + dt*Minv.*(Fh - K*xb - C*xda);
	else
		xd1 = xd0 + dt*dM\(Fh - K*xb - C*xda);
	end

    % Compute the displacement vector
    x1 = 1/alpha*xa - (1/alpha-1/2)*dt*xd1;

    % Save variables
    tvec(i+1)   = tE;
    uvec(i+1,:) = x1(DOFsave);
	
	if nUAll > 0
        if mod(i,nSave) == 0
            c = i/nSave;
            UAll{c} = x1;
		elseif i == NumTimeStep
            UAll{nUAll} = x1;
        end
	end
	
    % Update the variables
    x0 = x1;
    xd0 = xd1;
    
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/(NumTimeStep+1),h,['Zhang-Xing method 2 (explicit): '...
			num2str(i/(NumTimeStep+1)*100,'%5.1f') '%'])
		end
	end
end

if flag_wb
	% Delete the waitbar
	close(h)
end

% Results
if NumDoFActive == NumDoFSave
    U = zeros(NumTimeStep+1,NumDoF);
    U(:,ActiveDoFs) = uvec;
    T = tvec;
else
    U = uvec;
    T = tvec;
end

if nUAll > 0
    save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
              filesep 'UAll.mat'],'-v7.3')
end

end