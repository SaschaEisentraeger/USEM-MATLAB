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
% F_dt          Force vector (time-dependent)
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
% Kim, W. & Lee, J. H.: An Improved Explicit Time Integration Method for
% Linear and Nonlinear Structural Dynamics,Computers & Structures,
% Elsevier BV, 2018, 206, 42-53
% 
% In this article, a new explicit time integration method is developed to
% analyze linear and nonlinear problems of structural dynamics. Like 
% recently developed explicit time integration methods, the new explicit
% method can also control the amount of numerical dissipation in the high
% frequency range. The method is explicit in the presence of the damping
% matrix, if the mass matrix is diagonal. Due to the unconventional
% approximations of the displacement vector, the new method does not
% require evaluation of the initial acceleration vector and other
% acceleration vectors. Linear and nonlinear problems of structural
% dynamics can be tackled in a consistent manner, and iterative solution
% finding procedures are not required. Various illustrative problems are
% used to investigate improved performance of the new explicit method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = KimMethodExplicit(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
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
% Specific procedures of the explicit method are given by the three
% parameters tau1, tau2, and alpha; The first variant shows no dissipation,
% the second alternative exhibits moderate dissipation (this setting is 
% recommended by the authors as the standard), while the third option
% introduces maximum dissipation.
IntMethod = 2;

% Set-up of the parameter of the time integrator
tau1_vec  = [0.25; 0.2684; 0.2831];
tau2_vec  = [0.5;  0.5;    0.5];
alpha_vec = [0.5;  0.4219; 0.337];

tau1  = tau1_vec(IntMethod);
tau2  = tau2_vec(IntMethod);
alpha = alpha_vec(IntMethod);

% Inverted mass matrix
if isdiag(M)
    flag_diag = true;
    Minv = 1./diag(M);
else
    flag_diag = false;
	% Matrix decomposition for solving linear systems
    dM = decomposition(M);
end

% Stability limit of the Noh-Bathe method compared to the CDM 
% (dt ~= 1.8725*dt_CDM)
dtmax(2) = 1.8725*dtmax(2);
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
	h = waitbar(0,'Kim method (explicit): ');
end

disp('Solve the equations of motion using the Kim method (explicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time integration scheme
for i = 1:NumTimeStep
    % Simulation time
    tE = i*dt;          % Time at the end of the current interval
    tS = tE - dt;       % Time at the start of the current interval
    
    % First sub-step
    x_tau1 = x0 + tau1*dt*xd0;
    xd_tau1 = xd0;
	ttau1 = tS+tau1*dt;
    ftau1 = F                                                                + ...
            Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,ttau1)       + ...
            Fdt.K*ExcitationSignal(SignalType,ExciteParam,ttau1)             + ...
            Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,ttau1,1) + ...
            Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,ttau1,2) ;
    % Compute the internal force vector
    ntau1 = C*xd_tau1 + K*x_tau1;
    % Compute the displacement vector for the first sub-step
	if flag_diag
		x1 = (alpha*dt)^2/2*Minv.*(ftau1 - ntau1) + x0 + alpha*dt*xd0;
	else
		x1 = (alpha*dt)^2/2*dM\(ftau1 - ntau1) + x0 + alpha*dt*xd0;
	end
    
    % Second sub-step
    x_tau2 = (tau2/alpha)^2*x1 + ...
             ((alpha-tau2)*(alpha+tau2))/alpha^2*x0 + ...
             (tau2*(alpha-tau2)*dt)/alpha*xd0;
    xd_tau2 = (2*tau2)/(alpha^2*dt)*(x1-x0) + (alpha-2*tau2)/alpha*xd0;
    ttau2 = tS+tau2*dt;
    ftau2 = F                                                                + ...
            Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,ttau2)       + ...
            Fdt.K*ExcitationSignal(SignalType,ExciteParam,ttau2)             + ...
            Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,ttau2,1) + ...
            Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,ttau2,2) ;
    % Compute the internal force vector
    ntau2 = C*xd_tau2 + K*x_tau2;
    % Compute the displacement and velocity vectors for the second sub-step
	if flag_diag
		x2  = (alpha-1)*dt^2/(2*(alpha-3*tau2))*Minv.*(ftau2 - ntau2) - ...
			  (3*tau2-1)/(alpha^2*(alpha-3*tau2))*x1 + ...
			  ((alpha-1)*(alpha^2+alpha+1-3*tau2*alpha-3*tau2))/...
			  (alpha^2*(alpha-3*tau2))*x0 + ...
			  ((alpha-1)*(alpha+1-3*tau2)*dt)/(alpha*(alpha-3*tau2))*xd0;
	else
		x2  = (alpha-1)*dt^2/(2*(alpha-3*tau2))*dM\(ftau2 - ntau2) - ...
			  (3*tau2-1)/(alpha^2*(alpha-3*tau2))*x1 + ...
			  ((alpha-1)*(alpha^2+alpha+1-3*tau2*alpha-3*tau2))/...
			  (alpha^2*(alpha-3*tau2))*x0 + ...
			  ((alpha-1)*(alpha+1-3*tau2)*dt)/(alpha*(alpha-3*tau2))*xd0;
	end
    xd2 = (2*alpha-3)/(dt*(alpha-1))*x2 + ...
           1/(alpha^2*(alpha-1)*dt)*x1 - ...
           ((2*alpha+1)*(alpha-1))/(alpha^2*dt)*x0 - ...
           (alpha-1)/alpha*xd0;
    
    % Save variables
    tvec(i+1)   = tE;
    uvec(i+1,:) = x2(DOFsave);
	
	if nUAll > 0
        if mod(i,nSave) == 0
            c = i/nSave;
            UAll{c} = x2;
		elseif i == NumTimeStep
            UAll{nUAll} = x2;
        end
	end
	
    % Update the variables
    x0 = x2;
    xd0 = xd2;
    
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/NumTimeStep,h,['Kim method (explicit): '...
			num2str(i/NumTimeStep*100,'%5.1f') '%'])
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