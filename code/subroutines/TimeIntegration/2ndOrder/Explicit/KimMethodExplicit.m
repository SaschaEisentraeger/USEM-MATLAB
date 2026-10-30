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

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

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
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim,1.8725);

% Initial conditions (the method needs no initial acceleration)
x0  = U0;
xd0 = V0;

% Result variables
tvec = zeros(NumTimeStep+1,1);
uvec = zeros(NumTimeStep+1,numel(DOFsave));
% Assign the initial conditions to the results
tvec(1) = 0;
uvec(1,:) = x0(DOFsave);

pb = ProgressBar('Kim method (explicit)',NumTimeStep);

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
    [s,ds,dds] = SignalHistory(SignalType,ExciteParam,ttau1);
    ftau1 = LoadAtTime(F,Fdt,s,ds,dds);
    % Compute the internal force vector
    ntau1 = C*xd_tau1 + K*x_tau1;
    % Compute the displacement vector for the first sub-step
	if flag_diag
		x1 = (alpha*dt)^2/2*Minv.*(ftau1 - ntau1) + x0 + alpha*dt*xd0;
	else
		x1 = (alpha*dt)^2/2*(dM\(ftau1 - ntau1)) + x0 + alpha*dt*xd0;
	end
    
    % Second sub-step
    x_tau2 = (tau2/alpha)^2*x1 + ...
             ((alpha-tau2)*(alpha+tau2))/alpha^2*x0 + ...
             (tau2*(alpha-tau2)*dt)/alpha*xd0;
    xd_tau2 = (2*tau2)/(alpha^2*dt)*(x1-x0) + (alpha-2*tau2)/alpha*xd0;
    ttau2 = tS+tau2*dt;
    [s,ds,dds] = SignalHistory(SignalType,ExciteParam,ttau2);
    ftau2 = LoadAtTime(F,Fdt,s,ds,dds);
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
		x2  = (alpha-1)*dt^2/(2*(alpha-3*tau2))*(dM\(ftau2 - ntau2)) - ...
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

    % Update the variables
    x0 = x2;
    xd0 = xd2;
    
    pb.update(i);
end
pb.close();

% Results (all DOFs if every active DOF was saved)
U = ExpandToAllDoFs(NumDoF,ActiveDoFs,uvec);
T = tvec;

end