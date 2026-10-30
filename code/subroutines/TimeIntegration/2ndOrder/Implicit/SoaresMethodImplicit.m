%% Implicit time integration: Soares method
% 12.04.2017
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
% Soares, D.: A Simple and Effective New Family of Time Marching Procedures
% for Dynamics, Computer Methods in Applied Mechanics and Engineering,
% 2015, 283, 1138-1166
%
% In this work, a new family of time-marching procedures for dynamics is
% introduced. The proposed algorithm is based on displacements-velocities
% relations, requiring no computation of accelerations. Thus, not only it
% is a very simple technique, but also it is truly self-starting.
% Unconditionally stable implicit and conditionally stable explicit
% approaches are focused here, and optimized procedures are discussed for
% these analyses. The new technique is second-order accurate with
% controllable period elongation and algorithmic dissipation. Numerical
% results are presented, illustrating the good performance of the method.
% The new family is simple, efficient and highly accurate, resulting in a
% very effective time-marching technique.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = SoaresMethodImplicit(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

% Initial conditions
% Displacement
x0  = U0;
% Velocity
xd0 = V0;

% Determine the time increment used for the time integration (not via
% TimeStepSize: the parameters below use this step before its adjustment)
N_dt = numel(unique(dtmax));
if N_dt == 3
	% Deliberately described time step
	dt = dtmax(3);
else
	dt = dtmax(1);
end

% Parameter set
ParamSet = 1;

if ParamSet == 1
    % Parameter for the time integration scheme
    d = 1/2;
    a = 0.01; % controls the damping
    omega = sqrt(min(diag(K)./diag(M)));
    delta1 = 1/2 + d*omega*dt;
    alpha1 = 1/2 + 3/2*tanh(a*omega*dt);
    delta2 = 1 - delta1;
    alpha2 = 2*sqrt(2*alpha1) - alpha1 - 2*delta1;
%-------------------------------------------------------------------------%
elseif ParamSet == 2
    % We can also obtain an explicit algorithm for the following choices of
    % the integration parameter
    d = -1/4;
    omega = sqrt(min(diag(K)./diag(M)));
    delta1 = 1/2 + d*omega*dt;
    alpha1 = 0;
    delta2 = 1 - delta1;
    alpha2 = 2*delta2;
%-------------------------------------------------------------------------%
elseif ParamSet == 3
    % Test the high-order scheme by Soares
    delta1 = 1.0/2.0;
    alpha1 = 1.0/3.0;
    delta2 = 1.0/2.0;
    alpha2 = 2.0/3.0;

    % delta1 >= 1.0/2.0 and alpha1 < delta1 the method is only
    % conditionally stable and the critical time increment has to be
    % determined
    if delta1 >= 0.5 && alpha1 < delta1
        % Damping ratio
        xi = 0;
        % Maximum (angular) eigenfrequency --> computed from the critical
        % time step size of the centra difference method:
        % dt_CDM = 2/omega_max
        omega_max = 2/dtmax(4);
        % Critical time step size
        dtc = ((delta1 - 0.5)*xi + sqrt(0.5*(delta1-alpha1) + ...
               (delta1-0.5)^2*xi^2))/(0.5*(delta1-alpha1)*omega_max);
    else
        dtc = dtmax(1);
    end

    % Determine the time increment used for the time integration
    dt = min([dtc, dtmax(1:3)]);
end
%-------------------------------------------------------------------------%

% Number of time steps
NumTimeStep = ceil(round(tsim/dt,3,'decimal'));

% Adjusted time increment
dt = tsim/NumTimeStep;

fprintf('Time step size used: %6.3e\n\n',dt);

% Parameter for the interpolation of the force vector (trapezoidal rule)
beta1 = 0.5;
beta2 = 0.5;
% For beta1 = 1 and beta2 = 0 an explicit algorithm is also obtained for
% the load term.

% Effective (modified) stiffness matrix
K_eff = M + delta1*dt*C + 1/2*alpha1*dt^2*K;

% K_eff is unsymmetric in general: LU decomposition (chosen by MATLAB),
% computed once and reused in every time step
dK_eff = decomposition(K_eff);

% Result variables, starting with the initial displacements
uvec = zeros(NumTimeStep+1,numel(DOFsave));
uvec(1,:) = x0(DOFsave);

pb = ProgressBar('Soares method (implicit)',NumTimeStep);

disp('Solve the equations of motion using the Soares method (implicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector
T = (0:NumTimeStep)'*dt;

% Pre-compute the force vector
[Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,T);

% Time integration scheme
for i = 1:NumTimeStep
    % Compute the effective force vector at the time t + 1/2*dt
    F0 = LoadAtTime(F,Fdt,Amp(i),dAmp(i),ddAmp(i));
    F1 = LoadAtTime(F,Fdt,Amp(i+1),dAmp(i+1),ddAmp(i+1));
    Fstar = beta1*dt*F0 + beta2*dt*F1;
% 	Fstar = 0.5*dt*(F + Fdt*ExcitationSignal(SignalType,ExciteParam,t+0.5*dt));
    % Compute the effective load vector at the time t + dt
    F_eff = Fstar + M*xd0 - delta2*dt*C*xd0 - ...
            K*(dt*x0 + 1/2*alpha2*dt^2*xd0);
    % Solve the system of equations at the time t + dt to compute the
    % velocity
    xd1 = dK_eff\F_eff;
    % Compute the displacement at the time t + dt
    x1 = x0 + delta2*dt*xd0 + delta1*dt*xd1;
    
    % Save variables
    uvec(i+1,:) = x1(DOFsave);
    
    % Reset variables for the next time step
    x0   = x1;
	xd0  = xd1;
    
    pb.update(i);
end
pb.close();

% Results (all DOFs if every active DOF was saved)
U = ExpandToAllDoFs(NumDoF,ActiveDoFs,uvec);

end
