%% Implicit time integration: Newmark method
% 11.04.2017
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
% V:            Velocities
% A:			Accelerations
%-------------------------------------------------------------------------%
% Trapezoidal rule (average acceleration, beta = 1/4, gamma = 1/2):
% unconditionally stable, second-order accurate, no numerical damping.
% The calling convention and the shared helpers (TimeIntegration/Common)
% are described in TimeIntegrationMethods; this file is the template for
% new implicit methods.
%-------------------------------------------------------------------------%
function [T,U,V,A] = NewmarkMethod(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                     ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

% Time step (prescribed implicit step) and number of steps
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim);

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

% Initial conditions; acceleration from the equation of motion
x0  = U0;
xd0 = V0;
if isdiag(M)
    Minv = 1./diag(M);
    xdd0 = InitialAcceleration(@(r) Minv.*r,F,Fdt,K,C,x0,xd0,...
                               SignalType,ExciteParam);
else
    xdd0 = InitialAcceleration(@(r) M\r,F,Fdt,K,C,x0,xd0,...
                               SignalType,ExciteParam);
end

% Parameter for the time integration scheme:
beta_newmark = 1./4;
gamma_newmark = 1./2;

% Effective (modified) stiffness matrix
K_eff = K + 1/(beta_newmark*dt^2)*M + gamma_newmark/(beta_newmark*dt)*C;

% Parameter matrices
A1 = 1/(dt^2*beta_newmark)*M + gamma_newmark/(dt*beta_newmark)*C;
A2 = 1/(dt*beta_newmark)*M + (gamma_newmark/beta_newmark-1)*C;
A3 = (1/(2*beta_newmark)-1)*M + dt*(gamma_newmark/(2*beta_newmark)-1)*C;

% K_eff is unsymmetric in general: LU decomposition (chosen by MATLAB),
% computed once and reused in every time step
dK_eff = decomposition(K_eff);

% Result variables, starting with the initial values
uvec = zeros(NumTimeStep+1,numel(DOFsave));
vvec = zeros(NumTimeStep+1,numel(DOFsave));
avec = zeros(NumTimeStep+1,numel(DOFsave));
uvec(1,:) = x0(DOFsave);
vvec(1,:) = xd0(DOFsave);
avec(1,:) = xdd0(DOFsave);

pb = ProgressBar('Newmark method (implicit)',NumTimeStep);

disp('Solve the equations of motion using the Newmark method');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector and excitation signal
T = (0:NumTimeStep)'*dt;
[Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,T);

% Time integration scheme
for i = 2:NumTimeStep+1

    % Effective load at t + dt
    Fstar = LoadAtTime(F,Fdt,Amp(i),dAmp(i),ddAmp(i)) + ...
            A1*x0 + A2*xd0 + A3*xdd0;

    % Displacement at t + dt
    x1 = dK_eff\Fstar;

    % Velocity at t + dt
    xd1 = gamma_newmark/(beta_newmark*dt)*(x1-x0)   + ...
         (1-(gamma_newmark/beta_newmark))*xd0 	    + ...
         (1-gamma_newmark/(2*beta_newmark))*dt*xdd0;

    % Acceleration at t + dt
    xdd1 = 1/(beta_newmark*dt^2)*(x1-x0-xd0*dt)-(1/(2*beta_newmark)-1)*xdd0;

    % Save variables
    uvec(i,:) = x1(DOFsave);
    vvec(i,:) = xd1(DOFsave);
    avec(i,:) = xdd1(DOFsave);

    % Reset variables for the next time step
    x0   = x1;
	xd0  = xd1;
	xdd0 = xdd1;

    pb.update(i);
end
pb.close();

% Results (all DOFs if every active DOF was saved)
[U,V,A] = ExpandToAllDoFs(NumDoF,ActiveDoFs,uvec,vvec,avec);

end
