%% Implicit time integration: Soares method
% 12/08/2020
% UNSW Sydney
% Sascha Eisentraeger
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
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Soares, D.: A straightforward high-order accurate time-marching procedure
% for dynamic analyses, Engineering with Computers, 2020
% https://doi.org/10.1007/s00366-020-01129-1
%
% In this paper, a novel time-marching procedure is proposed for dynamic
% analysis. The proposed new technique is fourth-order accurate for
% undamped models; third-order accurate for numerically damped models;
% and second-order accurate for physically damped models. The method is
% very straightforward to implement and to apply, allowing providing
% higher-order accurate time-domain solutions taking into account basically
% the same computational effort of standard single-step first- or second-
% order accurate time-marching procedures. Numerical examples are presented
% at the end of the paper, illustrating the excellent accuracy of the novel
% approach.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U,V] = SoaresMethodImplicit_ho(K,M,C,F,F_dt,U0,V0,NumDoF,...
                   DOFsave,ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

% Time integration parameter:
% alpha = 1 (no numerical damping/dissipation is introduced)
% alpha > 1 (numerical damping/dissipation is introduced)
alpha = 1;

% Critical sampling frequency
Omega_c = sqrt(2/alpha*(sqrt(1+12*alpha)-1));

% Determine the time increment used for the time integration: the method
% is only conditionally stable, dt <= Omega_c/omega_max, i.e. Omega_c/2
% times the stability limit 2/omega_max of the CDM (dtmax(2))
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim,Omega_c/2);

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

% Initial conditions
% Displacement
x0  = U0;
% Velocity
xd0 = V0;

% Effective (modified) stiffness matrix
K_eff = M + 1.0/2.0*dt*C + 1.0/6.0*dt^2*K;

% K_eff is unsymmetric in general: LU decomposition (chosen by MATLAB),
% computed once and reused in every time step
dK_eff = decomposition(K_eff);

% Determine if the mass matrix is lumped (diagonal)
if isdiag(M)
    flag_lumped = 1;
    M = full(diag(M));
    Minv = 1./M;
else
    flag_lumped = 0;
end

% Auxiliary vectors
if flag_lumped == 1
    R0 = K_eff*xd0;
    S0 = 1.0/2.0*dt*C*xd0;
    V0 = 1.0/2.0*Minv.*(R0 - S0) - 1.0/2.0*xd0;
else
    R0 = K_eff*xd0;
    S0 = 1.0/2.0*dt*C*xd0;
    V0 = 1.0/2.0*(M\(R0 - S0)) - 1.0/2.0*xd0;
end

% Result variables
uvec = zeros(NumTimeStep+1,numel(DOFsave));
vvec = zeros(NumTimeStep+1,numel(DOFsave));
% Assign the initial conditions (displacements, velocities) to the results
uvec(1,:) = x0(DOFsave);
vvec(1,:) = xd0(DOFsave);

pb = ProgressBar('High-order Soares method (implicit)',NumTimeStep);

disp('Solve the equations of motion using the high-order Soares method (implicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector and excitation signal
T = (0:NumTimeStep)'*dt;
[Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,T);

% Time integration scheme
for i = 1:NumTimeStep
    % Compute the effective force vector at the time t + 1/2*dt
    % (trapezoidal rule for the load integral over the time step)
    F0 = LoadAtTime(F,Fdt,Amp(i),dAmp(i),ddAmp(i));
    F1 = LoadAtTime(F,Fdt,Amp(i+1),dAmp(i+1),ddAmp(i+1));
    Fbar = 1.0/2.0*dt*F0 + 1.0/2.0*dt*F1;
	
	% Compute the effective load vector at the time t + dt
    if flag_lumped == 1
        R1 = Fbar - S0 + M.*xd0 - ...
             K*(dt*x0 + 1.0/3.0*dt^2*xd0 + 1.0/2.0*alpha*dt^2*V0);
    else
        R1 = Fbar - S0 + M*xd0 - ...
             K*(dt*x0 + 1.0/3.0*dt^2*xd0 + 1.0/2.0*alpha*dt^2*V0);
    end
			
    % Solve the system of equations at the time t + dt to compute the
    % velocity
    xd1 = dK_eff\R1;

	% Auxiliary vectors at the time t + dt
	S1 = 1.0/2.0*dt*C*xd1;
    if flag_lumped == 1
        V1 = 1.0/2.0*Minv.*(R1 - S1) - 1.0/2.0*xd1;
    else
        V1 = 1.0/2.0*(M\(R1 - S1)) - 1.0/2.0*xd1;
    end
	
    % Compute the displacement at the time t + dt
    x1 = x0 + 1.0/2.0*dt*(xd0 + xd1 + V0 + V1);
    
    % Save variables
    uvec(i+1,:) = x1(DOFsave);
    vvec(i+1,:) = xd1(DOFsave);

    % Reset variables for the next time step
	S0 = S1;
	V0 = V1;
    x0  = x1;
	xd0 = xd1;

    pb.update(i);
end
pb.close();

% Results (all DOFs if every active DOF was saved)
[U,V] = ExpandToAllDoFs(NumDoF,ActiveDoFs,uvec,vvec);

end