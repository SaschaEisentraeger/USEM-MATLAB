%% Implicit time integration: Bathe-Baig method
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bathe, K. J. & Baig, M. M. I.: On a Composite Implicit Time Integration
% Procedure for Nonlinear Dynamics, Computers & Structures, 2005, 83,
% 2513-2524
%
% Transient analysis of nonlinear problems in structural and solid
% mechanics is mainly carried out using direct time integration of the
% equations of motion. For reliable solutions, a stable and efficient
% integration algorithm is desirable. Methods that are unconditionally
% stable in linear analyses appear to be a natural choice for use in
% nonlinear analyses, but unfortunately may not remain stable for a given
% time step size in large deformation and long time range response 
% solutions. A composite time integration scheme is proposed and tested in
% some example solutions against the trapezoidal rule and the Wilson
% theta-method, and found to be effective where the trapezoidal rule fails
% to produce a stable solution. These example results are indicative of the
% merits of the composite scheme.
%
% Bathe, K. J.: Conserving Energy and Momentum in Nonlinear Dynamics: A
% Simple Implicit Time Integration Scheme, Computers & Structures, 2007,
% 85, 437-445
%
% Bathe, K. J. & Noh, G.: Insight into an Implicit Time Integration Scheme
% for Structural Dynamics, Computers & Structures, 2012, 98-99, 1-6
%
% Noh, G.; Ham, S. & Bathe, K. J.: Performance of an Implicit Time
% Integration Scheme in the Analysis of Wave Propagation, Computers &
% Structures, 2013, 123, 93-105
%
% Noh, G. & Bathe, K. J.: Further insights into an implicit time
% integration scheme for structural dynamics, Computers & Structures, 2018,
% 202, 15-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U,V,A] = BatheBaigMethod(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                     ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

% Use in the first time step a different method to avoid overshoots
% == 1 --> yes; == 0 --> no
flag_overshoot = 0;

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

% Parameter for the time integration scheme
if flag_overshoot == 1
    p = 0.5;
else
    p  = 2 - sqrt(2); % for linear structural dynamics
end
% p  = 0.731; % for geometrically nonlinear structural dynamics

% Auxiliary variables
if flag_overshoot == 1
    alpha = 1;
    delta = 0.75;
    a11 = 1/(alpha*p^2*dt^2);
    a21 = 1/(alpha*p*dt);
    a31 = 1/(2*alpha)-1;
    a41 = delta/(alpha*p*dt);
    a51 = delta/alpha-1;
    a61 = ((delta-2*alpha)*p*dt)/(2*alpha);
end
a1n = 4/(p^2*dt^2);
a2n = 4/(p*dt);
a3n = 2/(p*dt);
c1 = (1-p)/(p*dt);
c2 = -1/((1-p)*p*dt);
c3 = (2-p)/((1-p)*dt);

% Effective stiffness matrices
if flag_overshoot == 1
    K_eff11 = K + a11*M + a41*C;
end
K_eff1n = K + a1n*M + a3n*C;
K_eff2 = K + c3*c3*M + c3*C;

% The effective matrices are unsymmetric in general: LU decompositions
% (chosen by MATLAB), computed once and reused in every time step
dK_eff1n = decomposition(K_eff1n);
dK_eff2 = decomposition(K_eff2);
if flag_overshoot == 1
    dK_eff11 = decomposition(K_eff11);
end

% Result variables, starting with the initial values
uvec = zeros(NumTimeStep+1,numel(DOFsave));
vvec = zeros(NumTimeStep+1,numel(DOFsave));
avec = zeros(NumTimeStep+1,numel(DOFsave));
uvec(1,:) = x0(DOFsave);
vvec(1,:) = xd0(DOFsave);
avec(1,:) = xdd0(DOFsave);

pb = ProgressBar('Bathe-Baig method (implicit)',NumTimeStep);

disp('Solve the equations of motion using the Bathe-Baig method (implicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector
T  = (0:NumTimeStep)'*dt;
% Steps in between
Tp = (0:NumTimeStep)'*dt - (1-p)*dt;

% Pre-compute the amplitude values for the force vectors
[Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,T);
[Amp_p,dAmp_p,ddAmp_p] = SignalHistory(SignalType,ExciteParam,Tp);

% Time integration scheme
for i = 1:NumTimeStep
    % First part of the algorithm
    % Compute the load vector at the time t + p*dt
    F_t_pdt = LoadAtTime(F,Fdt,Amp_p(i+1),dAmp_p(i+1),ddAmp_p(i+1));
    % Compute the load vector at the time t + dt
    F_t_dt  = LoadAtTime(F,Fdt,Amp(i+1),dAmp(i+1),ddAmp(i+1));
% 	% Compute the load vector at the time t
% 	F_t = LoadAtTime(F,Fdt,Amp(i),dAmp(i),ddAmp(i));
% 	% Linear interpolation of the force vector between time steps
% 	F_t_pdt = (1-p)*F_t + p*F_t_dt;
    if flag_overshoot == 1 && i == 1
        % Compute the effective load vector at the time t + p*dt
        F_eff1 = F_t_pdt + M*(a11*x0 + a21*xd0 + a31*xdd0) + ...
                           C*(a41*x0 + a51*xd0 + a61*xdd0);
        % Solve the system of equations at the time t + p*dt to compute the
        % displacement: xp = K_eff1\F_eff1;
        xp = dK_eff11\F_eff1;
        % Compute the velocity at the time t + p*dt
        xdp = xd0 + p*dt*(1-delta)*xdd0 + ...
                    delta*p*dt*(a11*(xp - x0) - a21*xd0 - a31*xdd0);
    else
        % Compute the effective load vector at the time t + p*dt
        F_eff1 = F_t_pdt + M*(a1n*x0 + a2n*xd0 + xdd0) + C*(a3n*x0 + xd0);
        % Solve the system of equations at the time t + p*dt to compute the
        % displacement: xp = K_eff1\F_eff1;
        xp = dK_eff1n\F_eff1;
        % Compute the velocity at the time t + p*dt
        xdp = (xp - x0)*a3n - xd0;
    end
    
    % Second part of the algorithm
    % Compute the effective load vector at the time t + dt
    F_eff2 = F_t_dt - M*(c1*xd0 + c2*xdp + c3*c1*x0 + c3*c2*xp) - ...
             C*(c1*x0 + c2*xp);
    % Solve the system of equations at the time t + dt to compute the
    % displacement: x1 = K_eff2\F_eff2;
    x1 = dK_eff2\F_eff2;
    % Compute the velocity at the time t + dt
    xd1 = c1*x0 + c2*xp + c3*x1;
    % Compute the acceleration at the time t + dt
    xdd1 = c1*xd0 + c2*xdp + c3*xd1;
    
    % Save variables
    uvec(i+1,:) = x1(DOFsave);
    vvec(i+1,:) = xd1(DOFsave);
    avec(i+1,:) = xdd1(DOFsave);
	
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
