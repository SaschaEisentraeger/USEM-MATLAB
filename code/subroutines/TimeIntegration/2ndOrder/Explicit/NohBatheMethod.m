%% Explicit time integration: Noh-Bathe method
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
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Noh, G. & Bathe, K. J.: An Explicit Time Integration Scheme for the
% Analysis of Wave Propagations, Computers & Structures, 2013, 129, 178-193
%
% A new explicit time integration scheme is presented for the solution of
% wave propagation problems. The method is designed to have small solution
% errors in the frequency range that can spatially be represented and to
% cut out high spurious frequencies. The proposed explicit scheme is
% second-order accurate for systems with and without damping, even when
% used with a non-diagonal damping matrix. The stability, accuracy and
% numerical dispersion are analyzed, and solutions to problems are given
% that illustrate the performance of the scheme.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = NohBatheMethod(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,....
                 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

% Inverted mass matrix
if isdiag(M)
    flag_diag = true;
    Minv = 1./diag(M);
else
    flag_diag = false;
	% Matrix decomposition for solving linear systems
    dM = decomposition(M);
end

% Initial conditions; acceleration from the equation of motion
x0  = U0;
xd0 = V0;
if flag_diag
    xdd0 = InitialAcceleration(@(r) Minv.*r,F,Fdt,K,C,x0,xd0,...
                               SignalType,ExciteParam);
else
    xdd0 = InitialAcceleration(@(r) dM\r,F,Fdt,K,C,x0,xd0,...
                               SignalType,ExciteParam);
end

% Parameter for the time integration scheme
p  = 0.54; % for p = 1 the Noh-Bathe method is identical to the CDM
q1 = (1-2*p)/(2*p*(1-p));
q2 = 1/2-p*q1;
q0 = -q1-q2+1/2;

% Stability limit of the Noh-Bathe method compared to the CDM 
% (dt ~= 1.8725*dt_CDM)
gamma = 1/4 - 1/2*(1-p)*q1;
val = 1/2*sqrt(1/(gamma*p*(1-p)));
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim,val);

% Parameter for the time integration scheme
a0 = p*dt;
a1 = 1/2*(p*dt)^2;
a2 = a0/2;
a3 = (1-p)*dt;
a4 = 1/2*((1-p)*dt)^2;
a5 = q0*a3;
a6 = (1/2+q1)*a3;
a7 = q2*a3;

% Result variables
uvec = zeros(NumTimeStep+1,numel(DOFsave));
% Assign the initial conditions to the results
uvec(1,:) = x0(DOFsave);

pb = ProgressBar('Noh-Bathe method (explicit)',NumTimeStep);

disp('Solve the equations of motion using the Noh-Bathe method (explicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector
T  = (0:NumTimeStep)'*dt;
% Steps in between
Tp = (0:NumTimeStep)'*dt - (1-p)*dt;

% Pre-compute the amplitude values for the force vectors (at the time
% steps and at the intermediate times t + p*dt)
[Amp,dAmp,ddAmp]       = SignalHistory(SignalType,ExciteParam,T);
[Amp_p,dAmp_p,ddAmp_p] = SignalHistory(SignalType,ExciteParam,Tp);

% Time integration scheme
for i = 1:NumTimeStep
	% First part of the algorithm
	% Compute the displacement at the time t + p*dt
    xp = x0 + a0*xd0 + a1*xdd0;
    % Compute the load vector at the time t + p*dt
    F_t_pdt = LoadAtTime(F,Fdt,Amp_p(i+1),dAmp_p(i+1),ddAmp_p(i+1));
    % Compute the load vector at the time t + dt
    F_t_dt  = LoadAtTime(F,Fdt,Amp(i+1),dAmp(i+1),ddAmp(i+1));
%     % Compute the load vector at the time t
%     F_t = LoadAtTime(F,Fdt,Amp(i),dAmp(i),ddAmp(i));
%     % Linear interpolation of the force vector between time steps
%     F_t_pdt = (1-p)*F_t + p*F_t_dt;
    % Compute the effective load vector at the time t + p*dt
    F_t_pdt_circ    = F_t_pdt - K*xp - C*(xd0 + a0*xdd0);
    % Solve the system of equations at the time t + p*dt to compute the
    % acceleration
    if flag_diag
        xddp = Minv.*F_t_pdt_circ;
    else
        xddp = dM\F_t_pdt_circ;
    end
    % Compute the velocity at the time t + p*dt
    xdp = xd0 + a2*(xdd0 + xddp);
    
    % Second part of the algorithm
    % Compute the displacement at the time t + dt
    x1 = xp + a3*xdp + a4*xddp;
    % Compute the load vector at the time t + dt
    F_t_dt_circ = F_t_dt - K*x1 - C*(xdp + a3*xddp);
 	% Solve the system of equations at the time t + dt to compute the
    % acceleration
    if flag_diag
        xdd1 = Minv.*F_t_dt_circ;
    else
        xdd1 = dM\F_t_dt_circ;
    end
    % Compute the velocity at the time t + dt
    xd1 = xdp + a5*xdd0 + a6*xddp + a7*xdd1;
    
 	% Save variables
    uvec(i+1,:) = x1(DOFsave);
    
    % Reset variables for the next time step
    x0   = x1;
	xd0  = xd1;
	xdd0 = xdd1;
    
    pb.update(i);
end
pb.close();

% Results (all DOFs if every active DOF was saved)
U = ExpandToAllDoFs(NumDoF,ActiveDoFs,uvec);

end