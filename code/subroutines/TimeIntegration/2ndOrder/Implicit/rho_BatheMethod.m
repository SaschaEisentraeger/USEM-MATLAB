%% Implicit time integration: rho-Bathe method
% 04/07/2019
% University of Magdeburg
% Sascha Duczek
%-------------------------------------------------------------------------%
% Input
% K:            Stiffness matrix
% M:            Mass matrix
% C:            Damping matrix
% F:            Force vector
% F_dt:     	Force vector (time-dependent)
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
% Noh, G. & Bathe, K.-J.: The Bathe Time Integration Method with 
% Controllable Spectral Radius: The rho_infinity-Bathe Method, Computers & 
% Structures, Elsevier BV, 2019, 212, 299-310
%
% We consider the Bathe implicit time integration method and focus on the
% time step splitting ratio and the spectral radius at large time steps to
% improve and generalize the scheme. The objective is to be able to
% prescribe the amplitude decay (dissipation) and period elongation
% (dispersion) for the numerical integration, and to achieve this aim in a
% direct and optimum manner with the minimum number of parameters. We show
% that the use of the time step splitting ratio and spectral radius is
% effective to prescribe in a smooth manner no amplitude decay to very
% large amplitude decays, with correspondingly small period elongation to
% very large period elongations while maintaining second-order accuracy. We
% analyze the effects of the splitting ratio and spectral radius on the
% stability and accuracy of the scheme and illustrate the use of these
% parameters in comparison with previously published methods. Furthermore,
% we show that with a proper setting of these parameters more accurate
% results may be obtained in some analyses.
%
% Noh, G. & Bathe, K.-J.: For Direct Time Integrations: A Comparison of the
% Newmark and rho_infinity-Bathe Schemes, Computers & Structures,
% Elsevier BV, 2019, 225, 106079
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U,V,A] = rho_BatheMethod(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
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

% Parameter for the time integration scheme
% This method includes as special cases also the standard Bathe method and
% the beta1-beta2-Bathe method.
% Increasing rho results in less amplitude decay and period elongation.
rho = 0.0;
% Largest numerical dissipation for rho = 0 (p = 2 - sqrt(2) = 0.586)
if rho == 1.0
	gamma = 1/2;
else
	gamma = @(x) (2-sqrt(2+2*x))/(1-x);
end
p = gamma(rho);
q1 = (rho+1)/(2*p*(rho-1)+4);
q0 = (p-1)*q1 + 0.5;
q2 = -p*q1 + 0.5;
a0 = 4/(p^2*dt^2);
a1 = 2/(p*dt);
a2 = 1/(dt^2*q2^2);
a3 = 1/(dt*q2);
a4 = 4/(p*dt);
a5 = (q0+q2)/(dt*q2^2);
a6 = q1/(dt*q2^2);
a7 = q0/q2;
a8 = q1/q2;

% Effective stiffness matrices
K_eff1 = K + a1*C + a0*M;
K_eff2 = K + a3*C + a2*M;

% The effective matrices are unsymmetric in general: LU decompositions
% (chosen by MATLAB), computed once and reused in every time step
dK_eff1 = decomposition(K_eff1);
dK_eff2 = decomposition(K_eff2);

% Result variables, starting with the initial values
uvec = zeros(NumTimeStep+1,numel(DOFsave));
vvec = zeros(NumTimeStep+1,numel(DOFsave));
avec = zeros(NumTimeStep+1,numel(DOFsave));
uvec(1,:) = x0(DOFsave);
vvec(1,:) = xd0(DOFsave);
avec(1,:) = xdd0(DOFsave);

pb = ProgressBar('rho-Bathe method (implicit)',NumTimeStep);

disp('Solve the equations of motion using the rho-Bathe method (implicit)');
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
    % Compute the effective load vector at the time t + p*dt
    F_eff1 = F_t_pdt + M*(a0*x0 + a4*xd0 + xdd0) + C*(a1*x0 + xd0);
 	% Solve the system of equations at the time t + p*dt to compute the
    % displacement: xp = K_eff1\F_eff1;
    xp = dK_eff1\F_eff1;
    % Compute the velocity at the time t + p*dt
    xdp  = a1*(xp - x0) - xd0;
	% Compute the acceleration at the time t + p*dt
    xddp = a1*(xdp - xd0) - xdd0;
    
    % Second part of the algorithm
    % Compute the effective load vector at the time t + dt
    F_eff2 = F_t_dt + ...
	         M*(a2*x0 + a5*xd0 + a6*xdp + a7*xdd0 + a8*xddp) + ...
			 C*(a3*x0 + a7*xd0 + a8*xdp);
    % Solve the system of equations at the time t + dt to compute the
    % displacement: x1 = K_eff2\F_eff2;
    x1 = dK_eff2\F_eff2;
    % Compute the velocity at the time t + dt
    xd1  = a3*(x1-x0) - a7*xd0  - a8*xdp;
    % Compute the acceleration at the time t + dt
    xdd1 = a3*(xd1-xd0) - a7*xdd0 - a8*xddp;
    
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
