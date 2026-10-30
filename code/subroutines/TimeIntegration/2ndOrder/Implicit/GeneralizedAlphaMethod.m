%% Implicit time integration: Generalized-alpha method
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
% Chung, J. & Hulbert, G. M.: A Time Integration Algorithm for Structural
% Dynamics with Improved Numerical Damping: The Generalized-alpha Method 
% Journal of Applied Mechanics, 1993, 60, 371-375
%
% A new family of time integration algorithms is presented for solving
% structural dynamics problems. The new method, denoted as the generalized-
% alpha method, possesses numerical dissipation that can be controlled by
% the user. In particular, it is shown that the generalized-a method
% achieves high-frequency dissipation while minimizing unwanted low-
% frequency dissipation. Comparisons are given of the generalized-alpha
% method with other numerically dissipative time integration methods; these
% results highlight the improved performance of the new algorithm. The new
% algorithm can be easily implemented into programs that already include
% the Newmark and Hilber-Hughes-Taylor-alpha time integration methods.

% Implementation:
% Schauer, M. M.: Ein effizienter gekoppelter FEM-SBFEM Ansatz zur Analyse
% von Boden-Bauwerk-Interaktionen im Zeitbereich, Technical University
% Carolo-Wilhelmina of Braunschweig, 2015.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U,V,A] = GeneralizedAlphaMethod(K,M,C,F,F_dt,U0,V0,NumDoF,...
                     DOFsave,ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

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
% If we choose alpha_m to be zero and 
% alpha_f = (1 - p_infinity)/(1 + p_infinity), we obtain the HHT-alpha 
% method. On the other hand, if alpha_f is set to zero  and 
% alpha_m = (p_infinity - 1)/(p_infinity + 1) the generalized alpha method 
% reproduces the WBZ-alpha method. However, if both alpha_m and alpha_f are
% set to zero the generalized alpha method is identical to the standard 
% Newmark method.
% p_infinity: User specified algorithmic dissipation control parameter
% p_infinity = 1 --> no dissipation
% p_infinity = 0 --> asymptotic annihilation (the high frequency response
% is annihilated after one time step; maximum algorithmic dissipation)
p_infinity    =  0.8;
alpha_m       = (2*p_infinity - 1)/(p_infinity + 1);
alpha_f       =  p_infinity/(p_infinity + 1);
beta_newmark  =  1.0/4.0*(1 - alpha_m + alpha_f)^2;
gamma_newmark =  1.0/2.0 - alpha_m + alpha_f;

% Effective (modified) stiffness matrix
K_eff = (1 - alpha_m)/(beta_newmark*dt^2)*M + ...
        (1 - alpha_f)*gamma_newmark/(beta_newmark*dt)*C + ...
        (1 - alpha_f)*K;
    
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

pb = ProgressBar('Generalized-alpha method (implicit)',NumTimeStep);

disp('Solve the equations of motion using the generalized-alpha method (implicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector
T = (0:NumTimeStep)'*dt;
% Steps in between (time vector)
Talpha = (0:NumTimeStep)'*dt - alpha_f*dt;

% Pre-compute the amplitude values for the force vectors
[Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,Talpha);

% Time integration scheme
for i = 1:NumTimeStep
	% The force vector is computed according to the original publication by
	% Chung & Hulbert
    f = LoadAtTime(F,Fdt,Amp(i+1),dAmp(i+1),ddAmp(i+1));
    % Effective force vector
    Fstar = f + (1 - alpha_m)*M*(1/(beta_newmark*dt^2)*x0            + ...
             1/(beta_newmark*dt)*xd0                                 + ...
            (1/(2*beta_newmark) - 1)*xdd0)                           + ...
            (1 - alpha_f)*C*(gamma_newmark/(beta_newmark*dt)*x0      + ...
            (gamma_newmark/beta_newmark - 1)*xd0                     + ...
            (gamma_newmark*dt/(2*beta_newmark) - dt)*xdd0)           - ...
            alpha_m*M*xdd0                                           - ...
            alpha_f*C*xd0                                            - ...
            alpha_f*K*x0;
        
    % Solve the system of equations at the time t + dt to compute the
    % displacement
    x1 = dK_eff\Fstar;
    
	% The velocity and acceleration vectors are updated according to the
    % Newmark method: Compute the velocity at the time t + dt
    xd1 = gamma_newmark/(beta_newmark*dt)*(x1-x0) + ...
         (1-(gamma_newmark/beta_newmark))*xd0 	  + ...
         (1-gamma_newmark/(2*beta_newmark))*dt*xdd0;

  	% Compute the acceleration at the time t + dt
    xdd1 = 1/(beta_newmark*dt^2)*(x1-x0-xd0*dt)-(1/(2*beta_newmark)-1)*xdd0;
    
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
