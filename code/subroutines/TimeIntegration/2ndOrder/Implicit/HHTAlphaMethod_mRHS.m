%% Implicit time integration: HHT-alpha method
%__________________________________________________________________________
% Numerical solution of the equation of motion by means of the HHT-alpha
% method with multiple right-hand sides (mRHS), i.e., multiple load cases
% are simultaneously considered in the simulation.
%__________________________________________________________________________
%
% 20/04/2022
% UNSW Sydney
% Sascha Eisentraeger
%-------------------------------------------------------------------------%
% Input
% All input matrices and vectors contain only active DOFs! Reduced system!
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
% Hilber, H. M. & Hughes, T. J. R. & Taylor, R. L.: Improved Numerical
% Dissipation for Time Integration Algorithms in Structural Dynamics 
% Earthquake Engineering and Structural Dynamics, 1977, 5, 283-292
%
% A new family of unconditionally stable one-step methods for the direct
% integration of the equations of structural dynamics is introduced and is
% shown to possess improved algorithmic damping properties which can be
% continuously controlled. The new methods are compared with members of the
% Newmark family, and the Houbolt and Wilson methods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U,V,A] = HHTAlphaMethod_mRHS(K,M,C,F,F_dt,U0,V0,NumDoF,...
                     DOFsave,ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

% Parameter for the time integration scheme: -1/3 <= alpha <= 0
% If alpha = 0 --> Newmark method!
% If alpha = -0.05 --> Default value in Abaqus/Standard!
alpha = -0.05;
beta_newmark = (1-alpha)^2/4;
gamma_newmark = 1./2-alpha;

% Time step (prescribed implicit step) and number of steps
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim);

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

% Check for multiple load cases
nLC = size(Fdt.Neumann,2);

% Adjust the size of the force vectors
nF = size(F,2);
if nF ~= nLC
    F = repmat(F,1,nLC);
end
nF = size(Fdt.K,2);
if nF ~= nLC
    Fdt.K = repmat(Fdt.K,1,nLC);
end
nF = size(Fdt.C,2);
if nF ~= nLC
    Fdt.C = repmat(Fdt.C,1,nLC);
end
nF = size(Fdt.M,2);
if nF ~= nLC
    Fdt.M = repmat(Fdt.M,1,nLC);
end

% Initial conditions; acceleration from the equation of motion
x0  = repmat(U0,1,nLC);
xd0 = repmat(V0,1,nLC);
if isdiag(M)
    Minv = 1./diag(M);
    xdd0 = InitialAcceleration(@(r) Minv.*r,F,Fdt,K,C,x0,xd0,...
                               SignalType,ExciteParam);
else
    xdd0 = InitialAcceleration(@(r) M\r,F,Fdt,K,C,x0,xd0,...
                               SignalType,ExciteParam);
end

% Effective (modified) stiffness matrix
K_eff = 1/(beta_newmark*dt^2)*M + gamma_newmark/(beta_newmark*dt)*C + ...
        (1+alpha)*K;

% Parameter matrices
A1 = 1/(dt^2*beta_newmark)*M + gamma_newmark/(dt*beta_newmark)*C + alpha*K;
A2 = 1/(dt*beta_newmark)*M + (gamma_newmark/beta_newmark-1)*C;
A3 = (1/(2*beta_newmark)-1)*M + dt*(gamma_newmark/(2*beta_newmark)-1)*C;

% K_eff is unsymmetric in general: LU decomposition (chosen by MATLAB),
% computed once and reused in every time step
dK_eff = decomposition(K_eff);

% Result variables (one page per load case), starting with the initial
% values
uvec = zeros(NumTimeStep+1,numel(DOFsave),nLC);
vvec = zeros(NumTimeStep+1,numel(DOFsave),nLC);
avec = zeros(NumTimeStep+1,numel(DOFsave),nLC);
uvec(1,:,:) = x0(DOFsave,:);
vvec(1,:,:) = xd0(DOFsave,:);
avec(1,:,:) = xdd0(DOFsave,:);

pb = ProgressBar('HHT-alpha method (implicit)',NumTimeStep);

disp('Solve the equations of motion using the HHT-alpha method - mRHS');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector
T = (0:NumTimeStep)'*dt;

% Pre-compute the force vector
[Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,T);

% Time integration scheme
for i = 2:NumTimeStep+1

    % Compute the force vector at the time t
    Fstar  = LoadAtTime(F,Fdt,Amp(i),dAmp(i),ddAmp(i)) + ...
             A1*x0 + A2*xd0 + A3*xdd0;
    
    % Solve the system of equations at the time t + dt to compute the
    % displacement
    x1 = dK_eff\Fstar;
    
    % Compute the velocity at the time t + dt
    xd1 = gamma_newmark/(beta_newmark*dt)*(x1-x0)   + ...
         (1-(gamma_newmark/beta_newmark))*xd0 	    + ...
         (1-gamma_newmark/(2*beta_newmark))*dt*xdd0;
    
    % Compute the acceleration at the time t + dt
    xdd1 = 1/(beta_newmark*dt^2)*(x1-x0-xd0*dt)-(1/(2*beta_newmark)-1)*xdd0;
    
    % Save variables
    uvec(i,:,:) = x1(DOFsave,:);
    vvec(i,:,:) = xd1(DOFsave,:);
    avec(i,:,:) = xdd1(DOFsave,:);
    
    % Reset variables for the next time step
    x0   = x1;
	xd0  = xd1;
	xdd0 = xdd1;
    
    pb.update(i);
end
pb.close();

% Results (all DOFs if every active DOF was saved; ExpandToAllDoFs only
% handles 2-D histories)
if numel(DOFsave) == numel(ActiveDoFs)
    U = zeros(NumTimeStep+1,NumDoF,nLC);
    V = zeros(NumTimeStep+1,NumDoF,nLC);
    A = zeros(NumTimeStep+1,NumDoF,nLC);
    U(:,ActiveDoFs,:) = uvec;
    V(:,ActiveDoFs,:) = vvec;
    A(:,ActiveDoFs,:) = avec;
else
    U = uvec;
    V = vvec;
    A = avec;
end

end
