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

% (De-)activate the waitbar
flag_wb = true;

% Decomposition method: (1) Cholesky decomposition; (2) Matrix
% decomposition (Matlab decides on an appropriate method)
flag_decomp = 2;

% Determine the time increment used for the time integration
N_dt = numel(unique(dtmax));
if N_dt == 3
	% Deliberately described time step
	dt = dtmax(3);
else
	dt = dtmax(1);
end

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

% Determine if the mass matrix is in diagonal form
if isdiag(M)
    flag_diag = true;
    Minv = 1./diag(M);
else
    flag_diag = false;
end

% Initial conditions
% Displacement
x0  = U0;
% Velocity
xd0 = V0;

% Acceleration
if flag_diag
	xdd0 = Minv.*(F - C*xd0 - K*x0 + ...
				  Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,0)        + ...
				  Fdt.K*ExcitationSignal(SignalType,ExciteParam,0)              + ...
				  Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,0,1)  + ...
				  Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,0,2));
else
	xdd0 = M\(F - C*xd0 - K*x0 + ...
              Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,0)        + ...
              Fdt.K*ExcitationSignal(SignalType,ExciteParam,0)              + ...
			  Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,0,1)  + ...
			  Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,0,2));
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
    
% Cholesky-decomposition of the stiffness matrix
if flag_decomp == 1
    % Cholesky decomposition of the stiffness matrix
    [R,~,order] = chol(K_eff,'vector');
elseif flag_decomp == 2
    % Matrix decomposition for solving linear systems
    dK_eff = decomposition(K_eff);
end

% Result variables
uvec = zeros(NumTimeStep+1,NumDoFSave);
vvec = zeros(NumTimeStep+1,NumDoFSave);
avec = zeros(NumTimeStep+1,NumDoFSave);
% Assign the initial conditions (displacements, velocities) and the initial
% acceleration to the results
uvec(1,:) = x0(DOFsave);
vvec(1,:) = xd0(DOFsave);
avec(1,:) = xdd0(DOFsave);

if flag_wb
	% Waitbar
	h = waitbar(0,'Generalized-alpha method (implicit): ');
end

disp('Solve the equations of motion using the generalized-alpha method (implicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector
T = (0:NumTimeStep)'*dt;
% Steps in between (time vector)
Talpha = (0:NumTimeStep)'*dt - alpha_f*dt;

% Pre-compute the amplitude values for the force vectors
Amp   = ExcitationSignal(SignalType,ExciteParam,Talpha);
dAmp  = ExcitationSignalDerivative(SignalType,ExciteParam,Talpha,1);
ddAmp = ExcitationSignalDerivative(SignalType,ExciteParam,Talpha,2);

% Parts of the force vector
FdtN = Fdt.Neumann;
FdtK = Fdt.K;
FdtC = Fdt.C;
FdtM = Fdt.M;

% Time integration scheme
for i = 1:NumTimeStep
	% The force vector is computed according to the original publication by
	% Chung & Hulbert
    f = F + FdtN*Amp(i+1) + FdtK*Amp(i+1) + FdtC*dAmp(i+1) + FdtM*ddAmp(i+1);
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
        
    % In each time step a system of equations has to be solved K*x = F.
    % Since the inversion of K will result in a fully populated matrix a
    % Cholesky-decomposition is performed instead. The solution is computed
    % by means of forward and backward substitution.
    % Solve the system of equations at the time t + dt to compute the
    % displacement
	if flag_decomp == 1
        % Cholesky decomposition
        x1 = zeros(size(x0));
        x1(order) = R\(R'\Fstar(order));
    elseif flag_decomp == 2
        % Matrix decomposition
        x1 = dK_eff\Fstar;
	end
    
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
	
	if nUAll > 0
        if mod(i,nSave) == 0
            c = i/nSave;
            UAll{c} = x1;
		elseif i == NumTimeStep
            UAll{nUAll} = x1;
        end
	end
    
    % Reset variables for the next time step
    x0   = x1;
	xd0  = xd1;
	xdd0 = xdd1;
    
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/NumTimeStep,h,['Generalized-alpha method (implicit): ' ...
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
    V = zeros(NumTimeStep+1,NumDoF);
    A = zeros(NumTimeStep+1,NumDoF);
    U(:,ActiveDoFs) = uvec;
    V(:,ActiveDoFs) = vvec;
    A(:,ActiveDoFs) = avec;
else
    U = uvec;
    V = vvec;
    A = avec;
end

if nUAll > 0
    save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
              filesep 'UAll.mat'],'-v7.3')
end

end