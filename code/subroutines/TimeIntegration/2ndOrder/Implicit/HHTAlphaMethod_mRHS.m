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

% (De-)activate the waitbar
flag_wb = true;

% Parameter for the time integration scheme: -1/3 <= alpha <= 0
% If alpha = 0 --> Newmark method!
% If alpha = -0.05 --> Default value in Abaqus/Standard!
alpha = -0.05;
beta_newmark = (1-alpha)^2/4;
gamma_newmark = 1./2-alpha;

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

% Determine if the mass matrix is in diagonal form
if isdiag(M)
    flag_diag = true;
    Minv = 1./diag(M);
else
    flag_diag = false;
end

% Initial conditions
% Displacement
x0  = repmat(U0,1,nLC);
% Velocity
xd0 = repmat(V0,1,nLC);

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

% Effective (modified) stiffness matrix
K_eff = 1/(beta_newmark*dt^2)*M + gamma_newmark/(beta_newmark*dt)*C + ...
        (1+alpha)*K;

% Parameter matrices
A1 = 1/(dt^2*beta_newmark)*M + gamma_newmark/(dt*beta_newmark)*C + alpha*K;
A2 = 1/(dt*beta_newmark)*M + (gamma_newmark/beta_newmark-1)*C;
A3 = (1/(2*beta_newmark)-1)*M + dt*(gamma_newmark/(2*beta_newmark)-1)*C;

if flag_decomp == 1
    % Cholesky decomposition of the stiffness matrix
    [R,~,order] = chol(K_eff,'vector');
elseif flag_decomp == 2
    % Matrix decomposition for solving linear systems
    dK_eff = decomposition(K_eff);
end

% Result variables
uvec = zeros(NumTimeStep+1,NumDoFSave,nLC);
vvec = zeros(NumTimeStep+1,NumDoFSave,nLC);
avec = zeros(NumTimeStep+1,NumDoFSave,nLC);
% Assign the initial conditions (displacements, velocities) and the initial
% acceleration to the results
uvec(1,:,:) = x0(DOFsave,:);
vvec(1,:,:) = xd0(DOFsave,:);
avec(1,:,:) = xdd0(DOFsave,:);

if flag_wb
	% Waitbar 
	h = waitbar(0,'HHT-alpha method (implicit): ');
end

disp('Solve the equations of motion using the HHT-alpha method - mRHS');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector
T = (0:NumTimeStep)'*dt;

% Pre-compute the force vector
Amp   = ExcitationSignal(SignalType,ExciteParam,T);
dAmp  = ExcitationSignalDerivative(SignalType,ExciteParam,T,1);
ddAmp = ExcitationSignalDerivative(SignalType,ExciteParam,T,2);
FdtN = Fdt.Neumann;
FdtK = Fdt.K;
FdtC = Fdt.C;
FdtM = Fdt.M;
clear Fdt

% Time integration scheme
for i = 2:NumTimeStep+1

    % Compute the force vector at the time t
    Fstar  = F + FdtN*Amp(i) + FdtK*Amp(i) + FdtC*dAmp(i) + ...
             FdtM*ddAmp(i) + A1*x0 + A2*xd0 + A3*xdd0;
    
    % In each time step a system of equations has to be solved K*x = F.
    % Since the inversion of K will result in a fully populated matrix a
    % Cholesky decomposition is performed instead. The solution is computed
    % by means of forward and backward substitution.
    % Solve the system of equations at the time t + dt to compute the
    % displacement
    if flag_decomp == 1
        % Cholesky decomposition
        x1 = zeros(size(x0));
        x1(order,:) = R\(R'\Fstar(order,:));
    elseif flag_decomp == 2
        % Matrix decomposition
        x1 = dK_eff\Fstar;
    end
    
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
    
	if nUAll > 0
        if mod(i,nSave) == 0
            c = i/nSave;
            UAll{c} = x1;
		elseif i == NumTimeStep+1
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
			waitbar(i/NumTimeStep,h,['HHT-alpha method (implicit): ' ...
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

if nUAll > 0
    save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
              filesep 'UAll.mat'],'-v7.3')
end

end