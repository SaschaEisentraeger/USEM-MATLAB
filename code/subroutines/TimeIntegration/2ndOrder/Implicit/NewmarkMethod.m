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
function [T,U,V,A] = NewmarkMethod(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                     ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

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
beta_newmark = 1./4;
gamma_newmark = 1./2;

% Effective (modified) stiffness matrix
K_eff = K + 1/(beta_newmark*dt^2)*M + gamma_newmark/(beta_newmark*dt)*C;

% Parameter matrices
A1 = 1/(dt^2*beta_newmark)*M + gamma_newmark/(dt*beta_newmark)*C;
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
    h = waitbar(0,'Newmark-Method (implicit): ');
end

disp('Solve the equations of motion using the Newmark method');
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
        x1(order) = R\(R'\Fstar(order));
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
    uvec(i,:) = x1(DOFsave);
    vvec(i,:) = xd1(DOFsave);
    avec(i,:) = xdd1(DOFsave);
    
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
            waitbar(i/NumTimeStep,h,['Newmark method (implicit): ' ...
            num2str(i/NumTimeStep*100,'%5.1f') '%'])
        end
    end
end

% save('Contour.mat','UAll');

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