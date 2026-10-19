%% Explicit time integration: Central difference method (CDM)
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
function [T,U,V,A] = CentralDifferenceMethod(K,M,C,F,F_dt,U0,V0,NumDoF,...
                     DOFsave,ActiveDoFs,ExciteParam,SignalType,dtmax,...
                     tsim,flag_wb)

% Determine the time increment used for the time integration
dt = min(dtmax);

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

% Effective matrix (dynamic stiffness matrix)
Keff = 1/dt^2*M + 1/(2*dt)*C;
% Check whether mass and damping matrices are diagonal
if isdiag(M)
    flag_diag_M = true;
    Minv = 1./diag(M);
	if isdiag(C)
		flag_diag_C = true;
		Keffinv = 1./diag(Keff);
	else
		flag_diag_C = false;
		% Matrix decomposition for solving linear systems
		dKeff = decomposition(Keff);
	end
else
    flag_diag_M = false;
	flag_diag_C = false;
    % Matrix decomposition for solving linear systems
	dKeff = decomposition(Keff);
    dM = decomposition(M);
end

% Initial conditions
% Displacement
x0  = U0;
% Velocity
xd0 = V0;

% Acceleration
if flag_diag_M
	xdd0 = Minv.*(F  - C*xd0 - K*x0                                             + ...
				  Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,0)        + ...
				  Fdt.K*ExcitationSignal(SignalType,ExciteParam,0)              + ...
                  Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,0,1)  + ...
				  Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,0,2));
else
	xdd0 = dM\(F  - C*xd0 - K*x0                                            + ...
			  Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,0)        + ...
			  Fdt.K*ExcitationSignal(SignalType,ExciteParam,0)              + ...
              Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,0,1)  + ...
			  Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,0,2));
end

% Start up computation
xm1 = x0 - dt*xd0 + dt^2/2*xdd0;

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
	h = waitbar(0,'Central Difference Method (explicit): ');
end

disp('Solve the equations of motion using the CDM (Central Difference Method)');
disp('');
disp('-----------------------------------------------------------------------');
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

% Auxiliary matrices
A = 1/dt^2*M - 1/(2*dt)*C;
B = K - 2/dt^2*M;

% Time integration scheme
for i = 1:NumTimeStep
	% Force vector
	Fhat = F + FdtN*Amp(i) + FdtK*Amp(i) + FdtC*dAmp(i) + FdtM*ddAmp(i);
	Ft = Fhat - A*xm1 - B*x0;
	
	% Displacement
	if flag_diag_M && flag_diag_C
		x1 = Keffinv.*Ft;
	else
		x1 = dKeff\Ft;
	end
	
	% Velocity
	xd1 = 1/(2*dt)*(x1 - xm1);
	
	% Acceleration
	xdd1 = 1/dt^2*(x1 - 2*x0 + xm1);
    
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
    % Displacement
	xm1 = x0;
    x0 = x1;
    
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/NumTimeStep,h,['Central Difference Method '...
			'(explicit): ' num2str(i/NumTimeStep*100,'%5.1f') '%'])
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