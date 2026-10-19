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
function [T,U,V] = SoaresMethodImplicit_ho(K,M,C,F,Fdt,U0,V0,NumDoF,...
                   DOFsave,ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

% (De-)activate the waitbar
flag_wb = true;

% Decomposition method: (1) Cholesky decomposition; (2) Matrix
% decomposition (Matlab decides on an appropriate method)
flag_decomp = 2;

% Time integration parameter:
% alpha = 1 (no numerical damping/dissipation is introduced)
% alpha > 1 (numerical damping/dissipation is introduced)
alpha = 1;

% Critical sampling frequency
Omega_c = sqrt(2/alpha*(sqrt(1+12*alpha)-1));

% Determine the time increment used for the time integration
dtmax(2) = dtmax(2)/2*Omega_c;
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

% Initial conditions
% Displacement
x0  = U0;
% Velocity
xd0 = V0;

% Effective (modified) stiffness matrix
K_eff = M + 1.0/2.0*dt*C + 1.0/6.0*dt^2*K;

if flag_decomp == 1
    % Cholesky decomposition of the stiffness matrix
    [R,~,order] = chol(K_eff,'vector');
elseif flag_decomp == 2
    % Matrix decomposition for solving linear systems
    dK_eff = decomposition(K_eff);
end

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
uvec = zeros(NumTimeStep+1,NumDoFSave);
vvec = zeros(NumTimeStep+1,NumDoFSave);
% Assign the initial conditions (displacements, velocities) and the initial
% acceleration to the results
uvec(1,:) = x0(DOFsave);
vvec(1,:) = xd0(DOFsave);

if flag_wb
	% Waitbar
	h = waitbar(0,'Soares Method (implicit): ');
end

disp('Solve the equations of motion using the high-order Soares method (implicit)');
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
for i = 1:NumTimeStep
    % Compute the effective force vector at the time t + 1/2*dt
    F0 = F + FdtN*Amp(i)   + FdtK*Amp(i)   + FdtC*dAmp(i)   + FdtM*ddAmp(i);
    F1 = F + FdtN*Amp(i+1) + FdtK*Amp(i+1) + FdtC*dAmp(i+1) + FdtM*ddAmp(i+1);
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
	if flag_decomp == 1
        % Cholesky decomposition
        xd1 = zeros(size(x0));
		xd1(order) = R\(R'\R1(order));
    elseif flag_decomp == 2
        % Matrix decomposition
        xd1 = dK_eff\R1;
	end
	
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
    
	if nUAll > 0
        if mod(i,nSave) == 0
            c = i/nSave;
            UAll{c} = x1;
		elseif i == NumTimeStep
            UAll{nUAll} = x1;
        end
	end
	
    % Reset variables for the next time step
	S0 = S1;
	V0 = V1;
    x0  = x1;
	xd0 = xd1;
	
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/NumTimeStep,h,['High-order Soares method (implicit): ' ...
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
    U(:,ActiveDoFs) = uvec;
    V(:,ActiveDoFs) = vvec;
else
    U = uvec;
    V = vvec;
end

if nUAll > 0
    save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
              filesep 'UAll.mat'],'-v7.3')
end

end