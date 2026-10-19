%% Implicit time integration: Soares method
% 12.04.2017
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
% Soares, D.: A Simple and Effective New Family of Time Marching Procedures
% for Dynamics, Computer Methods in Applied Mechanics and Engineering,
% 2015, 283, 1138-1166
%
% In this work, a new family of time-marching procedures for dynamics is
% introduced. The proposed algorithm is based on displacements-velocities
% relations, requiring no computation of accelerations. Thus, not only it
% is a very simple technique, but also it is truly self-starting.
% Unconditionally stable implicit and conditionally stable explicit
% approaches are focused here, and optimized procedures are discussed for
% these analyses. The new technique is second-order accurate with
% controllable period elongation and algorithmic dissipation. Numerical
% results are presented, illustrating the good performance of the method.
% The new family is simple, efficient and highly accurate, resulting in a
% very effective time-marching technique.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = SoaresMethodImplicit(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

% (De-)activate the waitbar
flag_wb = true;

% Decomposition method: (1) Cholesky decomposition; (2) Matrix
% decomposition (Matlab decides on an appropriate method)
flag_decomp = 2;

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

% Determine the time increment used for the time integration
N_dt = numel(unique(dtmax));
if N_dt == 3
	% Deliberately described time step
	dt = dtmax(3);
else
	dt = dtmax(1);
end

% Parameter set
ParamSet = 1;

if ParamSet == 1
    % Parameter for the time integration scheme
    d = 1/2;
    a = 0.01; % controls the damping
    omega = sqrt(min(diag(K)./diag(M)));
    delta1 = 1/2 + d*omega*dt;
    alpha1 = 1/2 + 3/2*tanh(a*omega*dt);
    delta2 = 1 - delta1;
    alpha2 = 2*sqrt(2*alpha1) - alpha1 - 2*delta1;
%-------------------------------------------------------------------------%
elseif ParamSet == 2
    % We can also obtain an explicit algorithm for the following choices of
    % the integration parameter
    d = -1/4;
    omega = sqrt(min(diag(K)./diag(M)));
    delta1 = 1/2 + d*omega*dt;
    alpha1 = 0;
    delta2 = 1 - delta1;
    alpha2 = 2*delta2;
%-------------------------------------------------------------------------%
elseif ParamSet == 3
    % Test the high-order scheme by Soares
    delta1 = 1.0/2.0;
    alpha1 = 1.0/3.0;
    delta2 = 1.0/2.0;
    alpha2 = 2.0/3.0;

    % delta1 >= 1.0/2.0 and alpha1 < delta1 the method is only
    % conditionally stable and the critical time increment has to be
    % determined
    if delta1 >= 0.5 && alpha1 < delta1
        % Damping ratio
        xi = 0;
        % Maximum (angular) eigenfrequency --> computed from the critical
        % time step size of the centra difference method:
        % dt_CDM = 2/omega_max
        omega_max = 2/dtmax(4);
        % Critical time step size
        dtc = ((delta1 - 0.5)*xi + sqrt(0.5*(delta1-alpha1) + ...
               (delta1-0.5)^2*xi^2))/(0.5*(delta1-alpha1)*omega_max);
    else
        dtc = dtmax(1);
    end

    % Determine the time increment used for the time integration
    dt = min([dtc, dtmax(1:3)]);
end
%-------------------------------------------------------------------------%

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

% Parameter for the interpolation of the force vector (trapezoidal rule)
beta1 = 0.5;
beta2 = 0.5;
% For beta1 = 1 and beta2 = 0 an explicit algorithm is also obtained for
% the load term.

% Effective (modified) stiffness matrix
K_eff = M + delta1*dt*C + 1/2*alpha1*dt^2*K;

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
% Assign the initial conditions to the results
uvec(1,:) = x0(DOFsave);

if flag_wb
	% Waitbar
	h = waitbar(0,'Soares Method (implicit): ');
end

disp('Solve the equations of motion using the Soares method (implicit)');
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
    Fstar = beta1*dt*F0 + beta2*dt*F1;
% 	Fstar = 0.5*dt*(F + Fdt*ExcitationSignal(SignalType,ExciteParam,t+0.5*dt));
    % Compute the effective load vector at the time t + dt
    F_eff = Fstar + M*xd0 - delta2*dt*C*xd0 - ...
            K*(dt*x0 + 1/2*alpha2*dt^2*xd0);
    % Solve the system of equations at the time t + dt to compute the
    % velocity
	if flag_decomp == 1
        % Cholesky decomposition
        xd1 = zeros(size(x0));
		xd1(order) = R\(R'\F_eff(order));
    elseif flag_decomp == 2
        % Matrix decomposition
        xd1 = dK_eff\F_eff;
	end
    % Compute the displacement at the time t + dt
    x1 = x0 + delta2*dt*xd0 + delta1*dt*xd1;
    
    % Save variables
    uvec(i+1,:) = x1(DOFsave);
    
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
    
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/NumTimeStep,h,['Soares method (implicit): ' ...
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
    U(:,ActiveDoFs) = uvec;
else
    U = uvec;
end

if nUAll > 0
    save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
              filesep 'UAll.mat'],'-v7.3')
end

end