%% Explicit time integration: Soares method
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
% Soares, D.: A Novel Family of Explicit Time Marching Techniques for
% Structural Dynamics and Wave Propagation Models, Computer Methods in
% Applied Mechanics and Engineering, 2016, 311, 838-855
%
% In this work, a novel family of explicit time marching techniques is
% presented, for hyperbolic models. The main features of the proposed
% procedure are: (i) the family is based only on single-step displacement-
% velocity relations, being very simple to implement; (ii) as usual in
% explicit analyses, it requires no system of equations to be dealt with,
% being very efficient; (iii) it is truly self-starting; (iv) it is second-
% order accurate; (v) it enables controllable algorithmic dissipation in
% the higher modes; (vi) it provides reduced period elongation errors in
% the lower modes; (vii) it is conditionally stable and exhibits high
% stability limits, with its critical sampling frequency varying between
% 3.571 and 4 (i.e., its upper limit is twice that of the Central
% Difference Method). The proposed technique is very effective, providing
% the main positive attributes that are requested from an efficacious
% explicit time marching procedure. Numerical results are presented along
% the paper, illustrating the excellent performance of the novel method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = SoaresMethodExplicit(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

% (De-)activate the waitbar
flag_wb = true;

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

% Parameter for the time integration scheme
% Pre-selected parameters that optimize the performance of the algorithm
b1 = 8.567e-3;
b2 = 8.590e-1;
% The beta parameter control the algorithm's dissipative features. No
% numerical damping of high frequency (spurious) oscillations is introduced
% for beta = 0, while the maximal numerical dissipation is provided for
% beta = 1 (default choice in the current implementation).
beta = 1;
% For the current set-up an excellent performance of the algorithm for the
% relative temporal discretization dt/T < 0.1 is observed. In this case,
% almost negligible period elongation takes place and spurious high
% frequency modes are effectively damped. Nonetheless, the low order modes
% are captured accurately by this second order accurate time integration
% scheme.

% Determine the time increment used for the time integration
% The stability limit of the Soares explicit time integration method is
% higher than that of the central difference method (CDM) (which is used to
% determine the critical time increment).
% If no numerical damping is introduced (beta = 0) the stability limit is
% twice that of the CDM.
% In the case of maximal numerical dissipation the stability limit is still
% 1.7855 times that of the CDM.
% Stability limit dt = 1.7855*dt_CDM
dtmax(2) = 1.7855*dtmax(2);
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

% Parameter for the interpolation of the force vector (trapezoidal rule)
beta1 = 0.5;
beta2 = 0.5;
% For beta1 = 1 and beta2 = 0 an explicit algorithm is also obtained for
% the load term.

% Effective (modified) mass matrix
E = M + 1/2*dt*C;
% For an explicit algorithm it is assumed that both the mass and damping
% matrices are diagonal and therefore the inversion of the effective mass
% matrix becomes trivial.
if isdiag(E)
    flag_diag = true;
    Einv = 1./diag(E);
else
    flag_diag = false;
	% Matrix decomposition for solving linear systems
    dE = decomposition(E);
end

% Result variables
uvec = zeros(NumTimeStep+1,NumDoFSave);
% Assign the initial conditions to the results
uvec(1,:) = x0(DOFsave);

if flag_wb
	% Waitbar
	h = waitbar(0,'Soares Method (explicit): ');
end

disp('Solve the equations of motion using the Soares method');
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

    if flag_diag
        % Solve the system of equations at the time t + dt to compute the
        % velocity
        xd1 = Einv.*(Fstar + M*xd0 - 0.5*dt*C*xd0 - K*(dt*x0+0.5*dt^2*xd0));
        % Compute the displacement at the time t + dt
        x1 = x0 + 0.5*dt*(xd0 + xd1) - Einv.*(0.5*dt^2*C*xd1 + ...
             K*(beta*b1*b2*dt^3*xd0 + (1/16+beta*b1)*dt^3*xd1));
    else
        % Solve the system of equations at the time t + dt to compute the
        % velocity
        xd1 = dE\(Fstar + M*xd0 - 0.5*dt*C*xd0 - K*(dt*x0+0.5*dt^2*xd0));
        % Compute the displacement at the time t + dt
        x1 = x0 + 0.5*dt*(xd0 + xd1) - dE\(0.5*dt^2*C*xd1 + ...
             K*(beta*b1*b2*dt^3*xd0 + (1/16+beta*b1)*dt^3*xd1));
    end
    
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
			waitbar(i/NumTimeStep,h,['Soares method (explicit): ' ...
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