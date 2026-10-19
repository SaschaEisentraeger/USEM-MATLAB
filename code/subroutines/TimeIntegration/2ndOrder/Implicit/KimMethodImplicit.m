%% Implicit time integration: Kim method
% 04.09.2018
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
% Kim, W. & Choi, S. Y.: An Improved Implicit Time Integration Algorithm:
% The Generalized Composite Time Integration Algorithm  Computers &
% Structures, Elsevier BV, 2018, 196, 341-354
%
% The weighted residual method is employed to develop one- and two-step
% time integration schemes. Newly developed time integration schemes are
% combined to obtain a new second-order accurate implicit time integration
% algorithm whose computational structure is similar to the Bathe method
% (Bathe and Noh, 2012). The newly developed algorithm can control
% algorithmic dissipation in the high frequency limit through the optimized
% weighting parameters. It contains only one free parameter, and always
% provides an identical effective stiffness matrix to the first and second
% sub-steps in linear analyses, which is not provided in the algorithm
% proposed by Kim and Reddy (2016). Various nonlinear test problems are
% used to investigate performance of the new algorithm in nonlinear
% analyses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U,V,A] = KimMethodImplicit(K,M,C,F,F_dt,U0,V0,NumDoF,...
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

% Parameter for the time integration scheme
rho_infinity = 0.95;

% Set of weighting parameters
SetWP = 3;

if SetWP == 1
    % First admissible set of weighting parameters
    theta11 = 0.5;
    theta21 = 0.5;
    theta22 = (5+rho_infinity)/(4*(3+rho_infinity));
elseif SetWP == 2
    % Second admissible set of weighting parameters
    rho1_infinity = 1;
    theta11 = 1/(rho1_infinity+1);
    theta21 = 0.5;
    theta22 = (4*rho1_infinity+rho_infinity+1)/...
              (4*(2*rho1_infinity+rho_infinity+1));
elseif SetWP == 3
    % Third admissible set of weighting parameters
    theta11 = (-2+sqrt(2*rho_infinity+2))/(rho_infinity-1);
    theta21 = 0.5;
    theta22 = (rho_infinity-3+sqrt(2*rho_infinity+2))/(4*(rho_infinity-1));
elseif SetWP == 4
    % Non-dissipative case
    theta11 = 0.5;
    theta21 = 0.5;
    theta22 = 0.375;
end

% First sub-step
c1 = 2/(theta11*dt);
c2 = -c1;
c3 = (theta11-1)/theta11;
% Second sub-step
d1 = (1-4*theta21)/((theta21-2*theta22)*dt);
d2 = 4*(2*theta21-1)/((theta21-2*theta22)*dt);
d3 = (3-4*theta21)/((theta21-2*theta22)*dt);
d4 = 4*(theta21-theta22)/(theta21-2*theta22);
d5 = (1-3*theta21+2*theta22)/(theta21-2*theta22);

% Effective stiffness matrices (to increase the efficiency a Cholesky- or
% LU-decomposition of these matrices should be saved --> increases also the
% amount of memory that is required for the simulation)
if SetWP ~= 3
    K_eff1 = c1^2*M + c1*C + K;
    K_eff2 = d1^2*M + d1*C + K;
elseif SetWP == 3
    K_eff1 = c1^2*M + c1*C + K;
end

% Cholesky-decomposition of the effective stiffness matrices
if SetWP ~= 3
	if flag_decomp == 1
		% Cholesky decomposition of the stiffness matrix
		[R1,~,order1] = chol(K_eff1,'vector');
		[R2,~,order2] = chol(K_eff2,'vector');
	elseif flag_decomp == 2
		% Matrix decomposition for solving linear systems
		dK_eff1 = decomposition(K_eff1);
		dK_eff2 = decomposition(K_eff2);
	end
elseif SetWP == 3
	if flag_decomp == 1
		% Cholesky decomposition of the stiffness matrix
		[R1,~,order1] = chol(K_eff1,'vector');
		R2 = R1;
		order2 = order1;
	elseif flag_decomp == 2
		% Matrix decomposition for solving linear systems
		dK_eff1 = decomposition(K_eff1);
		dK_eff2 = dK_eff1;
	end
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
	h = waitbar(0,'Kim method (implicit): ');
end

disp('Solve the equations of motion using the Kim method (implicit)');
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
    % First sub-step
    % Compute the load vector
    f1 = F + FdtN*Amp(i)   + FdtK*Amp(i)   + FdtC*dAmp(i)   + FdtM*ddAmp(i);
    % Effective load vector
    f1_eff = f1 - M*(c1*c2*x0+(c1*c3+c2)*xd0+c3*xdd0) - C*(c2*x0+c3*xd0);
    % Compute the displacement, velocity and acceleration vectors
    % x1   = K_eff1\f1_eff;
	if flag_decomp == 1
        % Cholesky decomposition
        x1 = zeros(size(x0));
		x1(order1) = R1\(R1'\f1_eff(order1));
    elseif flag_decomp == 2
        % Matrix decomposition
        x1 = dK_eff1\f1_eff;
	end
    xd1  = c1*x1 + c2*x0 + c3*xd0;
    xdd1 = c1*xd1 + c2*xd0 + c3*xdd0;
    
    % Second sub-step
    % Compute the load vector
    f2 = F + FdtN*Amp(i+1) + FdtK*Amp(i+1) + FdtC*dAmp(i+1) + FdtM*ddAmp(i+1);
    % Effective load vector
    f2_eff = f2 - M*(d1*d3*x0+d1*d2*x1+(d1*d5+d3)*xd0+ ...
                    (d1*d4+d2)*xd1+d5*xdd0+d4*xdd1)  - ...
                  C*(d3*x0+d2*x1+d5*xd0+d4*xd1);
    % Compute the displacement, velocity and acceleration vectors
    % x2   = K_eff2\f2_eff;
	if flag_decomp == 1
        % Cholesky decomposition
        x2 = zeros(size(x0));
		x2(order2) = R2\(R2'\f2_eff(order2));
    elseif flag_decomp == 2
        % Matrix decomposition
        x2 = dK_eff2\f2_eff;
 	end
    xd2  = d1*x2 + d2*x1 + d3*x0 + d4*xd1 + d5*xd0;
    xdd2 = d1*xd2 + d2*xd1 + d3*xd0 + d4*xdd1 + d5*xdd0;
    
    % Save variables
	uvec(i+1,:) = x2(DOFsave);
    vvec(i+1,:) = xd2(DOFsave);
    avec(i+1,:) = xdd2(DOFsave);
	
	if nUAll > 0
        if mod(i,nSave) == 0
            c = i/nSave;
            UAll{c} = x2;
		elseif i == NumTimeStep
            UAll{nUAll} = x2;
        end
	end
    
    % Reset variables for the next time step
    x0   = x2;
	xd0  = xd2;
	xdd0 = xdd2;
    
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/NumTimeStep,h,['Kim method (implicit): ' ...
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