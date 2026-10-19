%% Implicit time integration: Kim method based on collocation
% 05.09.2018
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
% Kim, W. & Reddy, J. N.: An Improved Time Integration Algorithm: A
% Collocation Time Finite Element Approach, International Journal of
% Structural Stability and Dynamics, World Scientific Pub Co Pte Lt,
% 2017, 17, 1750024
%
% A time collocation finite element approach is employed to develop one-
% and two-step time integration schemes with algorithmic dissipation
% control capability. The newly developed time integration schemes are
% combined to obtain a new family of time integration algorithms using the
% concept employed by Baig and Bathe. The newly developed algorithm can
% effectively control the algorithmic dissipation by relating the
% collocation parameters with the spectral radius in the high frequency
% limit. The new algorithm provides better accuracy compared with the
% generalized-alpha method for highly dissipative cases and includes the
% Baig and Bathe method within its family as a special case.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U,V,A] = KimMethodImplicitCollocation(K,M,C,F,F_dt,U0,V0,...
                     NumDoF,DOFsave,ActiveDoFs,ExciteParam,SignalType,...
                     dtmax,tsim)

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

% User parameters
tau = 0.5;		% 0.5 <= tau  <  1.0
rho1 = 1.0;		% 0.0 <= rho1 <= 1.0
rho2 = 0.9;     % 0.0 <= rho2 <= 1.0  (ultimate spectral radius)

% Derived parameters
theta1 = 1/(1+rho1);	% Second order accuracy is obtained for theta1 = 1.0
theta2 = (tau^2*theta1*(rho2-1)+1)/(2*(1-tau*theta1*(1-rho2))) + ...
         sqrt(tau^4*theta1^2*(rho2-1)^2-4*tau^2*(1-tau)*theta1^2*(rho2-1)+...
		 2*tau^2*theta1*(rho2+1)-4*tau*theta1+1)/(2*(1-tau*theta1*(1-rho2)));
	
% Algorithm coefficients
% First sub-step
c1  =  1/(tau*theta1*dt);
c2  = -c1;
c3  =  (theta1-1)/theta1;
c4  =  c1^2;
c5  =  c1;
c6  =  c1^2;
c7  =  c1*(1-c3);
c8  = -c3;
c9  =  c1;
c10 = -c3;
% Second sub-step
d1  =  (tau-2*theta2)/(theta2*(tau-theta2)*dt);
d2  =  (2*theta2-1)/(tau*theta2*(tau-theta2)*dt);
d3  =  ((1-tau)*(tau+1-2*theta2))/(tau*theta2*(tau-theta2)*dt);
d4  =  (theta2-1)/(tau*(theta2-tau));
d5  =  ((theta2-1)*(tau-1))/(tau*theta2);
d6  =  d1^2;
d7  =  d1;
d8  = -d1*d3;
d9  = -d1*d2;
d10 = -(d1*d5+d3);
d11 = -(d1*d4+d2);
d12 = -d5;
d13 = -d4;
d14 = -d3;
d15 = -d2;
d16 = -d5;
d17 = -d4;

% Effective stiffness matrices (to increase the efficiency a Cholesky- or
% LU-decomposition of these matrices should be saved --> increases also the
% amount of memory that is required for the simulation)
K_eff1 = c4*M + c5*C + K;
K_eff2 = d6*M + d7*C + K; 

% Cholesky-decomposition of the effective stiffness matrices
if flag_decomp == 1
    % Cholesky decomposition of the stiffness matrix
    [R1,~,order1] = chol(K_eff1,'vector');
	[R2,~,order2] = chol(K_eff2,'vector');
elseif flag_decomp == 2
    % Matrix decomposition for solving linear systems
    dK_eff1 = decomposition(K_eff1);
	dK_eff2 = decomposition(K_eff2);
end

% Time vector
T = (0:NumTimeStep)'*dt;

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
	h = waitbar(0,'Kim method (collocation, implicit): ');
end

disp('Solve the equations of motion using the Kim method (collocation, implicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time integration scheme
for i = 1:NumTimeStep
    % Simulation time
    tE = i*dt;          % Time at the end of the current interval
    tS = tE - dt;       % Time at the start of the current interval
    
    % First sub-step
    % Compute the load vector
    f1 = F                                                                    + ...
         Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,tS+tau*dt)       + ...
         Fdt.K*ExcitationSignal(SignalType,ExciteParam,tS+tau*dt)             + ...
         Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,tS+tau*dt,1) + ...
         Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,tS+tau*dt,2) ;
    % Effective load vector
    f1_eff = f1 + M*(c6*x0+c7*xd0+c8*xdd0) + C*(c9*x0+c10*xd0);
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
    f2 = F                                                             + ...
         Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,tE)       + ...
         Fdt.K*ExcitationSignal(SignalType,ExciteParam,tE)             + ...
         Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,tE,1) + ...
         Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,tE,2) ;
    % Effective load vector
    f2_eff = f2 + M*(d8*x0+d9*x1+d10*xd0+d11*xd1+d12*xdd0+d13*xdd1) + ...
	         C*(d14*x0+d15*x1+d16*xd0+d17*xd1);
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
			waitbar(i/NumTimeStep,h,['Kim method (collocation, implicit): ' ...
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