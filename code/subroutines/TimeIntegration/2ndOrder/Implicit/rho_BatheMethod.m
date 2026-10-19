%% Implicit time integration: rho-Bathe method
% 04/07/2019
% University of Magdeburg
% Sascha Duczek
%-------------------------------------------------------------------------%
% Input
% K:            Stiffness matrix
% M:            Mass matrix
% C:            Damping matrix
% F:            Force vector
% F_dt:     	Force vector (time-dependent)
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
% Noh, G. & Bathe, K.-J.: The Bathe Time Integration Method with 
% Controllable Spectral Radius: The rho_infinity-Bathe Method, Computers & 
% Structures, Elsevier BV, 2019, 212, 299-310
%
% We consider the Bathe implicit time integration method and focus on the
% time step splitting ratio and the spectral radius at large time steps to
% improve and generalize the scheme. The objective is to be able to
% prescribe the amplitude decay (dissipation) and period elongation
% (dispersion) for the numerical integration, and to achieve this aim in a
% direct and optimum manner with the minimum number of parameters. We show
% that the use of the time step splitting ratio and spectral radius is
% effective to prescribe in a smooth manner no amplitude decay to very
% large amplitude decays, with correspondingly small period elongation to
% very large period elongations while maintaining second-order accuracy. We
% analyze the effects of the splitting ratio and spectral radius on the
% stability and accuracy of the scheme and illustrate the use of these
% parameters in comparison with previously published methods. Furthermore,
% we show that with a proper setting of these parameters more accurate
% results may be obtained in some analyses.
%
% Noh, G. & Bathe, K.-J.: For Direct Time Integrations: A Comparison of the
% Newmark and rho_infinity-Bathe Schemes, Computers & Structures,
% Elsevier BV, 2019, 225, 106079
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U,V,A] = rho_BatheMethod(K,M,C,F,Fdt,U0,V0,NumDoF,DOFsave,...
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

if flag_wb
	% Less updates of the waitbar (avoid slowing the simulation down)
	nUpdate = ceil(NumTimeStep/1000);
end

% Adjusted time increment
dt = tsim/NumTimeStep;

% Save the full displacement field for a few time steps
nUAll = 0;
if nUAll > 0
    nSave = ceil(NumTimeStep/nUAll);
    UAll = cell(nUAll,1);
end

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
% This method includes as special cases also the standard Bathe method and
% the beta1-beta2-Bathe method.
% Increasing rho results in less amplitude decay and period elongation.
rho = 0.0;
% For large numerical dissipation p > 1.3 (rho = 0)
if rho == 1.0
	gamma = 1/2;
else
	gamma = @(x) (2-sqrt(2+2*x))/(1-x);
end
p = gamma(rho);
q1 = (rho+1)/(2*p*(rho-1)+4);
q0 = (p-1)*q1 + 0.5;
q2 = -p*q1 + 0.5;
a0 = 4/(p^2*dt^2);
a1 = 2/(p*dt);
a2 = 1/(dt^2*q2^2);
a3 = 1/(dt*q2);
a4 = 4/(p*dt);
a5 = (q0+q2)/(dt*q2^2);
a6 = q1/(dt*q2^2);
a7 = q0/q2;
a8 = q1/q2;

% Effective stiffness matrices (to increase the efficiency a Cholesky- or
% LU-decomposition of these matrices should be saved --> increases also the
% amount of memory that is required for the simulation)
K_eff1 = K + a1*C + a0*M;
K_eff2 = K + a3*C + a2*M;

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
	h = waitbar(0,'rho-Bathe method (implicit): ');
end

disp('Solve the equations of motion using the rho-Bathe method (implicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector
T  = (0:NumTimeStep)'*dt;
% Steps in between
Tp = (0:NumTimeStep)'*dt - (1-p)*dt;

% Pre-compute the amplitude values for the force vectors
Amp   = ExcitationSignal(SignalType,ExciteParam,T);
dAmp  = ExcitationSignalDerivative(SignalType,ExciteParam,T,1);
ddAmp = ExcitationSignalDerivative(SignalType,ExciteParam,T,2);
Amp_p   = ExcitationSignal(SignalType,ExciteParam,Tp);
dAmp_p  = ExcitationSignalDerivative(SignalType,ExciteParam,Tp,1);
ddAmp_p = ExcitationSignalDerivative(SignalType,ExciteParam,Tp,2);

% Parts of the force vector
FdtN = Fdt.Neumann;
FdtK = Fdt.K;
FdtC = Fdt.C;
FdtM = Fdt.M;

% Time integration scheme
for i = 1:NumTimeStep
	% First part of the algorithm
    % Compute the load vector at the time t + p*dt
    F_t_pdt = F + FdtN*Amp_p(i+1) + FdtK*Amp_p(i+1) + FdtC*dAmp_p(i+1) + FdtM*ddAmp_p(i+1);
    % Compute the load vector at the time t + dt
    F_t_dt  = F + FdtN*Amp(i+1) + FdtK*Amp(i+1) + FdtC*dAmp(i+1) + FdtM*ddAmp(i+1);
% 	% Compute the load vector at the time t
% 	F_t = F + FdtN*Amp(i) + FdtK*Amp(i) + FdtC*dAmp(i) + FdtM*ddAmp(i);
% 	% Linear interpolation of the force vector between time steps
% 	F_t_pdt = (1-p)*F_t + p*F_t_dt;
    % Compute the effective load vector at the time t + p*dt
    F_eff1 = F_t_pdt + M*(a0*x0 + a4*xd0 + xdd0) + C*(a1*x0 + xd0);
 	% Solve the system of equations at the time t + p*dt to compute the
    % displacement: xp = K_eff1\F_eff1;
	if flag_decomp == 1
        % Cholesky decomposition
        xp = zeros(size(x0));
		xp(order1) = R1\(R1'\F_eff1(order1));
    elseif flag_decomp == 2
        % Matrix decomposition
        xp = dK_eff1\F_eff1;
	end
    % Compute the velocity at the time t + p*dt
    xdp  = a1*(xp - x0) - xd0;
	% Compute the acceleration at the time t + p*dt
    xddp = a1*(xdp - xd0) - xdd0;
    
    % Second part of the algorithm
    % Compute the effective load vector at the time t + dt
    F_eff2 = F_t_dt + ...
	         M*(a2*x0 + a5*xd0 + a6*xdp + a7*xdd0 + a8*xddp) + ...
			 C*(a3*x0 + a7*xd0 + a8*xdp);
    % Solve the system of equations at the time t + dt to compute the
    % displacement: x1 = K_eff2\F_eff2;
	if flag_decomp == 1
        % Cholesky decomposition
        x1 = zeros(size(x0));
		x1(order2) = R2\(R2'\F_eff2(order2));
    elseif flag_decomp == 2
        % Matrix decomposition
        x1 = dK_eff2\F_eff2;
	end
    % Compute the velocity at the time t + dt
    xd1  = a3*(x1-x0) - a7*xd0  - a8*xdp;
    % Compute the acceleration at the time t + dt
    xdd1 = a3*(xd1-xd0) - a7*xdd0 - a8*xddp;
    
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
			waitbar(i/(NumTimeStep+1),h,['rho-Bathe method (implicit): ' ...
			num2str(i/(NumTimeStep+1)*100,'%5.1f') '%'])
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