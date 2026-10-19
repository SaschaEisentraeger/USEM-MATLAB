%% Implicit time integration: beta1-beta2-Bathe method
% 15/03/2019
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
% Malakiyeh, M. M.; Shojaee, S. & Bathe, K.-J.: The Bathe Time Integration
% Method Revisited for Prescribing Desired Numerical Dissipation,
% Computers & Structures, 2019, 212, 289-298
%
% In this paper we further consider the Bathe method for the direct time
% integration in structural dynamics and wave propagations. The method uses
% two sub-steps per time step and is an unconditionally stable scheme
% frequently used without adjusting any parameter. In the first sub-step
% the trapezoidal rule is used and in the second sub-step the 3-point Euler
% backward method is employed. In this contribution we derive the method
% using, instead of the Euler scheme, the 3-point trapezoidal rule for the
% complete step with two Newmark parameters. The parameters can then be
% used to smoothly prescribe desired numerical dissipation, from zero to
% very significant dissipation. To highlight the performance of the method,
% the stability, accuracy and overshooting are studied and some
% illustrative problems are solved. The results are compared with those of
% some other methods that also use parameters to introduce numerical
% dissipation. We conclude that the use of the parameters in the Bathe
% method can be valuable but probably will require some numerical
% experimentation.
%
% Malakiyeh, M. M.; Shojaee, S. & Bathe, K.-J.: New insights into the
% b1/b2-Bathe time integration scheme when L-stable, Computers &
% Structures, 2021, 245, 106433
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U,V,A] = b1b2_BatheMethod(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
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

% Parameter for the time integration scheme (influence dispersion and
% dissipation); unconditional stability is achieved for
% (1) beta1 \in [1/3,1/2] and beta2 = 1-beta1 (with increasing beta1
% amplitude decay and period elongation decrease beta1=beta2=0.5 no decay;
% always second order accurate) --> A-stable scheme (i.e., rho <=0)
% or (2) beta1 \in [1/3,\infty] and beta2 = 2*beta1 (with increasing beta1
% amplitude decay and period elongation increase,; second order accurate
% only for beta1 = 1/3) --> L-stable scheme (i.e., rho <=0 and rho = 0 for
% dt/T --> infinity)
p = 0.5;
beta1 = 1.0/3.0;	% has to be >= 1/3!
beta2 = 1 - beta1;
% Using the following values for beta1 and beta2 the beta1-beta2-Bathe
% method is identical to the standard Bathe method independent of the
% choice of p
% p  = 2 - sqrt(2);
% beta1 = 1 + 1/(2*p*(p-2));
% beta2 = 1/(2-p);
% Using the following value for p results in the identity of Keff1 and
% Keff2 (only one factorization is required); independent of the choice of
% beta1 and beta2
% beta1 = 1.0/2.0;
% beta2 = 1.0/2.0;
% p = 2*beta2/(1+2*beta2);
% Derived parameters
a1 = 2/(p*dt);
a2 = 4/(p^2*dt^2);
a3 = 4/(p*dt);
b1 = beta2*(1-p)*dt;
b2 = b1*b1;
b3 = p*(1-beta1);
b4 = beta2*(1-p);
b5 = b4*b4*dt;
b6 = p*beta1;
b7 = (1-beta2)*(1-p);
c1 = 1/b1;
c2 = 1/b2;
c3 = (b3+b4)/b5;
c4 = (b6+b7)/b5;
c5 = b3/b4;
c6 = (b6+b7)/b4;

% Effective stiffness matrices (to increase the efficiency a Cholesky- or
% LU-decomposition of these matrices should be saved --> increases also the
% amount of memory that is required for the simulation)
K_eff1 = a2*M + a1*C + K;
K_eff2 = c2*M + c1*C + K;

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
	h = waitbar(0,'beta1-beta2-Bathe method (implicit): ');
end

disp('Solve the equations of motion using the beta1-beta2-Bathe-method (implicit)');
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
    F_eff1 = F_t_pdt + M*(a2*x0 + a3*xd0 + xdd0) + C*(a1*x0 + xd0);
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
    xdp = (xp - x0)*a1 - xd0;
	% Compute the acceleration at the time t + p*dt
    xddp = (xdp - xd0)*a1 - xdd0;
    
    % Second part of the algorithm
    % Compute the effective load vector at the time t + dt
    F_eff2 = F_t_dt + M*(c2*x0 + c3*xd0 + c4*xdp + c5*xdd0 + c6*xddp) + ...
	                  C*(c1*x0 + c5*xd0 + c6*xdp);
	
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
    xd1 = c1*(x1-x0) - c5*xd0 - c6*xdp;

    % Compute the acceleration at the time t + dt
    xdd1 = c2*(x1-x0) - c3*xd0 - c4*xdp - c5*xdd0 - c6*xddp;
    
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
			waitbar(i/NumTimeStep,h,['beta1-beta2-Bathe method (implicit): ' ...
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