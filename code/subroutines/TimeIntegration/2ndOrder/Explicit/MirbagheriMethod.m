%% Explicit time integration: Mirbagheri method
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
% Mirbagheri, Y.; Nahvi, H.; Parvizian, J. & Düster, A.: Reducing Spurious
% Oscillations in Discontinuous Wave Propagation Simulation Using High-
% Order Finite Elements, Computers & Mathematics with Applications, 2015,
% 70, 1640-1658
%
% In this paper, a modified explicit time integration scheme is proposed
% for simulating the propagation of discontinuous waves. The spatial domain
% is discretized by the finite element method. To obtain accurate results,
% the standard finite element method requires a very fine mesh, which is
% why the computational effort can be very time consuming. The use of high-
% order finite element methods – such as the spectral element method based
% on Lagrange polynomials through Gauss–Lobatto–Legendre points or the
% isogeometric analysis using non-uniform rational B-splines – can reduce
% the enormous computational costs significantly, compared to the standard
% finite element method. However, explicit time integration schemes such as
% the central difference method cannot eliminate the spurious oscillation
% near the front wave. The procedure proposed here is basically similar to
% the explicit method of Noh and Bathe, but the semi-discrete equation of
% motion is modified by introducing a damping parameter to suit the high-
% order FEM analysis of the discontinuous wave propagation. The performance
% due to this modification is tested for one-dimensional wave propagation
% problems using both lumped and consistent mass matrices. Then, the idea
% of a combination of the consistent and row sum lumped mass matrices is
% evaluated. The proposed method is studied also in two-dimensions by
% considering the Lamb problem of wave propagation. The results are
% promising enough to provide a better simulation of the discontinuous wave
% propagation problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = MirbagheriMethod(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,....
                 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

% (De-)activate the waitbar
flag_wb = true;

% The Mirbagheri method is a straightforward extension of the Noh-Bathe
% method where additional damping is introduced during the computation of
% the displacements. Therefore, in each sub-step the displacement vector is
% corrected. The difference between the internal forces in the Mirbagheri
% and in the central difference method for a free undamped system is used
% as additional damping force (Fdamp = K*U_M - K*U_CDM)

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

% Inverted mass matrix
if isdiag(M)
    flag_diag = true;
    Minv = 1./diag(M);
else
    flag_diag = false;
	% Matrix decomposition for solving linear systems
    dM = decomposition(M);
end

% Initial conditions
% Displacement
x0  = U0;
% Velocity
xd0 = V0;

% Acceleration
if flag_diag
	xdd0 = Minv.*(F - C*xd0 - K*x0                                              + ...
				  Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,0)        + ...
				  Fdt.K*ExcitationSignal(SignalType,ExciteParam,0)              + ...
				  Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,0,1)  + ...
				  Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,0,2));
else
	xdd0 = dM\(F - C*xd0 - K*x0                                             + ...
              Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,0)        + ...
              Fdt.K*ExcitationSignal(SignalType,ExciteParam,0)              + ...
			  Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,0,1)  + ...
			  Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,0,2));
end

% Parameter for the time integration scheme
p  = 0.54; % for p = 1 the Noh-Bathe method is identical to the CDM
q1 = (1-2*p)/(2*p*(1-p));
q2 = 1/2-p*q1;
q0 = -q1-q2+1/2;

% Stability limit of the Noh-Bathe method compared to the CDM
% (dt ~= 1.8725*dt_CDM)
gamma = 1/4 - 1/2*(1-p)*q1;
val = 1/2*sqrt(1/(gamma*p*(1-p)));
dtmax(2) = val*dtmax(2);
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

% Parameter for the time integration scheme
a0 = p*dt;
a1 = 1/2*(p*dt)^2;
a2 = a0/2;
a3 = (1-p)*dt;
a4 = 1/2*((1-p)*dt)^2;
a5 = q0*a3;
a6 = (1/2+q1)*a3;
a7 = q2*a3;
% Damping parameter (if eta is set to 1.0 this method is identical to the
% Noh-Bathe method)/recommendation: 1.0 <= eta <= 1.1; --> higher values of
% eta result in a higher numerical damping (dissipation)
eta = 1.025;

% Result variables
uvec = zeros(NumTimeStep+1,NumDoFSave);
% Assign the initial conditions to the results
uvec(1,:) = x0(DOFsave);

if flag_wb
	% Waitbar
	h = waitbar(0,'Mirbagheri method (explicit): ');
end

disp('Solve the equations of motion using the Mirbagheri method (explicit)');
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
    % Compute the displacement at the time t + p*dt
    xp = x0 + a0*xd0 + a1*eta*xdd0;
	% Compute the load vector at the time t + p*dt
    F_t_pdt = F + FdtN*Amp_p(i+1) + FdtK*Amp_p(i+1) + FdtC*dAmp_p(i+1) + FdtM*ddAmp_p(i+1);
    % Compute the load vector at the time t + dt
    F_t_dt  = F + FdtN*Amp(i+1) + FdtK*Amp(i+1) + FdtC*dAmp(i+1) + FdtM*ddAmp(i+1);
%     % Compute the load vector at the time t
%     F_t = F + FdtN*Amp(i) + FdtK*Amp(i) + FdtC*dAmp(i) + FdtM*ddAmp(i);
%     % Linear interpolation of the force vector between time steps
% 	F_t_pdt = (1-p)*F_t + p*F_t_dt;
    % Compute the effective load vector at the time t + p*dt
    F_t_pdt_circ    = F_t_pdt - K*xp - C*(xd0 + a0*xdd0);
    % Solve the system of equations at the time t + p*dt to compute the
    % acceleration
    if flag_diag
        xddp = Minv.*F_t_pdt_circ;
    else
        xddp = dM\F_t_pdt_circ;
    end
    % Compute the velocity at the time t + p*dt
    xdp = xd0 + a2*(xdd0 + xddp);
    
    % Second part of the algorithm
    % Compute the displacement at the time t + dt
    x1 = xp + a3*xdp + a4*eta*xddp;
    % Compute the load vector at the time t + dt
    F_t_dt_circ = F_t_dt - K*x1 - C*(xdp + a3*xddp);
 	% Solve the system of equations at the time t + dt to compute the
    % acceleration
    if flag_diag
        xdd1 = Minv.*F_t_dt_circ;
    else
        xdd1 = dM\F_t_dt_circ;
    end
    % Compute the velocity at the time t + dt
    xd1 = xdp + a5*xdd0 + a6*xddp + a7*xdd1;
    
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
	xdd0 = xdd1;
    
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/NumTimeStep,h,['Mirbagheri method (explicit): ' ...
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