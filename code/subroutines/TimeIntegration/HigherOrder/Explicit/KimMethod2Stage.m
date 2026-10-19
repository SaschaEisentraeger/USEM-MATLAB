%% Two-Stage Explicit Time Integration Scheme
% Article:
% Kim, W.: An Accurate Two-Stage Explicit Time Integration Scheme for
% Structural Dynamics and Various Dynamic Problems, International Journal
% for Numerical Methods in Engineering, Wiley, 2019, 120, 1-28
% In this article, a two-stage explicit time integration method is
% presented for the numerical analysis of various dynamic problems
% described by second-order ordinary differential equations in time. The
% newly proposed method have negligible damping and period errors in the
% important low-frequency range and can introduce moderate numerical
% damping into the spurious high-frequency range. The computational effort
% per step required in the new method is the same as the effort needed in
% the existing two-stage explicit methods. The proposed method can provide
% fourth-order accurate solutions for velocity independent problems and
% second-order accurate solutions for velocity dependent problems, which is
% not possible in the existing methods. The proposed new method is expected
% to be an excellent tool for numerical analyses of various dynamic
% problems including structural dynamics.
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
% tsim:         Simulation time
%-------------------------------------------------------------------------%
% Output
% U:            Displacement matrix (the displacement is save for each
%               degree of freedom in each time step)
% T:            Time vector
%-------------------------------------------------------------------------%
function [T,U] = KimMethod2Stage(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

% (De-)activate the waitbar
flag_wb = true;

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

% Adjust the time increment to the simulation time
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

% Splitting parameter (value gives best accuracy)
p = 1/3;

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
	xdd0 = Minv.*(F - C*xd0 - K*x0                                             + ...
				  Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,0)       + ...
				  Fdt.K*ExcitationSignal(SignalType,ExciteParam,0)             + ...
				  Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,0,1) + ...
				  Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,0,2));	
else
	xdd0 = dM\(F - C*xd0 - K*x0                                            + ...
 			  Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,0)       + ...
			  Fdt.K*ExcitationSignal(SignalType,ExciteParam,0)             + ...
			  Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,0,1) + ...
			  Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,0,2)) ;
end		

% Result variables
tvec = zeros(NumTimeStep+1,1);
uvec = zeros(NumTimeStep+1,NumDoFSave);
% Assign the initial conditions to the results
tvec(1) = 0;
uvec(1,:) = x0(DOFsave);

if flag_wb
	% Waitbar
	h = waitbar(0,'Two-Stage Kim Method-Method (explicit): ');
end

disp('Solve the equations of motion using the Two-Stage Kim Method-Method');
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
    % Simulation time
    tE = i*dt;          % Time at the end of the current interval
    % Linear interpolation of the force vector between time steps
	F_t    = F + FdtN*Amp(i)   + FdtK*Amp(i)   + FdtC*dAmp(i)   + FdtM*ddAmp(i);
    F_t_dt = F + FdtN*Amp(i+1) + FdtK*Amp(i+1) + FdtC*dAmp(i+1) + FdtM*ddAmp(i+1);
	F_t_pdt = (1-p)*F_t + p*F_t_dt;
	
	% First stage
    x13   = x0+p*dt*xd0+(p*dt)^2*xdd0/2;
    xd13  = xd0+p*dt*xdd0;
	if flag_diag
		xdd13 = Minv.*(F_t_pdt - K*x13 - C*xd13);
	else
		xdd13 = dM\(F_t_pdt - K*x13 - C*xd13);
	end
	
    % Second stage
    xm1  = x0+dt*xd0+(dt^2)*((3*p-1)*xdd0/(3*p)+xdd13/(3*p))/2;
    xmd1 = xd0+dt*((2*p-1)*xdd0/(2*p)-xdd13/(2*p));
	if flag_diag
		xdd1 = Minv.*(F_t_dt - K*xm1 - C*xmd1);
	else
		xdd1 = dM\(F_t_dt - K*xm1 - C*xmd1);
	end
	
    % Correction 
    x1  = x0+dt*xd0+(1/2)*(dt^2)*((4*p-1)*xdd0/(6*p)-xdd13/(6*p*(p-1))+(2*p-1)*xdd1/(6*(p-1)));
    xd1 = xd0+dt*((3*p-1)*xdd0/(6*p)-xdd13/(6*p*(p-1))+(3*p-2)*xdd1/(6*(p-1))); 
	
    % Save variables
    tvec(i+1)   = tE;
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
			waitbar(i/NumTimeStep,h,['Two-Stage Kim Method-Method (explicit): ' ...
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