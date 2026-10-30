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
% flag_wb:      Show a waitbar
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

% Time step: stability limit of the method itself (factor 1)
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim,1);

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

% Effective matrix (dynamic stiffness matrix)
Keff = 1/dt^2*M + 1/(2*dt)*C;
% Diagonal mass and damping matrices make the scheme fully explicit
% (no linear system per time step)
if isdiag(M)
    flag_diag_M = true;
    Minv = 1./diag(M);
	if isdiag(C)
		flag_diag_C = true;
		Keffinv = 1./diag(Keff);
	else
		flag_diag_C = false;
		dKeff = decomposition(Keff);
	end
else
    flag_diag_M = false;
	flag_diag_C = false;
	dKeff = decomposition(Keff);
    dM = decomposition(M);
end

% Initial conditions; acceleration from the equation of motion
x0  = U0;
xd0 = V0;
if flag_diag_M
    xdd0 = InitialAcceleration(@(r) Minv.*r,F,Fdt,K,C,x0,xd0,...
                               SignalType,ExciteParam);
else
    xdd0 = InitialAcceleration(@(r) dM\r,F,Fdt,K,C,x0,xd0,...
                               SignalType,ExciteParam);
end

% Start up computation
xm1 = x0 - dt*xd0 + dt^2/2*xdd0;

% Result variables
uvec = zeros(NumTimeStep+1,numel(DOFsave));
vvec = zeros(NumTimeStep+1,numel(DOFsave));
avec = zeros(NumTimeStep+1,numel(DOFsave));
% Assign the initial conditions (displacements, velocities) and the initial
% acceleration to the results
uvec(1,:) = x0(DOFsave);
vvec(1,:) = xd0(DOFsave);
avec(1,:) = xdd0(DOFsave);

pb = ProgressBar('Central difference method (explicit)',NumTimeStep,flag_wb);

disp('Solve the equations of motion using the CDM (Central Difference Method)');
disp('');
disp('-----------------------------------------------------------------------');
disp('');

% Time vector and excitation signal
T = (0:NumTimeStep)'*dt;
[Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,T);

% Auxiliary matrices
A = 1/dt^2*M - 1/(2*dt)*C;
B = K - 2/dt^2*M;

% Time integration scheme
for i = 1:NumTimeStep
	% Force vector
	Fhat = LoadAtTime(F,Fdt,Amp(i),dAmp(i),ddAmp(i));
	Ft = Fhat - A*xm1 - B*x0;
	
	% Displacement
	if flag_diag_M && flag_diag_C
		x1 = Keffinv.*Ft;
	else
		x1 = dKeff\Ft;
	end
	
	% Velocity at t_n (central difference)
	xd1 = 1/(2*dt)*(x1 - xm1);

	% Acceleration at t_n
	xdd1 = 1/dt^2*(x1 - 2*x0 + xm1);

    % Save variables: displacement at t_{n+1}, velocity and acceleration
    % at t_n (the initial values at t_0 are kept)
    uvec(i+1,:) = x1(DOFsave);
    if i > 1
        vvec(i,:) = xd1(DOFsave);
        avec(i,:) = xdd1(DOFsave);
    end
    
    % Reset variables for the next time step
    % Displacement
	xm1 = x0;
    x0 = x1;

    pb.update(i);
end
pb.close();

% Velocity and acceleration at the final time (linear extrapolation, the
% central differences would require the displacement at t_{N+1})
if NumTimeStep > 1
    vvec(end,:) = 2*vvec(end-1,:) - vvec(end-2,:);
    avec(end,:) = 2*avec(end-1,:) - avec(end-2,:);
else
    vvec(end,:) = vvec(1,:);
    avec(end,:) = avec(1,:);
end

% Results (all DOFs if every active DOF was saved)
[U,V,A] = ExpandToAllDoFs(NumDoF,ActiveDoFs,uvec,vvec,avec);

end