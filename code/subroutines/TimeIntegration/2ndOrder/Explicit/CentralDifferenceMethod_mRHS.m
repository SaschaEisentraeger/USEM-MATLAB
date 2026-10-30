%% Explicit time integration: Central difference method (CDM)
%__________________________________________________________________________
% Numerical solution of the equation of motion by means of the CDM with
% multiple right-hand sides (mRHS), i.e., multiple load cases are
% simultaneously considered in the simulation.
%__________________________________________________________________________
%
% 20/04/2022
% UNSW Sydney
% Sascha Eisentraeger
%-------------------------------------------------------------------------%
% Input
% K:            Stiffness matrix
% M:            Mass matrix
% C:            Damping matrix
% F:            Force vector/can be a matrix if multiple
%               load cases are considered
% F_dt          Force vector (time-dependent)/can be a matrix if multiple
%               load cases are considered
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
function [T,U,V,A] = CentralDifferenceMethod_mRHS(K,M,C,F,F_dt,U0,V0,...
                     NumDoF,DOFsave,ActiveDoFs,ExciteParam,SignalType,...
                     dtmax,tsim,flag_wb)

% Time step: stability limit of the method itself (factor 1)
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim,1);

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

% Check for multiple load cases
nLC = size(Fdt.Neumann,2);

% Adjust the size of the force vectors
nF = size(F,2);
if nF ~= nLC
    F = repmat(F,1,nLC);
end
nF = size(Fdt.K,2);
if nF ~= nLC
    Fdt.K = repmat(Fdt.K,1,nLC);
end
nF = size(Fdt.C,2);
if nF ~= nLC
    Fdt.C = repmat(Fdt.C,1,nLC);
end
nF = size(Fdt.M,2);
if nF ~= nLC
    Fdt.M = repmat(Fdt.M,1,nLC);
end

% Effective matrix (dynamic stiffness matrix)
Keff = 1/dt^2*M + 1/(2*dt)*C;

% Matrix decomposition for solving linear systems
dKeff = decomposition(Keff);
dM = decomposition(M);

% Initial conditions (one column per load case); acceleration from the
% equation of motion
x0  = repmat(U0,1,nLC);
xd0 = repmat(V0,1,nLC);
xdd0 = InitialAcceleration(@(r) dM\r,F,Fdt,K,C,x0,xd0,...
                           SignalType,ExciteParam);

% Start up computation
xm1 = x0 - dt*xd0 + dt^2/2*xdd0;

% Result variables
uvec = zeros(NumTimeStep+1,numel(DOFsave),nLC);
vvec = zeros(NumTimeStep+1,numel(DOFsave),nLC);
avec = zeros(NumTimeStep+1,numel(DOFsave),nLC);
% Assign the initial conditions (displacements, velocities) and the initial
% acceleration to the results
uvec(1,:,:) = x0(DOFsave,:);
vvec(1,:,:) = xd0(DOFsave,:);
avec(1,:,:) = xdd0(DOFsave,:);

pb = ProgressBar('Central difference method (explicit)',NumTimeStep,flag_wb);

disp('Solve the equations of motion using the CDM (Central Difference Method) - mRHS');
disp('');
disp('------------------------------------------------------------------------------');
disp('');

% Time vector and excitation signal
T = (0:NumTimeStep)'*dt;
[Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,T);

% Auxilliary matrices
A = 1/dt^2*M - 1/(2*dt)*C;
B = K - 2/dt^2*M;

% Time integration scheme
for i = 1:NumTimeStep
	% Force vector
	Fhat = LoadAtTime(F,Fdt,Amp(i),dAmp(i),ddAmp(i));
	Ft = Fhat - A*xm1 - B*x0;
	
	% Displacement
    x1 = dKeff\Ft;
	
	% Velocity at t_n (central difference)
	xd1 = 1/(2*dt)*(x1 - xm1);

	% Acceleration at t_n
	xdd1 = 1/dt^2*(x1 - 2*x0 + xm1);

    % Save variables: displacement at t_{n+1}, velocity and acceleration
    % at t_n (the initial values at t_0 are kept)
    uvec(i+1,:,:) = x1(DOFsave,:);
    if i > 1
        vvec(i,:,:) = xd1(DOFsave,:);
        avec(i,:,:) = xdd1(DOFsave,:);
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
    vvec(end,:,:) = 2*vvec(end-1,:,:) - vvec(end-2,:,:);
    avec(end,:,:) = 2*avec(end-1,:,:) - avec(end-2,:,:);
else
    vvec(end,:,:) = vvec(1,:,:);
    avec(end,:,:) = avec(1,:,:);
end

% Results (all DOFs if every active DOF was saved; ExpandToAllDoFs handles
% 2-D histories only)
if numel(ActiveDoFs) == numel(DOFsave)
    U = zeros(NumTimeStep+1,NumDoF,nLC);
    V = zeros(NumTimeStep+1,NumDoF,nLC);
    A = zeros(NumTimeStep+1,NumDoF,nLC);
    U(:,ActiveDoFs,:) = uvec;
    V(:,ActiveDoFs,:) = vvec;
    A(:,ActiveDoFs,:) = avec;
else
    U = uvec;
    V = vvec;
    A = avec;
end

end