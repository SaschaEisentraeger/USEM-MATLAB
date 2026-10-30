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

% Time step: stability limit of the central difference method (factor 1)
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim,1);

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

% Splitting parameter (value gives best accuracy)
p = 1/3;

% Solution with the mass matrix (the only linear system of the method;
% element-wise division for a diagonal mass matrix)
if isdiag(M)
    Minv = 1./diag(M);
    solveM = @(r) Minv.*r;
else
    dM = decomposition(M);
    solveM = @(r) dM\r;
end

% Initial conditions; acceleration from the equation of motion
x0   = U0;
xd0  = V0;
xdd0 = InitialAcceleration(solveM,F,Fdt,K,C,x0,xd0,SignalType,ExciteParam);

% Result variables
uvec = zeros(NumTimeStep+1,numel(DOFsave));
% Assign the initial conditions to the results
uvec(1,:) = x0(DOFsave);

pb = ProgressBar('Two-stage Kim method (explicit)',NumTimeStep);

disp('Solve the equations of motion using the Two-Stage Kim Method-Method');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector and excitation signal
T = (0:NumTimeStep)'*dt;
[Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,T);

% Time integration scheme
for i = 1:NumTimeStep
    % Linear interpolation of the force vector between time steps
	F_t     = LoadAtTime(F,Fdt,Amp(i),dAmp(i),ddAmp(i));
    F_t_dt  = LoadAtTime(F,Fdt,Amp(i+1),dAmp(i+1),ddAmp(i+1));
	F_t_pdt = (1-p)*F_t + p*F_t_dt;

	% First stage
    x13   = x0+p*dt*xd0+(p*dt)^2*xdd0/2;
    xd13  = xd0+p*dt*xdd0;
	xdd13 = solveM(F_t_pdt - K*x13 - C*xd13);

    % Second stage
    xm1  = x0+dt*xd0+(dt^2)*((3*p-1)*xdd0/(3*p)+xdd13/(3*p))/2;
    xmd1 = xd0+dt*((2*p-1)*xdd0/(2*p)-xdd13/(2*p));
	xdd1 = solveM(F_t_dt - K*xm1 - C*xmd1);

    % Correction 
    x1  = x0+dt*xd0+(1/2)*(dt^2)*((4*p-1)*xdd0/(6*p)-xdd13/(6*p*(p-1))+(2*p-1)*xdd1/(6*(p-1)));
    xd1 = xd0+dt*((3*p-1)*xdd0/(6*p)-xdd13/(6*p*(p-1))+(3*p-2)*xdd1/(6*(p-1))); 
	
    % Save variables
    uvec(i+1,:) = x1(DOFsave);

    % Reset variables for the next time step
    x0   = x1;
	xd0  = xd1;
	xdd0 = xdd1;

    pb.update(i);
end
pb.close();

% Results (all DOFs if every active DOF was saved)
U = ExpandToAllDoFs(NumDoF,ActiveDoFs,uvec);

end