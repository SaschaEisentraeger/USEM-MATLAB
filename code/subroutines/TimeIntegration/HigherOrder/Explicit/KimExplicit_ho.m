%% Higher-order explicit time integration: Kim
% 09/06/2019
% UNSW Sydney
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
% pInt:			Order of the time integration scheme
%-------------------------------------------------------------------------%
% Output
% T:            Time vector
% U:            Displacements
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kim, W.: Higher Order Explicit Time Integration Methods for Numerical
% Analyses of Structural Dynamics, Latin American Journal of Solids and
% Structures, 2019, online
%
% In this article, third and fourth order accurate explicit time
% integration methods are developed for effective analyses of various
% linear and nonlinear dynamic problems stated by second order ordinary
% differential equations in time. Two sets of the new methods are
% developed by employing the collocation approach in the time domain. To
% remedy some shortcomings of using the explicit Runge Kutta methods for
% second order ordinary differential equations in time, the new methods are
% designed to introduce small period and damping errors in the important
% low frequency range. For linear cases, the explicitness of the new
% methods is not affected by the presence of non diagonal damping matrix.
% For nonlinear cases, the new methods can handle velocity dependent
% problems explicitly without decreasing order of accuracy. The new methods
% do not have any undetermined algorithmic parameters. Improved numerical
% solutions are obtained when they are applied to various linear and
% nonlinear problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = KimExplicit_ho(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim,pInt)

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

% Order of accuracy of the time-integrator: 3 or 4
IntOrder = pInt;
if IntOrder ~= 3 && IntOrder ~= 4
	IntOrder = 3;
	disp('The order of the time integrator was set to pInt = 3!')
end
% Note that the third-order scheme is fourth-order accurate for the
% undamped case and third-order accurate if damping is present. The fourth-
% order scheme is, however, only fourth-order accurate for both cases.

% Determine the time increment used for the time integration
% Stability limit of the third-order method:  2*pi*0.574976/omega_max -->
% The stability limit is 1.806340 times higher than that of the CDM.
% Stability limit of the fourth-order method: 2*pi*0.474023/omega_max -->
% The stability limit is 1.489187 times higher than that of the CDM.
% dtmax(2) is the stability limit of the CDM: 2/omega_max
if IntOrder == 3
	beta = pi*0.574976;
else
	beta = pi*0.474023;
end
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim,beta);

%% Higher-order time integration algorithm
% Set-up of the parameter of the time integrator
if IntOrder == 3
	tau1    =  1.0/3.0;
	tau2    =  2.0/3.0;
	tau3    =  1.0;
	alpha11 =  1.0;
	alpha21 = -2.0/3.0;
	alpha22 =  2.0/3.0;
	alpha31 =  1.0/3.0;
	alpha32 = -2.0/3.0;
	alpha33 =  1.0/3.0;
else
	tau1    =  1.0/3.0;
	tau2    =  1.0/2.0;
	tau3    =  1.0;
	alpha11 =  1.0;
	alpha21 = -3.0/5.0;
	alpha22 =  3.0/5.0;
	alpha31 =  3.0/5.0;
	alpha32 = -3.0/5.0;
	alpha33 =  0.0;
	alpha41 =  0.0;
	alpha42 =  0.0;
	alpha43 =  0.0;
	alpha44 =  0.0;
end

% Approximation functions (0 < t < dt)
% Stage 1
phi11 = @(t)  1.0;
phi12 = @(t)  t;
% Stage 2
phi21 = @(t)  1.0;
phi22 = @(t) -t*(t-2*tau1*dt)/(2*tau1*dt);
phi23 = @(t)  t^2/(2*tau1*dt);
% Stage 3
phi31 = @(t)  1.0;
phi32 = @(t)  t*(2*t^2-3*tau1*dt*t-3*tau2*dt*t+6*tau1*tau2*dt^2)/...
                (6*tau1*tau2*dt^2);
phi33 = @(t) -t^2*(-2*t+3*tau2*dt)/(6*tau1*(tau1-tau2)*dt^2);
phi34 = @(t)  t^2*(-2*t+3*tau1*dt)/(6*tau2*(tau1-tau2)*dt^2);
if IntOrder == 4
	% Stage 4
	phi41 = @(t)  1.0;
	phi42 = @(t) -t*(3*t^3-4*dt*t^2*(tau1+tau2+tau3)+...
                  6*dt^2*t*(tau1*tau2+tau2*tau3+tau1*tau3))/...
                  (12*tau1*tau2*tau3*dt^3) + t;
	phi43 = @(t)  t^2*(3*t^2-4*tau2*dt*t-4*tau3*dt*t+6*tau2*tau3*dt^2)/...
                  (12*tau1*(-tau1*tau2+tau2*tau3-tau3*tau1+tau1^2)*dt^3);
	phi44 = @(t)  t^2*(3*t^2-4*tau1*dt*t-4*tau3*dt*t+6*tau1*tau3*dt^2)/...
                  (12*tau2*(-tau1*tau2-tau2*tau3+tau3*tau1+tau2^2)*dt^3);
	phi45 = @(t)  t^2*(3*t^2-4*tau1*dt*t-4*tau2*dt*t+6*tau1*tau2*dt^2)/...
                  (12*tau3*( tau1*tau2-tau2*tau3-tau3*tau1+tau3^2)*dt^3);
end

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
Tvec = zeros(NumTimeStep+1,1);
Umat = zeros(NumTimeStep+1,numel(DOFsave));
% Assign the initial conditions to the results
Tvec(1)   = 0;
Umat(1,:) = x0(DOFsave);

pb = ProgressBar('Higher-order Kim method (explicit)',NumTimeStep);

disp(['Solve the equations of motion using the higher-order explicit ' ...
      'time-integrator (Kim)']);
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time integration scheme
for i = 1:NumTimeStep
    % Simulation time
    t = i*dt;
	
	% Stage 1
	xtau1   = phi11(tau1*dt)*x0 + phi12(tau1*dt)*xd0 + ...
              1/2*(tau1*dt)^2*alpha11*xdd0;
	xdtau1  = phi11(tau1*dt)*xd0 + phi12(tau1*dt)*xdd0;
    xddtau1 = KimStageAcceleration(solveM,F,Fdt,K,C,xtau1,xdtau1,...
                                SignalType,ExciteParam,t-(1-tau1)*dt);
	% Stage 2
	xtau2   = phi21(tau2*dt)*x0 + phi22(tau2*dt)*xd0 + ...
              phi23(tau2*dt)*xdtau1 + ...
              1/2*(tau2*dt)^2*(alpha21*xdd0 + alpha22*xddtau1);
	xdtau2  = phi21(tau2*dt)*xd0 + phi22(tau2*dt)*xdd0 + ...
              phi23(tau2*dt)*xddtau1;
    xddtau2 = KimStageAcceleration(solveM,F,Fdt,K,C,xtau2,xdtau2,...
                                SignalType,ExciteParam,t-(1-tau2)*dt);
	% Stage 3
	xtau3   = phi31(tau3*dt)*x0 + phi32(tau3*dt)*xd0 + ...
              phi33(tau3*dt)*xdtau1 + phi34(tau3*dt)*xdtau2 + ...
              1/2*(tau3*dt)^2*(alpha31*xdd0 + alpha32*xddtau1 + ...
              alpha33*xddtau2);
	xdtau3  = phi31(tau3*dt)*xd0 + phi32(tau3*dt)*xdd0 + ...
              phi33(tau3*dt)*xddtau1 + phi34(tau3*dt)*xddtau2;
    xddtau3 = KimStageAcceleration(solveM,F,Fdt,K,C,xtau3,xdtau3,...
                                SignalType,ExciteParam,t-(1-tau3)*dt);
	if IntOrder == 4
		% Stage 4
		x1   = phi41(dt)*x0 + phi42(dt)*xd0 + phi43(dt)*xdtau1 + ...
               phi44(dt)*xdtau2 + phi45(dt)*xdtau3 + ...
               1/2*dt^2*(alpha41*xdd0 + alpha42*xddtau1 + ...
               alpha43*xddtau2 + alpha44*xddtau3);
		xd1  = phi41(dt)*xd0 + phi42(dt)*xdd0 + phi43(dt)*xddtau1 + ...
               phi44(dt)*xddtau2 + phi45(dt)*xddtau3;
        xdd1 = KimStageAcceleration(solveM,F,Fdt,K,C,x1,xd1,...
                                 SignalType,ExciteParam,t);
	else
		x1   = xtau3;
		xd1  = xdtau3;
		xdd1 = xddtau3;
	end
	
    % Save the results
 	Tvec(i+1) = t;
	Umat(i+1,:) = x1(DOFsave)';

    % Update the variables for the next time step
    x0   = x1;
    xd0  = xd1;
	xdd0 = xdd1;

    pb.update(i);
end
pb.close();

% Results (all DOFs if every active DOF was saved)
T = Tvec;
U = ExpandToAllDoFs(NumDoF,ActiveDoFs,Umat);

end