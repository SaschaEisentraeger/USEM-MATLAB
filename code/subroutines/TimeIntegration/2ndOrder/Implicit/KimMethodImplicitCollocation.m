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

% Time step (prescribed implicit step) and number of steps
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim);

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,numel(ActiveDoFs));

% Initial conditions; acceleration from the equation of motion
x0  = U0;
xd0 = V0;
if isdiag(M)
    Minv = 1./diag(M);
    xdd0 = InitialAcceleration(@(r) Minv.*r,F,Fdt,K,C,x0,xd0,...
                               SignalType,ExciteParam);
else
    xdd0 = InitialAcceleration(@(r) M\r,F,Fdt,K,C,x0,xd0,...
                               SignalType,ExciteParam);
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

% Effective stiffness matrices
K_eff1 = c4*M + c5*C + K;
K_eff2 = d6*M + d7*C + K; 

% The effective matrices are unsymmetric in general: LU decompositions
% (chosen by MATLAB), computed once and reused in every time step
dK_eff1 = decomposition(K_eff1);
dK_eff2 = decomposition(K_eff2);

% Time vector
T = (0:NumTimeStep)'*dt;

% Result variables, starting with the initial values
uvec = zeros(NumTimeStep+1,numel(DOFsave));
vvec = zeros(NumTimeStep+1,numel(DOFsave));
avec = zeros(NumTimeStep+1,numel(DOFsave));
uvec(1,:) = x0(DOFsave);
vvec(1,:) = xd0(DOFsave);
avec(1,:) = xdd0(DOFsave);

pb = ProgressBar('Kim method (collocation, implicit)',NumTimeStep);

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
    % Compute the load vector at the collocation point t + tau*dt
    [s,ds,dds] = SignalHistory(SignalType,ExciteParam,tS+tau*dt);
    f1 = LoadAtTime(F,Fdt,s,ds,dds);
    % Effective load vector
    f1_eff = f1 + M*(c6*x0+c7*xd0+c8*xdd0) + C*(c9*x0+c10*xd0);
    % Compute the displacement, velocity and acceleration vectors
    x1   = dK_eff1\f1_eff;
    xd1  = c1*x1 + c2*x0 + c3*xd0;
    xdd1 = c1*xd1 + c2*xd0 + c3*xdd0;
    
    % Second sub-step
    % Compute the load vector at the end of the interval
    [s,ds,dds] = SignalHistory(SignalType,ExciteParam,tE);
    f2 = LoadAtTime(F,Fdt,s,ds,dds);
    % Effective load vector
    f2_eff = f2 + M*(d8*x0+d9*x1+d10*xd0+d11*xd1+d12*xdd0+d13*xdd1) + ...
	         C*(d14*x0+d15*x1+d16*xd0+d17*xd1);
    % Compute the displacement, velocity and acceleration vectors
    x2   = dK_eff2\f2_eff;
    xd2  = d1*x2 + d2*x1 + d3*x0 + d4*xd1 + d5*xd0;
    xdd2 = d1*xd2 + d2*xd1 + d3*xd0 + d4*xdd1 + d5*xdd0;
    
    % Save variables
	uvec(i+1,:) = x2(DOFsave);
    vvec(i+1,:) = xd2(DOFsave);
    avec(i+1,:) = xdd2(DOFsave);
    
    % Reset variables for the next time step
    x0   = x2;
	xd0  = xd2;
	xdd0 = xdd2;
    
    pb.update(i);
end
pb.close();

% Results (all DOFs if every active DOF was saved)
[U,V,A] = ExpandToAllDoFs(NumDoF,ActiveDoFs,uvec,vvec,avec);

end
