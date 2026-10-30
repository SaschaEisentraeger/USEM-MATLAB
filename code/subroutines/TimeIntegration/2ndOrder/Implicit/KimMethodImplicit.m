%% Implicit time integration: Kim method
% 04.09.2018
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
% Kim, W. & Choi, S. Y.: An Improved Implicit Time Integration Algorithm:
% The Generalized Composite Time Integration Algorithm  Computers &
% Structures, Elsevier BV, 2018, 196, 341-354
%
% The weighted residual method is employed to develop one- and two-step
% time integration schemes. Newly developed time integration schemes are
% combined to obtain a new second-order accurate implicit time integration
% algorithm whose computational structure is similar to the Bathe method
% (Bathe and Noh, 2012). The newly developed algorithm can control
% algorithmic dissipation in the high frequency limit through the optimized
% weighting parameters. It contains only one free parameter, and always
% provides an identical effective stiffness matrix to the first and second
% sub-steps in linear analyses, which is not provided in the algorithm
% proposed by Kim and Reddy (2016). Various nonlinear test problems are
% used to investigate performance of the new algorithm in nonlinear
% analyses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U,V,A] = KimMethodImplicit(K,M,C,F,F_dt,U0,V0,NumDoF,...
                     DOFsave,ActiveDoFs,ExciteParam,SignalType,dtmax,tsim)

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

% Parameter for the time integration scheme
rho_infinity = 0.95;

% Set of weighting parameters
SetWP = 3;

if SetWP == 1
    % First admissible set of weighting parameters
    theta11 = 0.5;
    theta21 = 0.5;
    theta22 = (5+rho_infinity)/(4*(3+rho_infinity));
elseif SetWP == 2
    % Second admissible set of weighting parameters
    rho1_infinity = 1;
    theta11 = 1/(rho1_infinity+1);
    theta21 = 0.5;
    theta22 = (4*rho1_infinity+rho_infinity+1)/...
              (4*(2*rho1_infinity+rho_infinity+1));
elseif SetWP == 3
    % Third admissible set of weighting parameters
    theta11 = (-2+sqrt(2*rho_infinity+2))/(rho_infinity-1);
    theta21 = 0.5;
    theta22 = (rho_infinity-3+sqrt(2*rho_infinity+2))/(4*(rho_infinity-1));
elseif SetWP == 4
    % Non-dissipative case
    theta11 = 0.5;
    theta21 = 0.5;
    theta22 = 0.375;
end

% First sub-step
c1 = 2/(theta11*dt);
c2 = -c1;
c3 = (theta11-1)/theta11;
% Second sub-step
d1 = (1-4*theta21)/((theta21-2*theta22)*dt);
d2 = 4*(2*theta21-1)/((theta21-2*theta22)*dt);
d3 = (3-4*theta21)/((theta21-2*theta22)*dt);
d4 = 4*(theta21-theta22)/(theta21-2*theta22);
d5 = (1-3*theta21+2*theta22)/(theta21-2*theta22);

% Effective stiffness matrices
if SetWP ~= 3
    K_eff1 = c1^2*M + c1*C + K;
    K_eff2 = d1^2*M + d1*C + K;
elseif SetWP == 3
    K_eff1 = c1^2*M + c1*C + K;
end

% The effective matrices are unsymmetric in general: LU decompositions
% (chosen by MATLAB), computed once and reused in every time step
if SetWP ~= 3
    dK_eff1 = decomposition(K_eff1);
    dK_eff2 = decomposition(K_eff2);
elseif SetWP == 3
    % Both sub-steps share the same effective matrix
    dK_eff1 = decomposition(K_eff1);
    dK_eff2 = dK_eff1;
end

% Result variables, starting with the initial values
uvec = zeros(NumTimeStep+1,numel(DOFsave));
vvec = zeros(NumTimeStep+1,numel(DOFsave));
avec = zeros(NumTimeStep+1,numel(DOFsave));
uvec(1,:) = x0(DOFsave);
vvec(1,:) = xd0(DOFsave);
avec(1,:) = xdd0(DOFsave);

pb = ProgressBar('Kim method (implicit)',NumTimeStep);

disp('Solve the equations of motion using the Kim method (implicit)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector
T = (0:NumTimeStep)'*dt;

% Pre-compute the force vector
[Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,T);

% Time integration scheme
for i = 1:NumTimeStep
    % First sub-step
    % Compute the load vector
    f1 = LoadAtTime(F,Fdt,Amp(i),dAmp(i),ddAmp(i));
    % Effective load vector
    f1_eff = f1 - M*(c1*c2*x0+(c1*c3+c2)*xd0+c3*xdd0) - C*(c2*x0+c3*xd0);
    % Compute the displacement, velocity and acceleration vectors
    x1   = dK_eff1\f1_eff;
    xd1  = c1*x1 + c2*x0 + c3*xd0;
    xdd1 = c1*xd1 + c2*xd0 + c3*xdd0;
    
    % Second sub-step
    % Compute the load vector
    f2 = LoadAtTime(F,Fdt,Amp(i+1),dAmp(i+1),ddAmp(i+1));
    % Effective load vector
    f2_eff = f2 - M*(d1*d3*x0+d1*d2*x1+(d1*d5+d3)*xd0+ ...
                    (d1*d4+d2)*xd1+d5*xdd0+d4*xdd1)  - ...
                  C*(d3*x0+d2*x1+d5*xd0+d4*xd1);
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
