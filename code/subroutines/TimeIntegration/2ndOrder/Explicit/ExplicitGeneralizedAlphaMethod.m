%% Explicit time integration: Explicit Generalized-Alpha method (EG-Alpha method)
% 13.07.2020, Bolong Tao, UNSW Sydney
% 
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
% Hulbert, G. M. and Chung, J.: Explicit Time Integration Algorithms for
% Structural Dynamics with Optimal Numerical Dissipation, Computer Methods
% in Applied Mechanics and Engineering, 1996, 137, 175-188
%
% A new predictor-corrector explicit time integration algorithm is
% presented for solving structural dynamics problems. The basis of the
% algorithm is the implicit generalized-a method developed by the authors.
% Like its implicit parent, the explicit generalized-alpha method is a one-
% parameter family of algorithms in which the parameter defines the high-
% frequency numerical dissipation. For a given value of high-frequency
% dissipation, the explicit generalized-alpha method minimizes low-
% frequency dissipation. The algorithm can be utilized effectively for
% structural dynamics calculations in which numerical dissipation is needed
% to reduce spurious oscillations. Parameter values are given that enable
% physical damping to be treated explicitly while maintaining second-order
% accuracy without the need for an expensive second corrector pass.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = ExplicitGeneralizedAlphaMethod(K,M,C,F,F_dt,U0,V0,...
				 NumDoF,DOFsave,ActiveDoFs,ExciteParam,SignalType,dtmax,...
				 tsim)

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

% Parameter for the time integration scheme:
% one-parameter family
rho_b = 1;
alpha_f = 1;
alpha_m = (2*rho_b-1)/(1+rho_b);
beta = (5-3*rho_b)/((1+rho_b)^2*(2-rho_b));
gamma = 1.5-alpha_m;

% recover CD method
% rho_p = 1;
% rho_s = 0;
% alpha_f = 0;
% alpha_m = (2*rho_p*rho_s+rho_p-1)/((rho_p+1)*(rho_s+1));
% beta = 0;
% gamma = 1/2 - alpha_m + alpha_f;

% beta = alpha_f + (1-rho_p-2*rho_p*rho_s)/((rho_p+1)*(rho_s+1)) + ...
       %((1-rho_p)*(1-rho_p*rho_s)^2/(alpha_f*(1+alpha_p)^2*(1+alpha_s)*(2+alpha_s-alpha_p*alpha_s));

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
	h = waitbar(0,'Explicit generalized alpha method: ');
end

disp('Solve the equations of motion using the EG-alpha (Explicit Generalized Alpha Method)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time integration scheme
for i = 1:NumTimeStep
    % Simulation time
    tE = i*dt;          % Time at the end of the current interval

    % Predictor
    x1_P = x0 + dt*xd0 + (0.5 - beta)*dt^2*xdd0;
    xd1_P = xd0 + (1-gamma)*dt*xdd0;
    
    % sub time step
    x1malpha_f_P = (1 - alpha_f)*x1_P + alpha_f*x0;
    xd1malpha_f_P = (1 - alpha_f)*xd1_P + alpha_f*xd0;
	if flag_diag
		xdd1malpha_m = Minv.*(F - C*xd1malpha_f_P - K*x1malpha_f_P                                            + ...
				              Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,dt*((i-1)+1-alpha_f))       + ...
				              Fdt.K*ExcitationSignal(SignalType,ExciteParam,dt*((i-1)+1-alpha_f))             + ...
				              Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,dt*((i-1)+1-alpha_f),1) + ...
				              Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,dt*((i-1)+1-alpha_f),2));
	else
		xdd1malpha_m = dM\(F - C*xd1malpha_f_P - K*x1malpha_f_P                                           + ...
				          Fdt.Neumann*ExcitationSignal(SignalType,ExciteParam,dt*((i-1)+1-alpha_f))       + ...
				          Fdt.K*ExcitationSignal(SignalType,ExciteParam,dt*((i-1)+1-alpha_f))             + ...
				          Fdt.C*ExcitationSignalDerivative(SignalType,ExciteParam,dt*((i-1)+1-alpha_f),1) + ...
				          Fdt.M*ExcitationSignalDerivative(SignalType,ExciteParam,dt*((i-1)+1-alpha_f),2));
	end
          
    % Acceleration
    xdd1 = (xdd1malpha_m - alpha_m*xdd0)/(1 - alpha_m);
    
    % Velocity
    xd1 = xd0 + dt*((1 - gamma)*xdd0 + gamma*xdd1);
    
    % Displacement
    x1 = x0 + dt*xd0 + dt^2*((0.5 - beta)*xdd0 + beta*xdd1);
    
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
    % Displacement
    x0 = x1;
    % Velocity
    xd0 = xd1;
    % Acceleration
    xdd0 = xdd1;
  
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/NumTimeStep,h,['Explicit generalized alpha method: '...
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
    T = tvec;
else
    U = uvec;
    T = tvec;
end

if nUAll > 0
    save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
              filesep 'UAll.mat'],'-v7.3')
end

end