%% Higher-order implicit time integration (Lagrange poly.): Kim and Reddy
% 30.08.2018
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
% pInt:			Order of the time integration scheme
%-------------------------------------------------------------------------%
% Output
% T:            Time vector
% U:            Displacements
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kim, W. & Reddy, J. N.: A New Family of Higher-Order Time Integration
% Algorithms for the Analysis of Structural Dynamics, Journal of Applied
% Mechanics, ASME International, 2017, 84, 071008
%
% For the development of a new family of implicit higher-order time 
% integration algorithms, mixed formulations that include three time-
% dependent variables (i.e., the displacement, velocity, and acceleration
% vectors) are developed. Equal degree Lagrange type interpolation
% functions in time are used to approximate the dependent variables in the
% mixed formulations, and the time finite element method and the modified
% weighted-residual method are applied to the velocity–displacement and
% velocity–acceleration relations of the mixed formulations. Weight
% parameters are introduced and optimized to achieve preferable attributes
% of the time integration algorithms. Specific problems of structural
% dynamics are used in the numerical examples to discuss some fundamental
% limitations of the well-known second-order accurate algorithms as well as
% to demonstrate advantages of using the developed higher-order algorithms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = KimLagrange_ho(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
				 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim,pInt)

% (De-)activate the waitbar
flag_wb = true;

% Decomposition method: (1) LU decomposition; (2) Matrix
% decomposition (Matlab decides on an appropriate method)
flag_decomp = 2;

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

%% Higher-order time integration algorithm
% Set-up of the parameter of the time integrator
% In the derivation of the method mu is the "ultimate" spectral radius
% (reached for dt/T > 1) which is 0 <= mu <= 1. The choice mu = 1 provides
% an extra order of accuracy while no algorithmic damping is present (2n).
% The maximum amount of numerical damping is achieved for mu = 0 and the
% order of accuracy is in this case 2n-1.
mu = 1.0;

% Order of the interpolation polynomial
IntOrder = pInt;
if IntOrder > 5
	IntOrder = 5;
	disp('The order of the time integrator was set to the maximum value pInt = 5!')
end

% Determine the time increment used for the time integration (the factor
% IntOrder^2 stems from the fact that the system of equations is increased
% compared to the low order case. Due to the Kronecker tensor product the
% system size is increased from [NumDoFActive x NumDoFActive] to
% [NumDoFActive*IntOrder x NumDoFActive*IntOrder]. Therefore, the time step
% size needs to increased as well compared to other implicit algorithms.
% Actually, the proposed increase does not make up for the disadvantage in
% system size.
N_dt = numel(unique(dtmax));
if N_dt == 3
	% Deliberately described time step
	dt = dtmax(3);
else
	dt = dtmax(1);
end
dt = dt*IntOrder^2;

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

% First order accurate algorithm (if mu = 1 the present algorithm is
% identical to the trapezoidal rule --> 'standard' Newmark method)
if IntOrder == 1
    % Parameters for the time integration
    alpha = 1/dt*(1+mu);
    beta  = -1/dt*(1+mu);
    gamma = -mu;
    
% Third order accurate algorithm
elseif IntOrder == 2
    % Parameters for the time integration
    alpha = 1/dt*(mu*[1, -1/4; -4, 1]+[1, 3/4; -4, 3]);
    beta  = 1/dt*(mu*[-3/4; 3]+[-7/4; 1]);
    gamma = mu*[-1/4; 1]+[-1/4; 0];
    
% Fifth order accurate algorithm
elseif IntOrder == 3
    % Parameters for the time integration
    alpha = 1/dt*(mu*[1, (sqrt(5)-3)/2, (sqrt(5)-1)/10;             ...
                      (-sqrt(5)-3)/2, 1, (-sqrt(5)-1)/10;           ...
                      (5*sqrt(5)+5)/2, (-5*sqrt(5)+5)/2, 1] +       ...
                      [(sqrt(5)+3)/2, sqrt(5)-1, (-2*sqrt(5)+3)/5;  ...
                       -sqrt(5)-1, (-sqrt(5)+3)/2, (2*sqrt(5)+3)/5; ...
                       (5*sqrt(5)-5)/2, (-5*sqrt(5)-5)/2, 6]);
    beta  = 1/dt*(mu*[(-3*sqrt(5)+3)/5; (3*sqrt(5)+3)/5; -6] +      ...
                     [(-11*sqrt(5)-11)/10; (11*sqrt(5)-11)/10; -1]);
    gamma = mu*[(-sqrt(5)+1)/10; (sqrt(5)+1)/10; -1] +              ...
               [(-sqrt(5)-1)/10; (sqrt(5)-1)/10; 0];
    
% Seventh order accurate algorithm
elseif IntOrder == 4
    % Parameters for the time integration
    alpha = 1/dt*(mu*[1, (8*sqrt(21)-56)/49, (-sqrt(21)+5)/2, (3*sqrt(21)-21)/98;                    ...
                      (-7*sqrt(21)-49)/32, 1, (7*sqrt(21)-49)/32, 3/16;                              ...
                      (sqrt(21)+5)/2, (-8*sqrt(21)-56)/49, 1, (-3*sqrt(21)-21)/98;                   ...
                      (-7*sqrt(21)-49)/6, 16/3, (7*sqrt(21)-49)/6, 1] +                              ...
                      [(sqrt(21)+5)/2, (88*sqrt(21)-168)/147, (-sqrt(21)+3)/3, (-12*sqrt(21)+63)/49; ...
                       (-77*sqrt(21)-147)/96, 1, (77*sqrt(21)-147)/96, -9/16;                        ...
                       (sqrt(21)+3)/3, (-88*sqrt(21)-168)/147, (-sqrt(21)+5)/2, (12*sqrt(21)+63)/49; ...
                       (7*sqrt(21)-49)/6, 16/3, (-7*sqrt(21)-49)/6, 10]);
    beta  = 1/dt*(mu*[(15*sqrt(21)-105)/49; 15/8; (-15*sqrt(21)-105)/49; 10] +                       ...
                     [(-51*sqrt(21)-357)/98; 21/8; (51*sqrt(21)-357)/98; 1]);
    gamma = mu*[(3*sqrt(21)-21)/98; 3/16; (-3*sqrt(21)-21)/98; 1] +                                  ...
               [(-3*sqrt(21)-21)/98; 3/16; (3*sqrt(21)-21)/98; 0];

% Ninth order accurate algorithm
elseif IntOrder == 5
	% Parameters for the time integration
    alpha = 1/dt*(mu*[ 1.0,         -0.3979905070,  0.2213383872, -0.1331089586,  0.04930278012; ...
                      -2.512622744,  1.0,          -0.5561398661,  0.3344525970, -0.1238792867;  ...
                       4.517969126, -1.798108823,   1.0,          -0.6013821654,  0.2227484384;  ...
                      -7.512642352,  2.989960338,  -1.662836143,   1.0,          -0.3703941540;  ...
                       20.28283188, -8.072374540,   4.489369296,  -2.699826628,   1.0] +         ...
					 [ 7.512642350,  2.056893217,  -0.6428201742,  0.3070950148, -0.1051682020;  ...
                      -7.960483029,  1.798108823,   2.505923933,  -0.9713311794,  0.3166527833;  ...
                       4.085336090, -4.505923934,   0.5561398660,  3.108061311,  -0.8460228092;  ...
                      -2.307095015,  2.703646823,  -5.268191942,   0.1331089586,  3.523427116;   ...
                       2.699826629, -4.489369295,   8.072374541,  -20.28283186,   15.0]);
    beta  = 1/dt*(mu*[-0.7395417018; 1.8581893005; -3.341226576; 5.55591231;  -15.0] + ...
	                 [-9.128642207;  4.311128671;  -2.397590522; 1.215104058; -1.0]);
    gamma = mu*[-0.04930278012; 0.1238792867; -0.2227484384; 0.370394154;   -1.0] + ...
	           [-0.3703941540;  0.2227484384; -0.1238792867; 0.04930278012;  0.0];
end

% Time nodes
[tau,~,~] = lglnodes(IntOrder);
% Map tau from the interval [-1,1] to the interval [0,1]
tau = (tau + 1)/2;
% Change the order of the elements in tau
tau = flipud(tau);

% Initial conditions
% Displacement
x0  = U0;
% Velocity
xd0 = V0;

% Result variables
Tvec = zeros(IntOrder*NumTimeStep+1,1);
Umat = zeros(IntOrder*NumTimeStep+1,NumDoFSave);
% Assign the initial conditions to the results
Tvec(1) = 0;
Umat(1,:) = x0(DOFsave);

% To set up the final system new parameters based on alpha, beta, and gamma
% are defined
alpha_bar = alpha*alpha;
beta_bar  = alpha*beta;
gamma_bar = alpha*gamma;

% Effective (dynamic) stiffness matrix)
Keff = kron(alpha_bar,M) + kron(alpha,C) + kron(speye(IntOrder),K);

% LU-decomposition of the effective stiffness matrix
% LTM: Lower triangular matrix 
% UTM: Upper triangular matrix
% A Cholesky decomposition is not possible since the effective stiffness
% matrix is not symmetric due to the parameter matrices and vectors
if flag_decomp == 1
    % LU decomposition of the stiffness matrix
    [LTM,UTM,row_order,column_order] = lu(Keff,'vector');
elseif flag_decomp == 2
    % Matrix decomposition for solving linear systems
    dK_eff = decomposition(Keff);
end

if flag_wb
	% Waitbar
	h = waitbar(0,'Higher-order time FEM (Lagrange, implicit): ');
end

disp('Solve the equations of motion using the higher-order time FEM (Lagrange)');
disp('');
disp('---------------------------------------------------------------------');
disp('');

if NumDoFActive ~= NumDoFSave
    % Save degrees of freedom for the 'sub-steps' (increased system size)
    % The size of the system is increased by a factor of IntOrder and
    % therefore, additional degrees of freedom have to be saved
    addDOF = (0:IntOrder-1)*NumDoFActive;
    % Repeat the degrees of freedom that have to be saved
    DOFint = repmat(DOFsave',IntOrder,1);
    % Repeat the additional degrees of freedom
    addDOFmat = repmat(addDOF,NumDoFSave,1);
    % Save in a vector
    addDOFvec = reshape(addDOFmat,[],1) + DOFint;
    % Save only selected degrees of freedom
    flag_save = 1;
else
    % Save all degrees of freedom
    flag_save = 0;
end

% Time integration scheme
for i = 1:NumTimeStep
    % Simulation time
    t = (i-1)*dt;
    t_vec = t + tau*dt;
    
    % Loading
    amp_t   = ExcitationSignal(SignalType,ExciteParam,t_vec)';
    damp_t  = ExcitationSignalDerivative(SignalType,ExciteParam,t_vec,1)';
    ddamp_t = ExcitationSignalDerivative(SignalType,ExciteParam,t_vec,2)';
    f0 = F + amp_t(1)*Fdt.Neumann + amp_t(1)*Fdt.K + damp_t(1)*Fdt.C + ...
         ddamp_t(1)*Fdt.M ;
    f1 = kron(ones(IntOrder,1),F) + ...
         kron(amp_t(2:IntOrder+1),Fdt.Neumann) + ...
         kron(amp_t(2:IntOrder+1),Fdt.K)       + ...
         kron(damp_t(2:IntOrder+1),Fdt.C)      + ...
         kron(ddamp_t(2:IntOrder+1),Fdt.M)     ;

	% Effective dynamic force vector
	Feff = f1 - ...
          (kron(beta_bar,M) + kron(beta,C) - kron(gamma,K))*x0 - ...
  	       kron(gamma_bar+beta,M)*xd0 - ...
           kron(gamma,speye(NumDoFActive))*f0;
    
    % In each time step a system of equations has to be solved:
    % Keff*x1 = Feff --> x1 = Keff\Feff;
    % Since the inversion of Keff will result in a fully populated matrix a
    % LU-decomposition is performed instead. The solution is computed by
    % means of forward and backward substitution. Solve the system of
    % equations at the time t + tau*dt to compute the displacement
	if flag_decomp == 1
        % LU decomposition
        x1 = zeros(size(Feff));
		x1(column_order) = UTM\(LTM\Feff(row_order));
    elseif flag_decomp == 2
        % Matrix decomposition
        x1 = dK_eff\Feff;
	end
	
    % Compute the velocity components
    xd1 = kron(alpha,speye(NumDoFActive))*x1 + ...
          kron(beta ,speye(NumDoFActive))*x0 + ...
          kron(gamma,speye(NumDoFActive))*xd0;
    
	if flag_save == 1
        % Save the results - selected degrees of freedom
        Tvec((i-1)*IntOrder+1:i*IntOrder+1) = t_vec;
        Umat((i-1)*IntOrder+1:i*IntOrder+1,:) = ...
        [x0(DOFsave)';reshape(x1(addDOFvec),NumDoFSave,IntOrder)'];
    else
        % Save the results - all degrees of freedom
        Tvec((i-1)*IntOrder+1:i*IntOrder+1) = t_vec;
        Umat((i-1)*IntOrder+1:i*IntOrder+1,:) = ...
        [x0';reshape(x1,NumDoFActive,IntOrder)'];
	end
    
	if nUAll > 0
        if mod(i,nSave) == 0
            c = i/nSave;
            UAll{c} = x1;
		elseif i == NumTimeStep
            UAll{nUAll} = x1;
        end
	end
	
    % Update the variables for the next time step
    x0  = x1(1+(IntOrder-1)*NumDoFActive:IntOrder*NumDoFActive);
    xd0 = xd1(1+(IntOrder-1)*NumDoFActive:IntOrder*NumDoFActive);
    
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/NumTimeStep,h,['Higher-order time FEM (Lagrange, '...
			'implicit): ' num2str(i/NumTimeStep*100,'%5.1f') '%'])
		end
	end
end

if flag_wb
	% Delete the waitbar
	close(h)
end

% Results
if NumDoFActive == NumDoFSave
    U = zeros(IntOrder*NumTimeStep+1,NumDoF);
    U(:,ActiveDoFs) = Umat;
    T = Tvec;
else
    U = Umat;
    T = Tvec;
end

if nUAll > 0
    save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
              filesep 'UAll.mat'],'-v7.3')
end

end