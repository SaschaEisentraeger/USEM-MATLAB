%% Highly precise time integration method
% 12.12.2018
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
%               [fex, nex, fch1, fch2, tch, a, b]
% SignalType:   Type of the excitation signal
% dtmax:        Vector of different allowed maximum time increments
%               [dt_fex,dt_CDM,dt_imp]; the critical time increment of the
% 				central difference method (CDM) is estimated by 2/omega_max
% tsim:         Simulation time
% nSub:			Number of sub-steps
%-------------------------------------------------------------------------%
% Output
% T:            Time vector
% U:            Displacements
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Xing, Y.; Zhang, H. & Wang, Z.: Highly Precise Time Integration Method
% for Linear Structural Dynamic Analysis, International Journal for
% Numerical Methods in Engineering, Wiley, 2018, 116, 505-529
%
% To achieve high accuracy, conventional time integration methods require
% small time steps, especially when applied to large-scale finite element
% systems. The use of small time steps calls for a large number of
% computations. In this paper, a strategy is proposed to construct an
% efficient and highly precise time integration method for linear time-
% invariant systems. Adopting the multi-substep notion, the highly precise
% time integration method creates an integrated amplification matrix prior
% to recursive calculations by multiplying the amplification matrices of N
% equal substeps, where a 2m algorithm and the method of storing
% incremental matrix are respectively employed to reduce computational cost
% and rounding error. On the basis of the Newmark method, two highly
% precise symplectic schemes and one highly precise dissipative scheme are
% developed. In addition, the consistent stability of the used Newmark
% schemes is presented. The free vibration of a three-dimensional truss,
% the impact of a particle against a rod, and the forced vibrations of a
% rectangular thin plate and a stiff system are investigated to validate
% the performance of the proposed method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = HPTIM_Newmark(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim,nSub)

% (De-)activate the waitbar
flag_wb = true;

% Newmark family
% 1: Trapezoidal rule
% 2: Central difference scheme
% 3: Implicit dissipative method
flag_method = 2;

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

%% Highly precise time integration algorithm based on the Newmark family
% Suitable only for the analysis of a second-order homogeneous linear time-
% invariant dynamic systems --> Transform an inhomogeneous system of
% equations

%% Initial operations

%-------------------------------------------------------------------------%
% 1) Adjusted system matrices
%-------------------------------------------------------------------------%
% Harmonic excitation signal
if SignalType == 4
    % Number of basis functions
    q = 1;
    % Amplitude vector/matrix
    T = Fdt.Neumann;
    omega = 2*pi*ExciteParam(1);
    B = 0;
    if ExciteParam(6) ~= 0  % Amplitude of the cosine function
        D  = omega^2;
        f  = @(time) ExciteParam(6)*cos(omega*time) + F;
        fd = @(time) -omega*ExciteParam(6)*sin(omega*time);
    elseif ExciteParam(7) ~= 0 % Amplitude of the sine function
        D  = omega^2;
        f  = @(time) ExciteParam(7)*sin(omega*time) + F;
        fd = @(time) ExciteParam(7)*omega*cos(omega*time);
    end
    % Zero vector/matrix
    Z = spalloc(NumDoFActive,1,0);
    % Identity matrix
    I = speye(q);
    % Adjusted system matrices
    Meff = [M,  Z; Z', I];
    Ceff = [C,  Z; Z', B];
    Keff = [K, -T; Z', D];
else
    % Note that arbitrary excitation functions are almost impossible to
    % implement in this time integration scheme. Not useful for general
    % applications but the general idea might be extendable?!?
    error(['The excitation force needs to be harmonic. ' ...
           'General excitations are not implemented yet'])
end

%-------------------------------------------------------------------------%
% 2) Newmark parameters
%-------------------------------------------------------------------------%
% Method
if flag_method == 1
    % Trapezoidal rule
    beta = 1.0/4.0;
    gamma = 1.0/2.0;
elseif flag_method == 2
    % Central difference scheme
    beta = 0.0;
    gamma = 1.0/2.0;
elseif flag_method == 3
    % Implicit dissipative method
    beta = 3.0/10.0;
    gamma = 11.0/20.0;
end

% Time step size
if flag_method == 1 || flag_method == 3
	% Determine the time increment used for the time integration
	N_dt = numel(unique(dtmax));
	if N_dt == 3
		% Deliberately described time step
		dtN = dtmax(3);
	else
		dtN = dtmax(1);
	end
elseif flag_method == 2
    % Critical time step for the central difference scheme
    dtN = min(dtmax);
end

% Time step size for the HPTIM
dt = dtN*2^nSub;

fprintf('Time step size used: %6.3e\n\n',dt);

%-------------------------------------------------------------------------%
% 3) Initialization
%-------------------------------------------------------------------------%
% Initial conditions
% Displacement
x0  = U0;
% Velocity
xd0 = V0;

% Basis function for the excitation
f0 = f(0);
fd0 = fd(0);
% State space vector
z0 = [x0; f0; xd0; fd0];

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

% Result variables
tvec = zeros(NumTimeStep+1,1);
uvec = zeros(NumTimeStep+1,NumDoFSave);
% Assign the initial conditions to the results
tvec(1)   = 0;
uvec(1,:) = x0(DOFsave);

%-------------------------------------------------------------------------%
% 4) Auxiliary matrices
%-------------------------------------------------------------------------%
% The auxiliary matrices s exhibit a large fill-in and therefore, this
% scheme is not suitable for time-integration of systems with more than a
% couple hundred degrees of freedom.
s1 = Meff\Keff;
s2 = Meff\Ceff;
s3 = (Meff + dtN*gamma*Ceff + dtN^2*beta*Keff)\Meff;
% Number of active degrees of freedom
NoDOF = numel(ActiveDoFs);
% Number of non-zero components in the auxiliary matrix
NNZs1 = nnz(s1);
% If the ratio between the number of non-zero components of s1 and the
% total number of matrix components is higher than 20% it does not make
% sense to use the sparse format as it slows down the computation
% significantly. High demands on the RAM!
if NNZs1/NoDOF > 0.2
    s1 = full(s1);
    s2 = full(s2);
    s3 = full(s3);
end

%-------------------------------------------------------------------------%
% 5) Incremental matrix
%-------------------------------------------------------------------------%
% Identity matrix
I = speye(size(s1));
% Incremental matrix
s11 = s1*s1;
s12 = s1*s2;
s21 = s2*s1;
s22 = s2*s2;
SN11 = -dtN/2*s1 + dtN^2*beta*s3*(s21+dtN/2*s11);
SN12 = I - dtN/2*s2 - dtN^2*beta*s3*(s1-s22-dtN/2*s12);
SN21 = -s1 + dtN*gamma*s3*(s21+dtN/2*s11);
SN22 = -s2 - dtN*gamma*s3*(s1-s22-dtN/2*s12);
SN = dtN*[SN11, SN12; SN21, SN22];

%-------------------------------------------------------------------------%
% 6) Update of the incremental matrix
%-------------------------------------------------------------------------%
for i = 1:nSub
    SN = 2*SN + SN*SN;
end

%-------------------------------------------------------------------------%
% 7) Integrated amplification matrix
%-------------------------------------------------------------------------%
% Identity matrix
I = speye(size(SN));
% Amplification matrix
AN = I + SN;

%% Loop over the time steps
if flag_wb
	% Waitbar
	h = waitbar(0,'HPTIM: ');
end

disp('Solve the equations of motion using HPTIM');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time integration scheme
for i = 1:NumTimeStep
    % Simulation time
    tE = i*dt;
    
    % Calculate the displacement and velocity at step i+1
    z1 = AN*z0;
    
%     % Calculate the acceleration at step i+1
%     xdd1 = -(s1*z1(1:dimSoE/2) + s2*z1(dimSoE/2+1:dimSoE));
    
	if nUAll > 0
        if mod(i,nSave) == 0
            c = i/nSave;
            UAll{c} = z1;
		elseif i == NumTimeStep
            UAll{nUAll} = z1;
        end
	end
	
    % Reset variables for the next time step
    z0 = z1;
    
    % Save variables
    tvec(i+1)   = tE;
    uvec(i+1,:) = z0(DOFsave);
    
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/NumTimeStep,h,['HPTIM: ' ...
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