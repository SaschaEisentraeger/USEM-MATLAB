%% Soares: Adaptive Hybrid Explicit Method
% 25/08/2020
% UNSW Sydney
% Sascha Eisentraeger
%-------------------------------------------------------------------------%
% Input
% Mg:           Global mass matrix
% Cg:           Global damping matrix
% K:            Cell array containing the element stiffness matrices
% M:            Cell array containing the element mass matrices
% F:            Force vector
% F_dt          Force vector (time-dependent)
% U0:			Initial displacement vector
% V0:			Initial velocity vector
% TopologyN:    Topology/connectivity data
% Rayleigh:     Damping parameter
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
% Soares, D.: A novel time-marching formulation for wave propagation
% analysis: The adaptive hybrid explicit method, Computer Methods in
% Applied Mechanics and Engineering, 2020, 366, 113095
%
% This paper proposes an adaptive hybrid approach to explicitly analyse
% hyperbolic models. The new methodology considers the combination of
% simple and enhanced explicit time-marching procedures, which are both
% based on adaptive time integrators. By applying this combined approach,
% the hybrid formulation allows avoiding solver routines (explicit
% analysis), enables considering larger time-step values (higher stability
% limits are provided by the method) and permits reduced computational
% efforts, per time step, to be obtained. Thus, the novel technique may be
% considered very efficient. In addition, the proposed hybrid approach is
% entirely automated, providing automatic distributions of regular and
% enhanced subdomains along the discretized model, demanding no expertise
% and/or effort from the user. The time integrators of the method are
% adaptively evaluated, taking into account the local properties of the
% discretized domain and the behaviour of the computed responses, enabling
% a very versatile solution methodology. In this context, a link between
% the temporal and the spatial discretizations may be established, as well
% as optimized algorithmic dissipation can be introduced into the analyses,
% providing enhanced accuracy. Numerical results are presented at the end
% of the manuscript, illustrating the excellent performance of the proposed
% adaptive hybrid explicit method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = AdaptiveHybridExplicitMethod(Mg,Cg,K,M,F,Fdt,U0,V0,...
                 TopologyN,Rayleigh,DOFsave,ActiveDoFs,ExciteParam,...
				 SignalType,dtmax,tsim)

% (De-)activate the waitbar
flag_wb = true;

% Damping parameter
RA = Rayleigh.alpha;
RB = Rayleigh.beta;
if RA == 0 && RB == 0
    flag_damping = false;
else
    flag_damping = true;
end

% Maximum value for the time integration parameters alpha/beta (numerical
% damping)
alpha_max = 3;
beta_max = 2;
             
% Connectivity data
Topo = TopologyN.TopologyD;

% Number of elements
NumElem = size(K,1);

% Number of degrees of freedom
NumDoF = max(cellfun(@(x) max(x),Topo));

% Number of active degrees of freedom
NumDoFActive = numel(ActiveDoFs);

% Number of degrees of freedom that are saved
NumDoFSave = numel(DOFsave);

% Initial conditions
if numel(U0) > NumDoFActive
    % Displacement
    x0  = U0(ActiveDoFs);
    % Velocity
    xd0 = V0(ActiveDoFs);
else
    % Displacement
    x0  = U0;
    % Velocity
    xd0 = V0;
end


% Determine the time increment used for the time integration:
% The stability limit of the adaptive hybrid explicit method is twice that
% of the central difference method (CDM) -->  dt = 2.0*dt_CDM
dtmax(2) = 2.0*dtmax(2);
dt = min(dtmax);

% Initialize the oscillatory parameter
Phi = zeros(NumDoF,1);

% Damping related value
zeta = Rayleigh.alpha*dt/2;

% Data concerning the parallel pool
MatlabPool = gcp('nocreate');
if isempty(MatlabPool)
    % Number of cores is a power of 2 in this implementation
    PowerNC = 2;
    % Number of cores
    NCores = 2^PowerNC;
    % Control maximum number of computational threads
    maxNumCompThreads(NCores);
    % Start the parallel pool
    MatlabPool = parpool(NCores);
end

% Number of element per core
NumElemCore = NumElem/MatlabPool.NumWorkers;
% Number of elements per core when a parfor- is used instead of a for-loop
NELEM = 2e4;

% Topology/connectivity data (including only active degrees of freedom)
TopoGlobalActive = cell(NumElem,1);
% Maximum elemental eigenfrequency
Omega_max = zeros(NumElem,1);
% Constrained element matrices
KeC = cell(NumElem,1);
CeC = cell(NumElem,1);
parfor count = 1:NumElem
    topo = Topo{count};
    %% Topology data
    [Lia,Locb] = ismember(ActiveDoFs,topo);
    % Local index
    ID = find(Lia);
    % Ordering
    Locb(Locb == 0) = [];
    % Global degree of freedom numbers (only active DoFs)
    TopoGlobalActive{count} = zeros(size(topo'));
    TopoGlobalActive{count}(Locb) = ID;
    TopoGlobalActive{count}(TopoGlobalActive{count} == 0) = [];
    
    %% Eigenfrequency
    % Determine the active degrees of freedom per element
    ActiveDoFsElemLocal = find(ismember(topo,ActiveDoFs))';
    % Element matrices (incl. Dirichlet BCs)
    Ke = K{count}(ActiveDoFsElemLocal,ActiveDoFsElemLocal);
    Me = M{count}(ActiveDoFsElemLocal,ActiveDoFsElemLocal);
    % Largest eigenvalue
    lambda_max = eigs(Ke,Me,1,'largestabs');
    % Largest angular frequency
    omega_max = sqrt(lambda_max);
    % Maximum sampling frequency
    Omega_max(count) = omega_max*dt;
    
    %% Damping matrix
    % Compute the damping matrices by means of the Rayleigh hypothesis
  	Ce = RA*M{count}(ActiveDoFsElemLocal,ActiveDoFsElemLocal) + ...
         RB*K{count}(ActiveDoFsElemLocal,ActiveDoFsElemLocal);
    
    % Save constrained element matrices
    KeC{count} = Ke;
    CeC{count} = Ce;
end

% Global matrices are typically stored as a sparse triangular matrix
Mg = Mg + Mg';
Cg = Cg + Cg';

% Diagonal effective matrix E (for efficiency reasons Mg and Cg should be
% diagonal/lumped matrices)
E = Mg + dt/2*Cg;
if isdiag(E)
    flag_diag = true;
    Einv = 1./diag(E);
else
    flag_diag = false;
	% Matrix decomposition for solving linear systems
    dE = decomposition(E);
end

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

% Result variables
uvec = zeros(NumTimeStep+1,NumDoFSave);
% Assign the initial conditions to the results
uvec(1,:) = x0(DOFsave);

if flag_wb
	% Waitbar
	h = waitbar(0,'Adaptive hybrid explicit method - Soares: ');
end

disp('Solve the equations of motion using the Soares method');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Save the displacement values for the last three time steps --> to
% calculate the oscillatory criterion
xPHI = zeros(NumDoF,3);

% Indices to update vector A
IDvec = cell2mat(TopoGlobalActive);
% Number of DoF per element
ncols = cellfun('size', TopoGlobalActive,1);
ncols_sum = cumsum(ncols);
IDstart = [1; ncols_sum(1:NumElem-1)+1];
IDend   = ncols_sum;

% "Regular" elements, i.e., elements with a maximum sampling frequency of
% Omega_max < 2
ElemRegular = Omega_max < 2;
% Indices to update vector B
IDvecEnhanced = cell2mat(TopoGlobalActive(~ElemRegular));

% Auxiliary displacement and velocity vectors (save for each element)
x0Elem  = cell(NumElem,1);
xd0Elem = cell(NumElem,1);
xd1Elem = cell(NumElem,1);

% Integration Parameter
% Regular
alphaElem = 2*sqrt((1+zeta)*(2./Omega_max).^2) - 1 - zeta*(2./Omega_max).^2;
alphaElem(alphaElem > alpha_max) = alpha_max;
alphaElem(alphaElem < 1) = 1;
% Enhanced
betaElem = (16*(1-zeta*(zeta+4)) - (Omega_max.^2-8).^2 + ...
            4*(zeta+1)*sqrt(16-48*zeta*(zeta-2)+         ...
            2*(Omega_max.^2-8).^2)).*Omega_max.^-4;
betaElem(betaElem > beta_max) = beta_max;
betaElem(ElemRegular) = 0;

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
    %% Step 1
    % Compute the effective force vector at the time t + 1/2*dt
    F0 = F + FdtN*Amp(i)   + FdtK*Amp(i)   + FdtC*dAmp(i)   + FdtM*ddAmp(i);
    F1 = F + FdtN*Amp(i+1) + FdtK*Amp(i+1) + FdtC*dAmp(i+1) + FdtM*ddAmp(i+1);
    Fbar = 0.5*dt*F0 + 0.5*dt*F1;
    
    %% Step 2
    % Initialize auxiliary vectors
    A = Fbar;
    B = zeros(NumDoFActive,1);
    
    %% Step 3 - Regular formulation
    % Compute the oscillatory criterion for each element
    PHI = zeros(NumElem,1);
    % Update the displacement and velocity vectors for each element
    % (current time step)
    for j = 1:NumElem
        % Topology data (global)
        topoG = IDvec(IDstart(j):IDend(j));
        % Oscillatory criterion
        PHI(j) = sum(Phi(topoG));
        % Displacement and velocity vectors for each element
        x0Elem{j}  = x0(topoG);
        xd0Elem{j} = xd0(topoG);
    end
    
    % Force terms
    Fstar = cell(NumElem,1);
    if NumElemCore > NELEM
        parfor j = 1:NumElem
            % Oscillatory criterion
            PHIval = PHI(j);
            % Compute alpha value for each element
            if PHIval == 0 || Omega_max(j) > 2
                alpha = 1;
            else
                alpha = alphaElem(j);
            end
            
            % Element stiffness matrix (constrained)
            Ke = KeC{j};
            % Element damping matrix (constrained)
            if flag_damping
                Ce = CeC{j};
            else
                Ce = [];
            end
            % Element displacements and velocities
            x0e  = x0Elem{j};
            xd0e = xd0Elem{j};
            % Force term
            if flag_damping
                Fstar{j} = Ce*1/2*dt*xd0e + Ke*(dt*x0e+1/2*alpha*dt^2*xd0e);
            else
                Fstar{j} = Ke*(dt*x0e+1/2*alpha*dt^2*xd0e);
            end
        end
    else
        for j = 1:NumElem
            % Oscillatory criterion
            PHIval = PHI(j);
            % Compute alpha value for each element
            if PHIval == 0 || Omega_max(j) > 2
                alpha = 1;
            else
                alpha = alphaElem(j);
            end

            % Element stiffness matrix (constrained)
            Ke = KeC{j};
            % Element damping matrix (constrained)
            if flag_damping
                Ce = CeC{j};
            else
                Ce = [];
            end
            % Element displacements and velocities
            x0e  = x0Elem{j};
            xd0e = xd0Elem{j};
            % Force term
            if flag_damping
                Fstar{j} = Ce*dt*xd0e + Ke*(dt*x0e+1/2*alpha*dt^2*xd0e);
            else
                Fstar{j} = Ke*(dt*x0e+1/2*alpha*dt^2*xd0e);
            end
        end
    end
    
    % Update vector A --> assemble element contributions
    Valvec = cell2mat(Fstar);
    Fbar = accumarray(IDvec,Valvec,[NumDoFActive 1]);
    A = A - Fbar;
    
    %% Step 4
    % Compute the velocity vector
    if flag_diag
        xd1 = xd0 + Einv.*A;
    else
        xd1 = xd0 + dE\A;
    end
    
    %% Step 5 - Enhanced formulation
    % This part of the algorithm is only executed for "enhanced" elements,
    % i.e., elements with a maximum sampling frequency Omega_max > 2.
    % Update the velocity vectors for each element (next time step)
    for j = 1:NumElem
        % Topology data (global)
        topoG = IDvec(IDstart(j):IDend(j));
        % Velocity vectors for each element
        xd1Elem{j} = xd1(topoG);
    end
    
    % Force terms
    Fstar = cell(NumElem,1);
    if NumElemCore > NELEM
        parfor j = 1:NumElem
            % Execute only for "enhanced" elements
            if ElemRegular(j)
                continue
            end
            % Oscillatory criterion
            PHIval = PHI(j);
            % Compute beta value for each element
            if PHIval == 0
                beta = 0;
            else
                beta = betaElem(j)
            end

            % Element stiffness matrix (constrained)
            Ke = KeC{j};
            % Element damping matrix (constrained)
            if flag_damping
                Ce = CeC{j};
            else
                Ce = [];
            end
            % Element velocities
            xd0e = xd0Elem{j};
            xd1e = xd1Elem{j};
            % Force term
            if flag_damping
                Fstar{j} = Ce*xd1e + Ke*1/8*dt*(1/2*beta*xd0e+(1+beta)*xd1e);
            else
                Fstar{j} = Ke*1/8*dt*(1/2*beta*xd0e+(1+beta)*xd1e);
            end
        end
    else
        for j = 1:NumElem
            % Execute only for "enhanced" elements
            if ElemRegular(j)
                continue
            end
            % Oscillatory criterion
            PHIval = PHI(j);
            % Compute beta value for each element
            if PHIval == 0
                beta = 0;
            else
                beta = betaElem(j);
            end

            % Element stiffness matrix (constrained)
            Ke = KeC{j};
            % Element damping matrix (constrained)
            if flag_damping
                Ce = CeC{j};
            else
                Ce = [];
            end
            % Element velocities
            xd0e = xd0Elem{j};
            xd1e = xd1Elem{j};
            % Force term
            if flag_damping
                Fstar{j} = Ce*dt*xd1e + Ke*1/8*dt^2*(1/2*beta*xd0e+(1+beta)*xd1e);
            else
                Fstar{j} = Ke*1/8*dt^2*(1/2*beta*xd0e+(1+beta)*xd1e);
            end
        end
    end
    
    % Update vector B --> assemble element contributions
    Valvec = cell2mat(Fstar);
    Fbar = accumarray(IDvecEnhanced,Valvec,[NumDoFActive 1]);
    B = B - Fbar;
    
    %% Step 6
    % Compute the displacement vector
    if flag_diag
        x1 = x0 + 1/2*dt*(xd0+xd1) + 1/2*dt*Einv.*B;
    else
        x1 = x0 + 1/2*dt*(xd0+xd1) + 1/2*dt*dE\B;
    end
    
    %% Step 7
    % Displacement field for the last three time steps
    xPHI(:,1)          = xPHI(:,2);
    xPHI(:,2)          = xPHI(:,3);
    xPHI(ActiveDoFs,3) = x1;
    
    % Update the oscillatory criterion
    Phi = abs(abs(xPHI(:,3)-xPHI(:,1)) - ...
              abs(xPHI(:,3)-xPHI(:,2)) - ...
              abs(xPHI(:,2)-xPHI(:,1)));
    
    %% Step 8 - Update variables
    % Reset variables for the next time step
    x0  = x1;
	xd0 = xd1;
	
	if nUAll > 0
        if mod(i,nSave) == 0
            c = i/nSave;
            UAll{c} = x1;
		elseif i == NumTimeStep
            UAll{nUAll} = x1;
        end
	end
    
    %% Step 9
	% Save variables
    uvec(i+1,:) = x1(DOFsave);
    
	if flag_wb
		% Waitbar
		if mod(i,nUpdate) == 0
			waitbar(i/NumTimeStep,h,['Adaptive hybrid explicit method - ' ... 
			'Soares: ' num2str(i/NumTimeStep*100,'%5.1f') '%'])
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