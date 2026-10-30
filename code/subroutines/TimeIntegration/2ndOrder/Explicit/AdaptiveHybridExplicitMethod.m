%% Soares: Adaptive Hybrid Explicit Method
% 25/08/2020
% UNSW Sydney
% Sascha Eisentraeger
%
% Adapted to the unsymmetric spectral element code (vectorized element
% operations, full unsymmetric matrices, element data from AssemblyFEM).
%-------------------------------------------------------------------------%
% Input
% K:            Stiffness matrix (active degrees of freedom)
% M:            Mass matrix (active degrees of freedom); lumped mass
%               matrices make the effective matrix E diagonal
% C:            Damping matrix (active degrees of freedom)
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
% ElemData:     Element data (see AssemblyFEM):
%               ElemData.Ke, ElemData.Me: element matrices (cell arrays)
%               ElemData.Topo: DOF connectivity [nElem x nDoF per element]
%               ElemData.damped: elements with Rayleigh damping
%               ElemData.alpha, ElemData.beta: Rayleigh parameters
%-------------------------------------------------------------------------%
% Output
% T:            Time vector
% U:            Displacements
% V:            Velocities
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
% Implementation of Table 1 of the paper:
%   regular elements  (Omega_e <= 2): Eqs. (5)-(8), parameter alpha_e
%   enhanced elements (Omega_e >  2): Eqs. (9)-(12), parameter beta_e
% with Omega_e = omega_max,e*dt (highest eigenfrequency of the element).
% The critical sampling frequency is 4, i.e. dt can be twice the critical
% time step of the central difference method.
%
% Element-wise operations sum_e w_e*K_e*x_e are evaluated with the block-
% diagonal matrix of all element matrices Kb and the gather operator G
% (x_e = G*x, assembly = G'): sum_e w_e*K_e*x_e = G'*(w.*(Kb*(G*x))).
% Because alpha_e*K_e = K_e + (alpha_e-1)*K_e, only the deviation from the
% non-dissipative scheme (alpha_e = 1, beta_e = 0) needs element weights.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U,V] = AdaptiveHybridExplicitMethod(K,M,C,F,F_dt,U0,V0,...
                   NumDoF,DOFsave,ActiveDoFs,ExciteParam,SignalType,...
                   dtmax,tsim,ElemData)

% Maximum value of the integration parameter alpha (Section 2: alpha -> inf
% for Omega -> 0; alpha_max ~ 3 avoids excessive algorithmic damping)
alpha_max = 3;

% Number of active degrees of freedom
NumDoFActive = numel(ActiveDoFs);

[K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,NumDoFActive);

% Initial conditions
if numel(U0) > NumDoFActive
    x0  = U0(ActiveDoFs);
    xd0 = V0(ActiveDoFs);
else
    x0  = U0;
    xd0 = V0;
end

%% Time increment
% The stability limit of the adaptive hybrid explicit method is twice that
% of the central difference method (CDM) -->  dt = 2.0*dt_CDM
[dt,NumTimeStep] = TimeStepSize(dtmax,tsim,2.0);

%% Element data
NumElem = numel(ElemData.Ke);
nde = size(ElemData.Topo,2);
if isempty(ElemData.Me)
    error('AdaptiveHybridExplicitMethod:mass', ...
          'The element mass matrices are required (ElemData.Me).');
end

% Gather operator G: element DOFs <-- active global DOFs (constrained DOFs
% are not active; their rows of G are zero)
vecID = zeros(NumDoF,1);
vecID(ActiveDoFs) = 1:NumDoFActive;
TopoActive = vecID(ElemData.Topo);
TopoActive = reshape(TopoActive,NumElem,nde);
colG = reshape(TopoActive',[],1);
rowG = (1:NumElem*nde)';
G = sparse(rowG(colG > 0),colG(colG > 0),1,NumElem*nde,NumDoFActive);

% Block-diagonal matrices of the element stiffness and mass matrices
[ii,jj] = ndgrid(1:nde,1:nde);
offset = (0:NumElem-1)*nde;
I = ii(:) + offset;
J = jj(:) + offset;
Kvals = cell2mat(cellfun(@(A) A(:),ElemData.Ke(:),'UniformOutput',false));
Mvals = cell2mat(cellfun(@(A) A(:),ElemData.Me(:),'UniformOutput',false));
Kb = sparse(I(:),J(:),Kvals,NumElem*nde,NumElem*nde);
Mb = sparse(I(:),J(:),Mvals,NumElem*nde,NumElem*nde);

% Element damping matrices (Rayleigh hypothesis; only damped elements)
damped = logical(ElemData.damped(:));
Dd = spdiags(repelem(double(damped),nde),0,NumElem*nde,NumElem*nde);
Cb = Dd*(ElemData.alpha*Mb + ElemData.beta*Kb);

% Maximum sampling frequency of each element: Omega_e = omega_max,e*dt
% (generalized eigenvalue problem of the element; the eigenvalues of the
% unsymmetric element matrices can be complex)
omega_max = zeros(NumElem,1);
for e = 1:NumElem
    lambda = eig(full(ElemData.Ke{e}),full(ElemData.Me{e}));
    lambda = lambda(isfinite(lambda));
    omega_max(e) = sqrt(max(abs(lambda)));
end
Omega = omega_max*dt;
if any(Omega > 4)
    warning('AdaptiveHybridExplicitMethod:stability', ...
           ['The time step exceeds the stability limit Omega <= 4 ' ...
            'for %d element(s) (max. Omega = %.3f).'], ...
            nnz(Omega > 4), max(Omega));
end

% Damping ratio zeta_e = varsigma_e*dt/(2*rho_e) (Eq. (8b)); for Rayleigh
% damping the mass-proportional part corresponds to varsigma/rho = alpha
zeta = ElemData.alpha*dt/2*double(damped);

% "Regular" (Omega_e <= 2) and "enhanced" (Omega_e > 2) elements
ElemRegular = Omega <= 2;
ElemEnhanced = ~ElemRegular;

% Integration parameters for oscillating elements
% Regular elements, Eq. (8b) (alpha >= 1, alpha <= alpha_max)
alphaElem = 2*sqrt((1+zeta).*(2./Omega).^2) - 1 - zeta.*(2./Omega).^2;
alphaElem = min(max(alphaElem,1),alpha_max);
alphaElem(~ElemRegular) = 1;
% Enhanced elements, Eq. (12b) (beta >= 0)
betaElem = (16*(1 - zeta.*(zeta+4)) - (Omega.^2-8).^2 + ...
            4*(zeta+1).*sqrt(16 - 48*zeta.*(zeta-2) + ...
            2*(Omega.^2-8).^2)).*Omega.^-4;
betaElem = max(betaElem,0);
betaElem(ElemRegular) = 0;

% Global matrices of the enhanced elements (step 5)
De = spdiags(repelem(double(ElemEnhanced),nde),0,NumElem*nde,NumElem*nde);
Kenh = G'*(De*Kb)*G;
Cenh = G'*(De*Cb)*G;
flag_enhanced = any(ElemEnhanced);

% Effective matrix E = M + dt/2*C (diagonal for lumped matrices)
E = M + dt/2*C;
Einv = [];
dE = [];
if isdiag(E)
    flag_diag = true;
    Einv = 1./full(diag(E));
else
    flag_diag = false;
	dE = decomposition(E);
end
solveE = @(b) AhemSolveEffective(b,flag_diag,Einv,dE);

%% Time integration
% Result variables
uvec = zeros(NumTimeStep+1,numel(DOFsave));
vvec = zeros(NumTimeStep+1,numel(DOFsave));
uvec(1,:) = x0(DOFsave);
vvec(1,:) = xd0(DOFsave);

pb = ProgressBar('Adaptive hybrid explicit method - Soares',NumTimeStep);

disp('Solve the equations of motion using the adaptive hybrid explicit method');
disp('');
disp('---------------------------------------------------------------------');
disp('');

% Time vector and excitation signal
T = (0:NumTimeStep)'*dt;
[Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,T);

% Oscillatory criterion of the degrees of freedom, Eq. (7): phi_i is not
% zero if and only if the last two increments of U_i have opposite signs
% (||a+b| - |a| - |b|| = 0 <=> a*b >= 0); this exact test avoids round-off
% in the evaluation of Eq. (7)
xm1 = x0;
oscDoF = false(NumDoFActive,1);

for i = 1:NumTimeStep
    %% Step 1: time integral of the force vector (trapezoidal rule)
    F0 = LoadAtTime(F,Fdt,Amp(i),dAmp(i),ddAmp(i));
    F1 = LoadAtTime(F,Fdt,Amp(i+1),dAmp(i+1),ddAmp(i+1));
    Fbar = 0.5*dt*(F0 + F1);

    % Oscillating elements: at least one of their DOFs oscillates
    oscElem = any(reshape(G*double(oscDoF),nde,NumElem) > 0,1)';

    %% Steps 2-3: regular formulation, vector A
    % alpha_e = 1 unless the element oscillates and is regular
    w = (alphaElem - 1).*oscElem;
    A = Fbar - dt*(C*xd0) - K*(dt*x0 + 0.5*dt^2*xd0);
    if any(w)
        A = A - 0.5*dt^2*(G'*(repelem(w,nde).*(Kb*(G*xd0))));
    end

    %% Step 4: velocity vector
    xd1 = xd0 + solveE(A);

    %% Steps 5-6: enhanced formulation (only enhanced elements), vector B
    if flag_enhanced
        % beta_e = 0 unless the element oscillates (and is enhanced)
        w = betaElem.*oscElem;
        B = -Cenh*xd1 - dt/8*(Kenh*xd1);
        if any(w)
            B = B - dt/8*(G'*(repelem(w,nde).*(Kb*(G*(0.5*xd0 + xd1)))));
        end
        x1 = x0 + 0.5*dt*(xd0 + xd1) + 0.5*dt^2*solveE(B);
    else
        x1 = x0 + 0.5*dt*(xd0 + xd1);
    end

    %% Step 7: oscillatory criterion for the next time step
    oscDoF = (x1 - x0).*(x0 - xm1) < 0;

    %% Step 8: update variables
    xm1 = x0;
    x0  = x1;
	xd0 = xd1;

    % Save variables
    uvec(i+1,:) = x1(DOFsave);
    vvec(i+1,:) = xd1(DOFsave);

    pb.update(i);
end
pb.close();

% Results (all DOFs if every active DOF was saved)
[U,V] = ExpandToAllDoFs(NumDoF,ActiveDoFs,uvec,vvec);

end
