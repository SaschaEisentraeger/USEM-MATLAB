function results = ModalAnalysis(model,settings)
%MODALANALYSIS Eigenfrequencies and mode shapes.
%
%   results = ModalAnalysis(model, settings) solves K*phi = omega^2*M*phi
%   for the settings.numEigen smallest eigenvalues (eigs with the options
%   settings.EVP_Solver; eig if more than half of the eigenvalues are
%   requested). With Rayleigh damping, the quadratic eigenvalue problem
%   (lambda^2*M + lambda*R + K)*phi = 0 is solved in state-space form.
%   Only homogeneous Dirichlet conditions are taken into account.
%
%   The element matrices are unsymmetric, so eigenvalues and mode shapes
%   can be complex; EF holds the moduli.
%
%   Fields of results:
%     EF         eigenfrequencies [Hz], sorted in ascending order
%     U          mode shapes [num_DoF x numEigen] (fixed DOFs are zero)
%     numEigen   number of computed eigenpairs
%     K, M, R    system matrices (all DOFs)
%     ActiveDoF  DOFs without homogeneous Dirichlet conditions
%     NOI        nodes of interest marked in the input file
%     t_mat      CPU time of the assembly [s]
%     t_cpu      CPU time of the eigenvalue solver [s]

num_DoF = model.num_DoF;

PrintStep(settings.verbose,'Assembly of the system matrices')
[K,M,R,~,t_mat] = AssembleSystem(model,settings,1);

PrintStep(settings.verbose,'Implement boundary conditions')
[BC_DOFs,NOI] = DirichletBCs_homogeneous(model.boundary_conditions);
ActiveDoF = setdiff((1:num_DoF)',unique(BC_DOFs));
NoDoFActive = numel(ActiveDoF);

PrintStep(settings.verbose,'Modal analysis - Solving the eigenvalue problem')

% Options of eigs; the default size of the Krylov subspace depends on the
% requested number of eigenvalues
numEigen = settings.numEigen;
EVP_Solver = settings.EVP_Solver;
if isempty(EVP_Solver.p)
    EVP_Solver.p = max(2*numEigen,20);
end

% The Krylov subspace must be larger than numEigen + 1
if NoDoFActive <= numEigen + 1
    numEigen = NoDoFActive - 2;
end
if numEigen == 0
    numEigen = NoDoFActive;
end
% If more than half of the spectrum is requested, compute all eigenvalues
% with eig (dense)
if 2*numEigen > NoDoFActive
    EVP_Solver.p = NoDoFActive;
    numEigen = NoDoFActive;
    flag_eig = true;
else
    flag_eig = false;
end

Ka = K(ActiveDoF,ActiveDoF);
Ma = M(ActiveDoF,ActiveDoF);

tic;
if nnz(R) == 0 || isempty(R)
    % Undamped problem
    type = EVP_Solver.type;
    opts = rmfield(EVP_Solver,'type');
    if opts.v0 == true
        opts.v0 = randn(NoDoFActive,1);
    else
        opts = rmfield(opts,'v0');
    end
    if flag_eig
        [MS,EV] = eig(full(Ka),full(Ma));
    else
        [MS,EV] = eigs(Ka,Ma,numEigen,type,opts);
    end
    EF = sqrt(abs(diag(EV))/(4*pi*pi));
else
    % Damped problem in state-space form with the eigenvector
    % [lambda*phi; phi]: keep the displacement part and convert |lambda|
    % (rad/s) to Hz
    Z = sparse(NoDoFActive,NoDoFActive);
    A = [-speye(NoDoFActive), Z; Z, Ka];
    B = [Z, speye(NoDoFActive); Ma, R(ActiveDoF,ActiveDoF)];
    [MS,EV] = eigs(A,-B,numEigen,'smallestabs');
    EF = abs(diag(EV))/(2*pi);
    MS = MS(NoDoFActive+1:end,:);
end
[EF,IdxSort] = sort(EF);
MS = MS(:,IdxSort);
t_cpu = toc;

% Mode shapes including the fixed degrees of freedom
U = zeros(num_DoF,numEigen);
U(ActiveDoF,:) = MS;

results.EF        = EF;
results.U         = U;
results.numEigen  = numEigen;
results.K         = K;
results.M         = M;
results.R         = R;
results.ActiveDoF = ActiveDoF;
results.NOI       = NOI;
results.t_mat     = t_mat;
results.t_cpu     = t_cpu;

end
