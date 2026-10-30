function results = FrequencyDomainAnalysis(model,settings)
%FREQUENCYDOMAINANALYSIS Steady-state response to a harmonic excitation.
%
%   results = FrequencyDomainAnalysis(model, settings) solves
%   (K - omega^2*M)*u = F + F_dt with omega = 2*pi*settings.fHarmonic.
%   Damping and non-homogeneous Dirichlet conditions are not taken into
%   account.
%
%   Fields of results:
%     Kdyn       dynamic stiffness matrix K - omega^2*M (all DOFs)
%     K, M, F, F_dt
%     u, ux, uy  displacement amplitudes
%     ActiveDoF  DOFs without homogeneous Dirichlet conditions
%     NOI        nodes of interest marked in the input file
%     t_mat      CPU time of the assembly [s]
%     t_cpu      CPU time of the solution [s]

num_DoF = model.num_DoF;

PrintStep(settings.verbose,'Assembly of the system matrices')
[K,M,~,~,t_mat] = AssembleSystem(model,settings,1);

PrintStep(settings.verbose,'Implement boundary conditions')
[BC_DOFs,NOI] = DirichletBCs_homogeneous(model.boundary_conditions);

PrintStep(settings.verbose,'Apply loads')
[F,F_dt] = LoadVectors(model);

PrintStep(settings.verbose,'Frequency domain analysis - Solving the linear system of equations')
ActiveDoF = setdiff((1:num_DoF)',unique(BC_DOFs));

% Dynamic stiffness matrix for the harmonic excitation frequency
Kdyn = K - (2*pi*settings.fHarmonic)^2*M;
Ka = Kdyn(ActiveDoF,ActiveDoF);
Fa = F(ActiveDoF) + F_dt(ActiveDoF);

u = zeros(num_DoF,1);
tic;
try
    u(ActiveDoF) = Ka\Fa;
catch
    % Fallback: sparse LU with row and column permutation vectors,
    % Ka(r,c) = L*U
    [L,U,r,c] = lu(Ka,'vector');
    uA = zeros(numel(ActiveDoF),1);
    uA(c) = U\(L\Fa(r));
    u(ActiveDoF) = uA;
end
t_cpu = toc;

results.Kdyn      = Kdyn;
results.K         = K;
results.M         = M;
results.F         = F;
results.F_dt      = F_dt;
results.u         = u;
results.ux        = u(1:2:end);
results.uy        = u(2:2:end);
results.ActiveDoF = ActiveDoF;
results.NOI       = NOI;
results.t_mat     = t_mat;
results.t_cpu     = t_cpu;

end
