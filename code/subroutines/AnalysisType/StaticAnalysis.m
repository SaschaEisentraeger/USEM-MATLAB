function results = StaticAnalysis(model,settings)
%STATICANALYSIS Linear static analysis K*u = F.
%
%   results = StaticAnalysis(model, settings) assembles the stiffness
%   matrix and the loads of the model (ReadModel plus the material data of
%   RunSimulation), applies the Dirichlet conditions and solves for the
%   displacements.
%
%   Homogeneous conditions remove their DOFs from the system; prescribed
%   (non-homogeneous) values are applied by elimination in
%   DirichletBCs_nonhomogeneous, which keeps the size of K.
%
%   Fields of results:
%     K, F       stiffness matrix and load vector after the Dirichlet
%                conditions
%     F_dt       time-dependent part of the loads (not used here)
%     u, ux, uy  displacements (all DOFs, x and y components per node)
%     ActiveDoF  DOFs without homogeneous Dirichlet conditions
%     BC_DOFs    DOFs with homogeneous Dirichlet conditions
%     NOI        nodes of interest marked in the input file
%     t_mat      CPU time of the assembly [s]
%     t_cpu      CPU time of the solution [s]

num_DoF = model.num_DoF;

PrintStep(settings.verbose,'Assembly of the system matrices')
[K,~,~,~,t_mat] = AssembleSystem(model,settings,0);

PrintStep(settings.verbose,'Implement homogeneous Dirichlet boundary conditions')
[BC_DOFs,NOI] = DirichletBCs_homogeneous(model.boundary_conditions);

PrintStep(settings.verbose,'Apply loads')
[F,F_dt] = LoadVectors(model);

PrintStep(settings.verbose,'Implement non-homogeneous Dirichlet boundary conditions')
[K,F] = DirichletBCs_nonhomogeneous(K,F,num_DoF,...
        model.boundary_conditions,model.xCoor,model.yCoor);

% Active degrees of freedom
ActiveDoF = setdiff((1:num_DoF)',unique(BC_DOFs));

PrintStep(settings.verbose,'Static analysis - Solving the linear system of equations')
u = zeros(num_DoF,1);
tic;
u(ActiveDoF) = K(ActiveDoF,ActiveDoF)\F(ActiveDoF);
t_cpu = toc;

results.K         = K;
results.F         = F;
results.F_dt      = F_dt;
results.u         = u;
results.ux        = u(1:2:end);
results.uy        = u(2:2:end);
results.ActiveDoF = ActiveDoF;
results.BC_DOFs   = BC_DOFs;
results.NOI       = NOI;
results.t_mat     = t_mat;
results.t_cpu     = t_cpu;

end
