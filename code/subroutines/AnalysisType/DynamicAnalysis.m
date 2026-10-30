function results = DynamicAnalysis(model,settings)
%DYNAMICANALYSIS Transient analysis by direct time integration.
%
%   results = DynamicAnalysis(model, settings) assembles K, M and the
%   Rayleigh damping matrix R, sets up the loads
%       F + F_dt*s(t)   (s: excitation signal, settings.SignalType)
%   and the lifting of prescribed displacements, and integrates
%       M*a + R*v + K*u = f(t)
%   from rest with the method settings.TimeIntegration.Method
%   (Solver_Transient).
%
%   Time step: the solver receives three candidates [dtf, dte, dtImp]:
%     dtf    0.05/f_max, f_max the highest frequency of the signal
%            (about 20 samples per period)
%     dte    stability limit of the central difference method from the
%            largest element eigenvalue (AssemblyFEM)
%     dtImp  settings.TimeIntegration.dtImp (implicit methods)
%   Each method chooses from them (see the method's documentation).
%
%   Fields of results:
%     t_comp     time instants [nt x 1]
%     y_comp     displacements [nt x nSaved]: all DOFs (save_var = 1) or
%                the DOFs of the nodes of interest (save_var = 0)
%     K, M, R    system matrices (all DOFs)
%     F, F_dt    load vectors; F_dt has the fields Neumann, K, C, M
%                (multiplied by s, s, ds/dt and d2s/dt2)
%     ActiveDoF, NOI, DOF_C, DOF_C_Dyn
%     dte, dtf   time step candidates (see above)
%     t_end      simulation time [s]
%     t_mat      CPU time of the assembly [s]
%     t_cpu      CPU time of the time integration [s]

num_DoF = model.num_DoF;
TimeIntegration = settings.TimeIntegration;

% Signal parameters (defaults if none are given)
ExciteParam = settings.ExciteParam;
if isempty(ExciteParam)
    ExciteParam = DefaultExcitationParameters(settings.SignalType,settings.tE);
end
% Maximum time increment of the excitation (at least 20 samples for the
% highest frequency of the signal)
dtf = 0.05/FFT_Signal(ExciteParam,settings.SignalType);

% Simulation time: from the settings or from the step of the input file
if isempty(settings.t_end)
    t_end = model.t_end;
else
    t_end = settings.t_end;
end

PrintStep(settings.verbose,'Assembly of the system matrices')
% The adaptive hybrid explicit method also needs the element matrices
needElemData = TimeIntegration.Method == 109 && ~settings.flag_mRHS;
if needElemData
    [K,M,R,dte,t_mat,ElemData] = AssembleSystem(model,settings,2);
else
    [K,M,R,dte,t_mat] = AssembleSystem(model,settings,2);
end

PrintStep(settings.verbose,'Implement homogeneous Dirichlet boundary conditions')
[BC_DOFs,NOI] = DirichletBCs_homogeneous(model.boundary_conditions);

PrintStep(settings.verbose,'Apply loads')
[F,FNeumann_dt] = LoadVectors(model,settings.flag_mRHS);

PrintStep(settings.verbose,'Implement non-homogeneous Dirichlet boundary conditions')
[DOF_C,DOF_C_Dyn,FDBCdyn,FDBCdyn_dt] = ...
    DirichletBCs_nonhomogeneous_dynamic(K,R,M,model.Topology_Node,...
    model.boundary_conditions,model.xCoor,model.yCoor);

% Loads: constant part and the parts multiplied by s(t) and its
% derivatives (lifting of the time-dependent prescribed displacements)
F = F + FDBCdyn;
F_dt.Neumann = FNeumann_dt;
F_dt.K       = FDBCdyn_dt.K;
F_dt.C       = FDBCdyn_dt.C;
F_dt.M       = FDBCdyn_dt.M;

% Active degrees of freedom
ActiveDoF = setdiff((1:num_DoF)',unique([BC_DOFs;DOF_C.BC;DOF_C_Dyn.BC]));

% Nodes of interest: the first three active nodes if none are marked
if isempty(NOI) && settings.save_var == 0
    NOI = unique(round(ActiveDoF/2));
    NOI = NOI(1:min(3,numel(NOI)));
end

%% Input of the transient solver
S.K = K;
S.M = M;
S.R = R;
S.F = F;
S.F_dt = F_dt;
S.U0 = zeros(num_DoF,1);
S.V0 = zeros(num_DoF,1);
S.ActiveDoF = ActiveDoF;
S.NOI = NOI;
S.num_DoF = num_DoF;
S.save_var = settings.save_var;
S.t_end = t_end;
S.dtf = dtf;
S.dte = dte;
S.SignalType = settings.SignalType;
S.ExciteParam = ExciteParam;
S.TimeIntegration = TimeIntegration;
S.flag_mRHS = settings.flag_mRHS;
if needElemData
    S.ElemData = ElemData;
end

% Keep the input so that Solver_Transient can be rerun on its own:
%   [T,U] = Solver_Transient(fullfile(settings.TempFolder,'Solver_input.mat'))
solverFile = fullfile(settings.TempFolder,'Solver_input.mat');
if TimeIntegration.SolverSave
    PrintStep(settings.verbose,'Save input data for the solver')
    if ~exist(settings.TempFolder,'dir')
        mkdir(settings.TempFolder);
    end
    save(solverFile,'-struct','S','-v7.3')
elseif isfile(solverFile)
    delete(solverFile);
end

%% Time integration
[t_comp,y_comp,t_cpu] = Solver_Transient(S);

% Add the prescribed displacements (only if all DOFs are saved)
if ~isempty(DOF_C_Dyn.BC) && settings.save_var ~= 0
    Upre = DOF_C_Dyn.Amp*ExcitationSignal(settings.SignalType,ExciteParam,t_comp);
    y_comp(:,DOF_C_Dyn.BC) = Upre';
end
if ~isempty(DOF_C.BC) && settings.save_var ~= 0
    Upre = repmat(DOF_C.Amp,1,numel(t_comp));
    y_comp(:,DOF_C.BC) = Upre';
end

results.t_comp      = t_comp;
results.y_comp      = y_comp;
results.K           = K;
results.M           = M;
results.R           = R;
results.F           = F;
results.F_dt        = F_dt;
results.ActiveDoF   = ActiveDoF;
results.NOI         = NOI;
results.DOF_C       = DOF_C;
results.DOF_C_Dyn   = DOF_C_Dyn;
results.dte         = dte;
results.dtf         = dtf;
results.ExciteParam = ExciteParam;
results.t_end       = t_end;
results.t_mat       = t_mat;
results.t_cpu       = t_cpu;

end
