function [model,results,settings] = RunSimulation(settings)
%RUNSIMULATION Run a complete simulation: input, analysis, output.
%
%   [model, results] = RunSimulation(settings) runs the analysis described
%   by the settings (see DefaultSettings for all parameters):
%
%     1. ReadModel        ABAQUS input file -> mesh, materials, loads and
%                         boundary conditions; lifting of the mesh to
%                         spectral elements of degree p
%     2. MaterialMatrices plane elasticity matrices; RayleighDamping
%     3. analysis         StaticAnalysis, DynamicAnalysis, ModalAnalysis
%                         or FrequencyDomainAnalysis
%     4. output           SaveResults (settings.saveWS), PlotResults
%
%   model    struct of ReadModel plus the material data C, rho, zeta and
%            the Rayleigh parameters alpha_Rayleigh, beta_Rayleigh
%   results  struct of the analysis function (see its help)
%
%   [model, results, settings] = RunSimulation(settings) also returns the
%   settings with the fields added during the run (e.g. ResultsFile).
%
%   The code folders must be on the MATLAB path (SetupPaths).

line = repmat('-',1,75);
if settings.verbose
    disp([repmat('-',1,28) 'Start of Simulation!' repmat('-',1,27)])
    disp('Spectral elements')
    disp(line)
end

if settings.parallel
    StartParallelPool();
end

%% Input file and spectral element mesh
PrintStep(settings.verbose,['Read ABAQUS input file ' settings.filename])
t0 = tic;
model = ReadModel(settings.filename,settings.p,settings.verbose);
if settings.verbose
    fprintf('Loading of %s is finished (%.2f s)\n%s\n',settings.filename,...
            toc(t0),line);
end

%% Material data and damping
[model.C,model.rho,model.zeta] = MaterialMatrices(model.material_properties,...
                                 settings.antype);
[model.alpha_Rayleigh,model.beta_Rayleigh] = RayleighDamping(...
    settings.omega_vec,settings.zeta_vec,settings.plot_damping);

if settings.plot_mesh
    % Mesh with loads and Dirichlet conditions (last argument: plot them)
    FE_Mesh(model.p,model.xCoor,model.yCoor,model.Topology_Node,...
            model.cload,model.dsload,model.body_force,model.gravity,...
            model.boundary_conditions,1)
end

%% Analysis
switch settings.analysis_type
    case 'Static_analysis'
        results = StaticAnalysis(model,settings);
        uPlot = results.u;
    case 'Dynamic_analysis'
        results = DynamicAnalysis(model,settings);
        uPlot = [];
    case 'Modal_analysis'
        results = ModalAnalysis(model,settings);
        % First mode (eigenvectors of the unsymmetric or damped problem
        % can be complex)
        uPlot = real(results.U(:,1));
    case 'FrequencyDomain_analysis'
        results = FrequencyDomainAnalysis(model,settings);
        uPlot = results.u;
    otherwise
        error('RunSimulation:analysis_type', ['Unknown analysis_type ' ...
              '"%s" (Static_analysis, Dynamic_analysis, Modal_analysis, ' ...
              'FrequencyDomain_analysis).'],settings.analysis_type);
end

%% Output
if settings.saveWS
    settings.ResultsFile = SaveResults(model,settings,results);
    PrintStep(settings.verbose,['Results saved to ' settings.ResultsFile])
end
if ~isempty(uPlot)
    PlotResults(model,settings,uPlot);
end

if settings.verbose
    disp([repmat('-',1,28) 'End of Simulation!' repmat('-',1,29)])
    disp(line)
end

end
