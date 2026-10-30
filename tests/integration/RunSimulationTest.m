classdef RunSimulationTest < matlab.unittest.TestCase
    %RUNSIMULATIONTEST The complete program flow of Start_Sim.m.
    %
    %   RunSimulation reads an example ABAQUS input file and runs one of the
    %   four analyses. The tests check the static result against the stored
    %   reference solution, the modal and frequency-domain results against
    %   direct computations with the returned matrices, the transient
    %   analysis for consistency, and the saved files.

    properties
        rootDir
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            tc.rootDir = root;
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'code'), 'IncludingSubfolders', true));
            old = get(groot, 'DefaultFigureVisible');
            set(groot, 'DefaultFigureVisible', 'off');
            tc.addTeardown(@() set(groot, 'DefaultFigureVisible', old));
        end
    end

    methods (Test)
        function staticAnalysisMatchesReference(tc)
            % Cubic Lee & Bathe patch test on five distorted elements, p = 3
            s = quietSettings('Static_analysis');
            [model, R] = run(s);
            ref = load(fullfile(tc.rootDir, 'tests', 'fixtures', 'reference', ...
                'PatchTestLeeCubic_Rajendran_adj_ndof104.mat'));
            tc.verifyEqual(model.num_DoF, 104);
            tc.verifyEqual(R.u, ref.u, 'AbsTol', 1e-8*max(abs(ref.u)));
            tc.verifyEqual([R.ux R.uy], [R.u(1:2:end) R.u(2:2:end)]);
        end

        function modalAnalysisMatchesDenseEigenvalues(tc)
            s = quietSettings('Modal_analysis');
            s.p.pxi = 2; s.p.peta = 2;
            s.numEigen = 5;
            [~, R] = run(s);
            a = R.ActiveDoF;
            lambda = eig(full(R.K(a,a)), full(R.M(a,a)));
            f = sort(sqrt(abs(lambda))/(2*pi));
            tc.verifySize(R.EF, [5 1]);
            tc.verifyTrue(issorted(R.EF));
            tc.verifyEqual(R.EF, f(1:5), 'RelTol', 1e-8);
            tc.verifySize(R.U, [numel(R.K(:,1)) 5]);
        end

        function frequencyDomainTendsToStaticSolution(tc)
            % omega -> 0: (K - omega^2 M) u = F gives the static solution
            s = quietSettings('FrequencyDomain_analysis');
            s.fHarmonic = 1e-6;
            [~, Rf] = run(s);
            [~, Rs] = run(quietSettings('Static_analysis'));
            tc.verifyEqual(Rf.u, Rs.u, 'AbsTol', 1e-9*max(abs(Rs.u)));
        end

        function transientAnalysisRuns(tc)
            % Newmark method from rest, sine burst, all DOFs saved
            s = quietSettings('Dynamic_analysis');
            s.p.pxi = 2; s.p.peta = 2;
            s.flag_mRHS = false;
            s.save_var = 1;
            s.t_end = 1e-6;
            s.TimeIntegration.Method = 301;
            s.TimeIntegration.SolverSave = false;
            [model, R] = run(s);
            tc.verifyEqual(R.t_comp(end), 1e-6, 'AbsTol', 1e-18);
            tc.verifySize(R.y_comp, [numel(R.t_comp) model.num_DoF]);
            tc.verifyTrue(all(isfinite(R.y_comp(:))));
            tc.verifyGreaterThan(max(abs(R.y_comp(:))), 0);
            % Constrained DOFs stay at rest
            fixed = setdiff(1:model.num_DoF, R.ActiveDoF);
            tc.verifyEqual(max(abs(R.y_comp(:,fixed)), [], 'all'), 0);
        end

        function resultsAreSaved(tc)
            tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            s = quietSettings('Dynamic_analysis');
            s.p.pxi = 2; s.p.peta = 2;
            s.saveWS = true;
            s.ResultsFolder = fullfile(pwd, 'Results');
            s.TempFolder = fullfile(pwd, 'Temp');
            s.flag_mRHS = false;
            s.t_end = 2e-7;
            s.TimeIntegration.Method = 301;
            [model, R, s] = run(s);
            tc.verifyTrue(isfile(s.ResultsFile));
            S = load(s.ResultsFile);
            tc.verifyEqual(sort(fieldnames(S)), {'model'; 'results'; 'settings'});
            tc.verifyEqual(S.results.y_comp, R.y_comp);
            % y-displacements of the nodes of interest sorted by x: only
            % written if both DOFs of every node of interest are active
            uy = strrep(s.ResultsFile, '_CMM.mat', '_Uy.mat');
            tc.verifyEqual(isfile(uy), size(R.y_comp, 2) == 2*numel(R.NOI));
            % Input of the transient solver (rerun without preprocessing)
            solverFile = fullfile(s.TempFolder, 'Solver_input.mat');
            tc.verifyTrue(isfile(solverFile));
            [~, T, U] = evalc('Solver_Transient(solverFile)');
            tc.verifyEqual(T, R.t_comp);
            tc.verifyEqual(U, R.y_comp);
            tc.verifyEqual(model.name, 'PatchTestLeeCubic_Rajendran_adj');
        end

        function unknownAnalysisTypeIsRejected(tc)
            s = quietSettings('Buckling_analysis');
            tc.verifyError(@() run(s), 'RunSimulation:analysis_type');
        end

        function missingInputFileIsReported(tc)
            p = struct('pxi', 2, 'peta', 2, 'ndist', 2);
            tc.verifyError(@() ReadModel('no_such_file.inp', p, false), ...
                'ReadModel:fileNotFound');
        end
    end
end

function s = quietSettings(analysis)
s = DefaultSettings();
s.analysis_type = analysis;
s.saveWS = false;
s.plot_disp = false;
s.parallel = false;
s.verbose = false;
s.numEigen = 10;
end

function [model, results, settings] = run(settings) %#ok<INUSD>
% RunSimulation without command window output
[~, model, results, settings] = evalc('RunSimulation(settings)');
end
