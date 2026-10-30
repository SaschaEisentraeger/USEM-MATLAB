classdef SolverTransientTest < matlab.unittest.TestCase
    %SOLVERTRANSIENTTEST Solver_Transient as called by DynamicAnalysis.
    %
    %   The solver input is the struct that DynamicAnalysis builds (or the
    %   MAT file it saves). The tests use the 2-DOF oscillator M = I,
    %   K = [2 -1; -1 2] (free vibration from U0 = [1; 0]).

    properties (TestParameter)
        % method code, order (higher-order methods), rho_inf (Pade), tolerance
        method = struct( ...
            'Newmark',      {{301, 2, 1,   1e-3}}, ...
            'rhoBathe',     {{306, 2, 1,   1e-3}}, ...
            'KimTwoStage',  {{202, 2, 1,   1e-3}}, ...
            'SoaresHO',     {{402, 2, 1,   1e-3}}, ...
            'Pade44',       {{403, 4, 1,   1e-6}}, ...
            'PadeDamped',   {{403, 4, 0.5, 1e-3}});
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'code', 'subroutines'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'tests', 'helpers')));
            old = get(groot, 'DefaultFigureVisible');
            set(groot, 'DefaultFigureVisible', 'off');
            tc.addTeardown(@() set(groot, 'DefaultFigureVisible', old));
        end
    end

    methods (Test)
        function freeVibrationIsAccurate(tc, method)
            S = solverInput(method{1}, method{2}, method{3}, 1);
            [~, T, U] = evalc('Solver_Transient(S)');
            tc.verifyEqual(T(end), 2, 'AbsTol', 1e-12);
            tc.verifySize(U, [numel(T) 2]);
            tc.verifyLessThan(max(abs(U(:) - reshape(exactFree(T), [], 1))), method{4});
        end

        function padeStartsFromRest(tc)
            % Zero initial conditions (the usual case in DynamicAnalysis)
            % used to fail in PadeTimeIntegrator (undefined ndof), and a
            % constant load without signal (SignalType 0) was ignored.
            % Step load F = [1; 0.5] from rest.
            S = solverInput(403, 4, 1, 1);
            S.U0 = zeros(2, 1);
            S.F = [1; 0.5];
            [~, T, U] = evalc('Solver_Transient(S)');
            [Phi, W2] = eig(S.K);
            Uex = zeros(numel(T), 2);
            for i = 1:2
                Uex = Uex + (1 - cos(sqrt(W2(i,i))*T)) * (Phi(:,i)*Phi(:,i)'*S.F)'/W2(i,i);
            end
            tc.verifyLessThan(max(abs(U(:) - Uex(:))), 1e-6);
        end

        function inputCanBeReadFromFile(tc)
            % DynamicAnalysis saves the input as Solver_input.mat
            tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            S = solverInput(301, 2, 1, 1);
            save('Solver_input.mat', '-struct', 'S');
            [~, T1, U1] = evalc('Solver_Transient(S)');
            [~, T2, U2] = evalc('Solver_Transient(''Solver_input.mat'')');
            tc.verifyEqual(T2, T1);
            tc.verifyEqual(U2, U1);
        end

        function unknownMethodFallsBackToNewmark(tc)
            S = solverInput(999, 2, 1, 1);
            [~, T, U] = evalc(['tc.verifyWarning(@() Solver_Transient(S), ' ...
                '''Solver_Transient:unknownMethod'')']);
            S.TimeIntegration.Method = 301;
            [~, Tn, Un] = evalc('Solver_Transient(S)');
            tc.verifyEqual(T, Tn);
            tc.verifyEqual(U, Un);
        end

        function everyMethodCodeIsUnique(tc)
            codes = [TimeIntegrationMethods().code];
            tc.verifyEqual(numel(unique(codes)), numel(codes));
        end

        function allDegreesOfFreedomAreSavedOnRequest(tc)
            % save_var = 1 stores all active DOFs (used to be only one)
            S = solverInput(301, 2, 1, 1, [], 3);
            [~, ~, U] = evalc('Solver_Transient(S)');
            tc.verifySize(U, [201 3]);
        end

        function adaptiveHybridMethodUsesElementData(tc)
            % Method 109 works element by element. Spring chain
            % ground-1-2-ground (assembled: M = I, K = [2 -1; -1 2]);
            % DOFs 1, 2 active, 3, 4 grounded (constrained).
            S = solverInput(109, 2, 1, 1);
            S.num_DoF = 4;
            S.U0 = [1; 0; 0; 0];
            S.V0 = zeros(4, 1);
            S.K = [2 -1 -1 0; -1 2 0 -1; -1 0 1 0; 0 -1 0 1];
            S.M = diag([1 1 0.5 0.5]);
            S.R = zeros(4);
            S.F = zeros(4, 1);
            S.F_dt = struct('Neumann', S.F, 'K', S.F, 'C', S.F, 'M', S.F);
            S.ElemData.Ke = {[1 -1; -1 1]; [1 -1; -1 1]; [1 -1; -1 1]};
            S.ElemData.Me = {0.5*eye(2); 0.5*eye(2); 0.5*eye(2)};
            S.ElemData.Topo = [3 1; 1 2; 2 4];
            S.ElemData.damped = false(3, 1);
            S.ElemData.alpha = 0;
            S.ElemData.beta = 0;
            [~, T, U] = evalc('Solver_Transient(S)');
            tc.verifySize(U, [numel(T) 4]);
            tc.verifyLessThan(max(abs(U(:,1:2) - exactFree(T)), [], 'all'), 2e-3);
        end

        function adaptiveHybridMethodNeedsElementData(tc)
            S = solverInput(109, 2, 1, 1);
            tc.verifyError(@() quietSolve(S), 'Solver_Transient:method109');
        end
    end
end

function S = solverInput(code, order, rho, save_var, NOI, n)
% Input of Solver_Transient as built by DynamicAnalysis
if nargin < 5, NOI = []; end
if nargin < 6, n = 2; end
S.K = full(spdiags(repmat([-1 2 -1], n, 1), -1:1, n, n));
S.M = eye(n);
S.R = zeros(n);
S.F = zeros(n, 1);
S.F_dt = struct('Neumann', zeros(n,1), 'K', zeros(n,1), 'C', zeros(n,1), 'M', zeros(n,1));
S.U0 = zeros(n, 1); S.U0(1) = 1;
S.V0 = zeros(n, 1);
S.ActiveDoF = (1:n)';
S.NOI = NOI;
S.num_DoF = n;
S.save_var = save_var;
S.t_end = 2;
S.dtf = 0.015;
S.dte = 0.017;
S.SignalType = 0;
S.ExciteParam = S.t_end;
S.TimeIntegration = struct('Method', code, 'Order', order, 'rho', rho, ...
    'dtImp', 0.01, 'nInterp', 1, 'SolverSave', true);
S.flag_mRHS = false;
end

function quietSolve(S) %#ok<INUSD>
% Solver_Transient without command window output
evalc('Solver_Transient(S);');
end

function U = exactFree(T)
% U0 = [1;0], V0 = 0: modes [1;1]/sqrt2 (w=1) and [1;-1]/sqrt2 (w=sqrt3)
U = 0.5*[cos(T) + cos(sqrt(3)*T), cos(T) - cos(sqrt(3)*T)];
end
