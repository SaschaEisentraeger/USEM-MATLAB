classdef SolverTransientTest < matlab.unittest.TestCase
    %SOLVERTRANSIENTTEST Solver_Transient as called by call_Dynamic.
    %
    %   The solver reads its input from subroutines/Temp/*.mat relative to
    %   the current folder. The tests write these files for the 2-DOF
    %   oscillator M = I, K = [2 -1; -1 2] (free vibration from
    %   U0 = [1; 0]) to a temporary folder and run the solver there.

    properties (TestParameter)
        % method code, order (higher-order methods), rho_inf (Pade), tolerance
        method = struct( ...
            'Newmark',      {{301, 2, 1,   1e-3}}, ...
            'rhoBathe',     {{306, 2, 1,   1e-3}}, ...
            'KimTwoStage',  {{202, 2, 1,   1e-3}}, ...
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

    methods (TestMethodSetup)
        function useTemporaryFolder(tc)
            tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
            mkdir(fullfile('subroutines', 'Temp'));
        end
    end

    methods (Test)
        function freeVibrationIsAccurate(tc, method)
            writeInput(method{1}, method{2}, method{3}, 1);
            [~, T, U] = evalc('Solver_Transient()');
            Uex = 0.5*[cos(T) + cos(sqrt(3)*T), cos(T) - cos(sqrt(3)*T)];
            tc.verifyEqual(T(end), 2, 'AbsTol', 1e-12);
            tc.verifySize(U, [numel(T) 2]);
            tc.verifyLessThan(max(abs(U(:) - Uex(:))), method{4});
        end

        function allDegreesOfFreedomAreSavedOnRequest(tc)
            % save_var = 1 stores all active DOFs (used to be only one)
            writeInput(301, 2, 1, 1, [], 3);
            [~, ~, U] = evalc('Solver_Transient()');
            tc.verifySize(U, [201 3]);
        end

        function adaptiveHybridMethodUsesElementData(tc)
            % Method 109 works element by element: call_Dynamic.m saves
            % the element matrices to subroutines/Temp/ElemData.mat. Spring
            % chain ground-1-2-ground (assembled: M = I, K = [2 -1; -1 2]).
            writeInput(109, 2, 1, 1);
            % 4 DOFs: 1, 2 active, 3, 4 grounded (constrained)
            d = fullfile('subroutines', 'Temp');
            S = load(fullfile(d, 'Solver_input.mat'));
            S.num_DoF = 4;
            S.U0 = [1; 0; 0; 0];
            S.V0 = zeros(4, 1);
            save(fullfile(d, 'Solver_input.mat'), '-struct', 'S');
            K = [2 -1 -1 0; -1 2 0 -1; -1 0 1 0; 0 -1 0 1];
            M = diag([1 1 0.5 0.5]);
            R = zeros(4);
            F = zeros(4, 1);
            F_dt = struct('Neumann', F, 'K', F, 'C', F, 'M', F); %#ok<NASGU>
            save(fullfile(d, 'K.mat'), 'K');
            save(fullfile(d, 'M.mat'), 'M');
            save(fullfile(d, 'R.mat'), 'R');
            save(fullfile(d, 'F.mat'), 'F');
            save(fullfile(d, 'F_dt.mat'), 'F_dt');
            ElemData.Ke = {[1 -1; -1 1]; [1 -1; -1 1]; [1 -1; -1 1]};
            ElemData.Me = {0.5*eye(2); 0.5*eye(2); 0.5*eye(2)};
            ElemData.Topo = [3 1; 1 2; 2 4];
            ElemData.damped = false(3, 1);
            ElemData.alpha = 0;
            ElemData.beta = 0; %#ok<STRNU>
            save(fullfile(d, 'ElemData.mat'), 'ElemData');
            [~, T, U] = evalc('Solver_Transient()');
            Uex = 0.5*[cos(T) + cos(sqrt(3)*T), cos(T) - cos(sqrt(3)*T)];
            tc.verifySize(U, [numel(T) 4]);
            tc.verifyLessThan(max(abs(U(:,1:2) - Uex), [], 'all'), 2e-3);
        end

        function adaptiveHybridMethodNeedsElementData(tc)
            writeInput(109, 2, 1, 1);
            tc.verifyError(@() evalc('Solver_Transient()'), ...
                'Solver_Transient:method109');
        end
    end
end

function writeInput(code, order, rho, save_var, NOI, n)
% Write the files that call_Dynamic.m saves for Solver_Transient.m
if nargin < 5, NOI = []; end
if nargin < 6, n = 2; end
K = full(spdiags(repmat([-1 2 -1], n, 1), -1:1, n, n));
M = eye(n);
R = zeros(n);
F = zeros(n, 1);
F_dt = struct('Neumann', zeros(n,1), 'K', zeros(n,1), 'C', zeros(n,1), 'M', zeros(n,1));
U0 = zeros(n, 1); U0(1) = 1;
V0 = zeros(n, 1);
num_DoF = n;
ActiveDoF = (1:n)';
t_end = 2;
dtf = 0.015; dte = 0.017;
ExciteParam = t_end; SignalType = 0;
TimeIntegration = struct('Method', code, 'Order', order, 'rho', rho, ...
    'dtImp', 0.01, 'nInterp', 1, 'SolverSave', true);
Topology_Node = [];
flag_mRHS = false; %#ok<*NASGU>
d = fullfile('subroutines', 'Temp');
save(fullfile(d, 'Solver_input.mat'), 'ActiveDoF', 'dte', 'dtf', 'ExciteParam', ...
    'NOI', 'num_DoF', 'save_var', 'SignalType', 't_end', 'TimeIntegration', ...
    'Topology_Node', 'U0', 'V0', 'flag_mRHS');
save(fullfile(d, 'K.mat'), 'K');
save(fullfile(d, 'M.mat'), 'M');
save(fullfile(d, 'R.mat'), 'R');
save(fullfile(d, 'F.mat'), 'F');
save(fullfile(d, 'F_dt.mat'), 'F_dt');
end
