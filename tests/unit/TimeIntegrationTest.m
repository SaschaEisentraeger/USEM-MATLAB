classdef TimeIntegrationTest < matlab.unittest.TestCase
    %TIMEINTEGRATIONTEST Direct time integration schemes on a 2-DOF
    %oscillator with a known (modal) solution.
    %
    %   M = I, K = [2 -1; -1 2], no damping. Free vibration from
    %   U0 = [1; 0] and the step response to a constant load are compared
    %   with the exact solution. All schemes share the calling convention
    %       [T, U, ...] = Method(K, M, C, F, F_dt, U0, V0, NumDoF, DOFsave,
    %                            ActiveDoFs, ExciteParam, SignalType,
    %                            dtmax, tsim, ...)

    properties (TestParameter)
        % handle, extra arguments, observed order of accuracy on this
        % problem (Mirbagheri: first order observed, see
        % docs/known-issues.md; rho-Bathe: see the KnownIssue test)
        scheme = struct( ...
            'Newmark',            {{@NewmarkMethod,                {}, 2}}, ...
            'HHTAlpha',           {{@HHTAlphaMethod,               {}, 2}}, ...
            'GeneralizedAlpha',   {{@GeneralizedAlphaMethod,       {}, 2}}, ...
            'BatheBaig',          {{@BatheBaigMethod,              {}, 2}}, ...
            'b1b2Bathe',          {{@b1b2_BatheMethod,             {}, 2}}, ...
            'KimImplicit',        {{@KimMethodImplicit,            {}, 2}}, ...
            'KimCollocation',     {{@KimMethodImplicitCollocation, {}, 2}}, ...
            'CentralDifference',  {{@CentralDifferenceMethod,      {false}, 2}}, ...
            'NohBathe',           {{@NohBatheMethod,               {}, 2}}, ...
            'SoaresExplicit',     {{@SoaresMethodExplicit,         {}, 2}}, ...
            'KimExplicit',        {{@KimMethodExplicit,            {}, 2}}, ...
            'Mirbagheri',         {{@MirbagheriMethod,             {}, 1}}, ...
            'ZhangXing1',         {{@ZhangXingMethod1,             {}, 2}}, ...
            'ZhangXing2',         {{@ZhangXingMethod2,             {}, 2}}, ...
            'ExplicitGenAlpha',   {{@ExplicitGeneralizedAlphaMethod, {}, 2}}, ...
            'KimExplicitHO3',     {{@KimExplicit_ho,               {3}, 3}}, ...
            'KimLagrangeHO3',     {{@KimLagrange_ho,               {3}, 3}});
    end

    properties (Constant)
        K = [2 -1; -1 2];
        M = eye(2);
        tsim = 2;
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'code', 'subroutines'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'tests', 'helpers')));
        end

        function hideFigures(tc)
            % Several schemes open a waitbar
            tc.applyFixture(matlab.unittest.fixtures.SuppressedWarningsFixture( ...
                'MATLAB:hg:AutoSoftwareOpenGL'));
            old = get(groot, 'DefaultFigureVisible');
            set(groot, 'DefaultFigureVisible', 'off');
            tc.addTeardown(@() set(groot, 'DefaultFigureVisible', old));
        end
    end

    methods (Test)
        function freeVibrationIsAccurate(tc, scheme)
            dt = 0.01;
            [T, U] = integrate(tc, scheme, dt, [0; 0], [1; 0]);
            Uex = exactFree(T);
            tc.verifyEqual(T(end), tc.tsim, 'AbsTol', 1e-12);
            tc.verifyLessThan(max(abs(U(:) - Uex(:))), 1e-3);
        end

        function convergesWithExpectedOrder(tc, scheme)
            dt = [0.02 0.01];
            e = zeros(size(dt));
            for k = 1:2
                [T, U] = integrate(tc, scheme, dt(k), [0; 0], [1; 0]);
                Uex = exactFree(T);
                e(k) = max(abs(U(:) - Uex(:)));
            end
            rate = log2(e(1)/e(2));
            tc.verifyGreaterThan(rate, scheme{3} - 0.3, ...
                sprintf('errors %s', mat2str(e, 3)));
        end

        function stepLoadResponseIsAccurate(tc, scheme)
            % u(t) = sum_i phi_i phi_i' F (1 - cos(w_i t)) / w_i^2
            F = [1; 0.5];
            [T, U] = integrate(tc, scheme, 0.01, F, [0; 0]);
            [Phi, W2] = eig(tc.K);
            Uex = zeros(numel(T), 2);
            for i = 1:2
                Uex = Uex + (1 - cos(sqrt(W2(i,i))*T)) * (Phi(:,i)*Phi(:,i)'*F)'/W2(i,i);
            end
            tc.verifyLessThan(max(abs(U(:) - Uex(:))), 2e-3);
        end

        function multipleRightHandSidesMatchSingleRuns(tc)
            % Default route of Solver_Transient (flag_mRHS = true)
            F = [1 0; 0.5 1];
            dtmax = 0.01*[1.5 1.7 1];
            for solver = {@CentralDifferenceMethod_mRHS, @HHTAlphaMethod_mRHS}
                if isequal(solver{1}, @CentralDifferenceMethod_mRHS)
                    extra = {false};
                else
                    extra = {};
                end
                [~, ~, U] = evalc(['solver{1}(tc.K, tc.M, zeros(2), F, zeros(2,2), ' ...
                    'zeros(2,1), zeros(2,1), 2, [1 2], [1;2], tc.tsim, 0, dtmax, tc.tsim, extra{:})']);
                for lc = 1:2
                    if isequal(solver{1}, @CentralDifferenceMethod_mRHS)
                        ref = {@CentralDifferenceMethod, false};
                    else
                        ref = {@HHTAlphaMethod};
                    end
                    [~, ~, Ur] = evalc(['ref{1}(tc.K, tc.M, zeros(2), F(:,lc), zeros(2,1), ' ...
                        'zeros(2,1), zeros(2,1), 2, [1 2], [1;2], tc.tsim, 0, dtmax, ' ...
                        'tc.tsim, ref{2:end})']);
                    tc.verifyEqual(squeeze(U(:,:,lc)), Ur, 'AbsTol', 1e-12, func2str(solver{1}));
                end
            end
        end
    end

    methods (Test, TestTags = {'KnownIssue'})
        function rhoBatheMethodRuns(tc)
            % rho_BatheMethod names its 5th argument Fdt but uses F_dt
            % (method 306, see docs/known-issues.md).
            [T, U] = integrate(tc, {@rho_BatheMethod, {}, 2}, 0.01, [0; 0], [1; 0]);
            tc.verifyLessThan(max(abs(U(:) - reshape(exactFree(T), [], 1))), 1e-3);
        end

        function mirbagheriMethodIsSecondOrder(tc)
            % Observed convergence rate ~1 on the free-vibration problem.
            dt = [0.02 0.01];
            e = zeros(size(dt));
            for k = 1:2
                [T, U] = integrate(tc, {@MirbagheriMethod, {}, 2}, dt(k), [0; 0], [1; 0]);
                e(k) = max(abs(U(:) - reshape(exactFree(T), [], 1)));
            end
            tc.verifyGreaterThan(log2(e(1)/e(2)), 1.7);
        end
    end
end

function [T, U] = integrate(tc, scheme, dt, F, U0)
% dtmax = [dtf, dte, dtImp]: implicit schemes use dtImp (three distinct
% values), explicit schemes use min(dtmax) after scaling dte by their
% stability factor.
dtmax = dt*[1.5 1.7 1];
fh = scheme{1};
extra = scheme{2};
[~, T, U] = evalc(['fh(tc.K, tc.M, zeros(2), F, zeros(2,1), U0, zeros(2,1), 2, ' ...
                   '[1 2], [1;2], tc.tsim, 0, dtmax, tc.tsim, extra{:})']);
end

function U = exactFree(T)
% U0 = [1;0], V0 = 0: modes [1;1]/sqrt2 (w=1) and [1;-1]/sqrt2 (w=sqrt3)
U = 0.5*[cos(T) + cos(sqrt(3)*T), cos(T) - cos(sqrt(3)*T)];
end
