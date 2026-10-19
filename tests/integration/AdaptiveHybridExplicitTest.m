classdef AdaptiveHybridExplicitTest < matlab.unittest.TestCase
    %ADAPTIVEHYBRIDEXPLICITTEST Adaptive hybrid explicit method (Soares
    %2020, time integration method 109).
    %
    %   - Spring chain ground-1-2-ground written as three 2-DOF "elements"
    %     (assembled: M = I, K = [2 -1; -1 2]): accuracy, convergence
    %     order, mass-proportional damping.
    %   - Spectral element cantilever (lumped mass): stable and accurate
    %     at twice the critical time step of the central difference
    %     method, which diverges there.

    properties
        chain       % element data of the spring chain
        beam        % spectral element model
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

        function buildModels(tc)
            ED.Ke = {[1 -1; -1 1]; [1 -1; -1 1]; [1 -1; -1 1]};
            ED.Me = {0.5*eye(2); 0.5*eye(2); 0.5*eye(2)};
            ED.Topo = [3 1; 1 2; 2 4];      % DOFs 3 and 4: ground
            ED.damped = false(3, 1);
            ED.alpha = 0;
            ED.beta = 0;
            tc.chain = ED;

            % Cantilever [0,8] x [0,1], 8 x 2 elements, p = 3, GLL nodal
            % quadrature (lumped mass), plane stress, tip load
            [nodes, topo] = usemTestMesh(8, 2, 8, 1, 1);
            m = usemBuildModel(nodes, topo, 3);
            Cmat = HookeMatrix2d(1, 0.3, 1);
            [~, K, M, R, dte, ~, ED] = evalc(['AssemblyFEM(m.Topology_Node, ' ...
                'm.p, m.xCoor, m.yCoor, m.num_Elem, m.num_DoF, {Cmat}, {1}, 1, 2, ' ...
                '{(1:m.num_Elem)''}, 0, 0, {0}, 2)']);
            fixed = find(abs(m.xCoor) < 1e-12);
            act = setdiff((1:m.num_DoF)', [2*fixed-1; 2*fixed]);
            tip = find(abs(m.xCoor - 8) < 1e-12 & abs(m.yCoor - 1) < 1e-12);
            F = zeros(m.num_DoF, 1);
            F(2*tip) = -1e-3;
            tc.beam = struct('K', K(act,act), 'M', M(act,act), 'R', R(act,act), ...
                'F', F(act), 'act', act, 'num_DoF', m.num_DoF, 'dte', dte, ...
                'dofTip', find(act == 2*tip), 'ED', ED);
        end
    end

    methods (Test)
        function freeVibrationConvergesWithSecondOrder(tc)
            dt = [0.02 0.01];
            e = zeros(size(dt));
            for k = 1:2
                [T, U] = runChain(tc, dt(k), tc.chain, zeros(2));
                Uex = 0.5*[cos(T) + cos(sqrt(3)*T), cos(T) - cos(sqrt(3)*T)];
                e(k) = max(abs(U(:,1:2) - Uex), [], 'all');
            end
            tc.verifyLessThan(e(2), 1e-3);
            tc.verifyGreaterThan(log2(e(1)/e(2)), 1.8, mat2str(e, 3));
        end

        function massProportionalDamping(tc)
            % C = a*M: modal solution with zeta_i = a/(2*omega_i)
            a = 0.2;
            ED = tc.chain;
            ED.damped = true(3, 1);
            ED.alpha = a;
            [T, U] = runChain(tc, 0.01, ED, a*eye(2));
            [Phi, W2] = eig([2 -1; -1 2]);
            Uex = zeros(numel(T), 2);
            for i = 1:2
                wd = sqrt(W2(i,i) - a^2/4);
                q = exp(-a*T/2).*(cos(wd*T) + a/(2*wd)*sin(wd*T));
                Uex = Uex + q*(Phi(:,i)*Phi(:,i)'*[1; 0])';
            end
            tc.verifyLessThan(max(abs(U(:,1:2) - Uex), [], 'all'), 1e-3);
        end

        function stableAtTwiceTheCentralDifferenceLimit(tc)
            b = tc.beam;
            tsim = 40;
            n = numel(b.act);
            z = zeros(n, 1);
            % reference: trapezoidal rule (Newmark) with a small time step
            [~, Tr, Ur] = evalc(['NewmarkMethod(b.K, b.M, b.R, b.F, z, z, z, ' ...
                'b.num_DoF, b.dofTip, b.act, tsim, 0, [1 1 0.1]*b.dte, tsim)']);
            % adaptive hybrid method: dt = 2*dte (dte = 0.99*2/omega_max)
            [~, T, U] = evalc(['AdaptiveHybridExplicitMethod(b.K, b.M, b.R, ' ...
                'b.F, z, z, z, b.num_DoF, b.dofTip, b.act, tsim, 0, ' ...
                '[1e9 b.dte 1e9], tsim, b.ED)']);
            tc.verifyEqual(T(2) - T(1), 2*b.dte, 'RelTol', 0.01);
            ref = interp1(Tr, Ur, T);
            tc.verifyLessThan(max(abs(U - ref))/max(abs(ref)), 0.02);
            % central difference method with the same time step diverges
            [~, ~, Uc] = evalc(['CentralDifferenceMethod(b.K, b.M, b.R, b.F, ' ...
                'z, z, z, b.num_DoF, b.dofTip, b.act, tsim, 0, ' ...
                '[1e9 2*b.dte 1e9], tsim, false)']);
            tc.verifyGreaterThan(max(abs(Uc)), 1e3*max(abs(ref)));
        end

        function lumpedMassIsDiagonal(tc)
            % GLL nodal quadrature gives an exactly diagonal mass matrix,
            % so the effective matrix E needs no solver
            tc.verifyTrue(isdiag(tc.beam.M));
        end

        function missingElementMassMatricesRaiseError(tc)
            ED = tc.chain;
            ED.Me = {};
            tc.verifyError(@() AdaptiveHybridExplicitMethod([2 -1; -1 2], ...
                eye(2), zeros(2), [0;0], zeros(2,1), [1;0], [0;0], 4, [1 2], ...
                [1;2], 1, 0, 0.01*[1 1 1], 1, ED), ...
                'AdaptiveHybridExplicitMethod:mass');
        end
    end
end

function [T, U] = runChain(~, dt, ED, C)
% dtmax = [dtf, dte, dtImp]; the method uses min(dtmax) after doubling dte
[~, T, U] = evalc(['AdaptiveHybridExplicitMethod([2 -1; -1 2], eye(2), C, ' ...
    '[0;0], zeros(2,1), [1;0], [0;0], 4, [1 2], [1;2], 2, 0, ' ...
    'dt*[1 1 1], 20, ED)']);
end
