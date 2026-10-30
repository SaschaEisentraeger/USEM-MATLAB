classdef DriverTest < matlab.unittest.TestCase
    %DRIVERTEST Building blocks of the program flow (RunSimulation) and the
    %shared helpers of the time integration methods.

    properties (Constant)
        E  = 70e9;
        nu = 0.3;
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'code', 'subroutines'), 'IncludingSubfolders', true));
        end
    end

    methods (Test)
        function isotropicMaterialGivesPlaneMatrices(tc)
            props = cell(2, 5);
            props(:, 1) = {2700; 7800};
            props(:, 3) = {[tc.E, tc.nu]; [2*tc.E, 0.25]};
            props(:, 5) = {[]; 0.02};
            for antype = 1:2
                [C, rho, zeta] = MaterialMatrices(props, antype);
                tc.verifyEqual(C{1}, HookeMatrix2d(tc.E, tc.nu, antype), 'RelTol', 1e-12);
                tc.verifyEqual(C{2}, HookeMatrix2d(2*tc.E, 0.25, antype), 'RelTol', 1e-12);
                tc.verifyEqual(rho, {2700; 7800});
                tc.verifyEqual(zeta, {[]; 0.02});
            end
        end

        function orthotropicConstantsAreMappedToVoigtOrder(tc)
            % ABAQUS order D1111 D1122 D2222 D1133 D2233 D3333 D1212 D1313
            % D2323; plane strain keeps the 11, 22 and 12 terms
            D = (1:9)';
            props = {1, [], D, [], []};
            C = MaterialMatrices(props, 2);
            tc.verifyEqual(C{1}, [1 2 0; 2 3 0; 0 0 7]);
        end

        function unknownPlaneStateIsRejected(tc)
            props = {1, [], [tc.E, tc.nu], [], []};
            tc.verifyError(@() MaterialMatrices(props, 3), 'MaterialMatrices:antype');
        end

        function alignedElementsKeepTheirDegrees(tc)
            % Shared edge 2-5: local edge 2 of both elements (eta direction)
            T = [1 2 5 4; 2 3 6 5];
            p = AlignPolynomialDegrees(T, degrees(2, 4), false);
            tc.verifyEqual([p.pxi p.peta], [2 4]);
        end

        function rotatedElementGetsEqualDegrees(tc)
            % Element 2 numbered from node 5: the shared edge is its local
            % edge 1 (xi direction) but local edge 2 of element 1 (eta)
            T = [1 2 5 4; 5 2 3 6];
            p = AlignPolynomialDegrees(T, degrees(2, 4), false);
            tc.verifyEqual([p.pxi p.peta], [2 2]);
        end

        function everySignalHasDefaultParameters(tc)
            t = linspace(0, 1e-3, 11);
            for type = 0:7
                P = DefaultExcitationParameters(type, 0.2);
                s = ExcitationSignal(type, P, t);
                tc.verifySize(s, size(t), sprintf('signal %d', type));
                tc.verifyTrue(all(isfinite(s)), sprintf('signal %d', type));
            end
            tc.verifyError(@() DefaultExcitationParameters(8, 1), ...
                'DefaultExcitationParameters:SignalType');
        end

        function implicitTimeStepIsThePrescribedOne(tc)
            % [dtf, dte, dtImp] all different: dtImp; tsim a multiple of dt
            [~, dt, n] = evalc('TimeStepSize([0.3 0.2 0.1], 1)');
            tc.verifyEqual([dt n], [0.1 10], 'AbsTol', 1e-15);
            [~, dt, n] = evalc('TimeStepSize([0.3 0.3 0.3], 1)');
            tc.verifyEqual(n, 4);
            tc.verifyEqual(dt, 0.25, 'AbsTol', 1e-15);
        end

        function explicitTimeStepIsScaledStabilityLimit(tc)
            % min(dtf, factor*dte, dtImp), reduced to divide tsim
            [~, dt, n] = evalc('TimeStepSize([0.3 0.1 0.25], 1, 2)');
            tc.verifyEqual(n, 5);
            tc.verifyEqual(dt, 0.2, 'AbsTol', 1e-15);
        end

        function loadAtTimeCombinesAllParts(tc)
            Fdt = struct('Neumann', [1; 0], 'K', [0; 2], 'C', [3; 0], 'M', [0; 4]);
            tc.verifyEqual(LoadAtTime([1; 1], Fdt, 2, 3, 5), [1+2+9; 1+4+20]);
        end

        function initialAccelerationSolvesEquationOfMotion(tc)
            K = [2 -1; -1 2]; C = 0.1*eye(2); M = diag([2 4]);
            Fdt = struct('Neumann', [0; 0], 'K', [0; 0], 'C', [0; 0], 'M', [0; 0]);
            a0 = InitialAcceleration(@(r) M\r, [1; 2], Fdt, K, C, [0.1; 0], [0; 1], 0, 1);
            tc.verifyEqual(M*a0, [1; 2] - C*[0; 1] - K*[0.1; 0], 'AbsTol', 1e-15);
        end

        function historiesAreExpandedOnlyIfComplete(tc)
            u = [1 2; 3 4];
            U = ExpandToAllDoFs(4, [1 3], u);
            tc.verifyEqual(U, [1 0 2 0; 3 0 4 0]);
            U = ExpandToAllDoFs(4, [1 2 3], u);
            tc.verifyEqual(U, u);
        end

        function triangularStorageIsRestored(tc)
            K = triu(magic(3)); M = eye(3); C = zeros(3);
            [K2, ~, ~, Fdt] = PrepareSystem(K, M, C, [1; 2; 3], 3);
            tc.verifyEqual(K2, K + K');
            tc.verifyEqual(Fdt.Neumann, [1; 2; 3]);
            tc.verifyEqual(Fdt.M, zeros(3, 1));
        end
    end
end

function p = degrees(pxi, peta)
p = struct('pxi', pxi, 'peta', peta, 'ndist', 2);
end
