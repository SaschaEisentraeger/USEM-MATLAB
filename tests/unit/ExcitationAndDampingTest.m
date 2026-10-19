classdef ExcitationAndDampingTest < matlab.unittest.TestCase
    %EXCITATIONANDDAMPINGTEST Time signals (ExcitationSignal,
    %ExcitationSignalDerivative) and Rayleigh damping parameters.

    properties (TestParameter)
        % SignalType and ExciteParam as set in Start_Sim.m
        signal = struct( ...
            'sineBurst',  {{1, [5e3, 3]}}, ...
            'gaussBurst', {{2, 5e3}}, ...
            'chirp',      {{3, [1e3, 5e3, 2e-3]}}, ...
            'harmonic',   {{4, [1, 0.5, 2e3, 3e3]}}, ...
            'ricker',     {{6, 5e3}}, ...
            'gaussPulses',{{7, [1, 1e-3, 1e-4, 1.5e-3, 3]}});
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'code', 'subroutines'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'tests', 'helpers')));
        end
    end

    methods (Test)
        function zeroSignalIsZero(tc)
            t = linspace(0, 1, 11);
            tc.verifyEqual(ExcitationSignal(0, 1, t), zeros(1, 11));
            tc.verifyEqual(ExcitationSignalDerivative(0, 1, t, 2), zeros(1, 11));
        end

        function signalIsRowVector(tc, signal)
            t = linspace(0, 4e-3, 50)';
            x = ExcitationSignal(signal{1}, signal{2}, t);
            tc.verifySize(x, [1 50]);
            tc.verifyTrue(all(isfinite(x)));
        end

        function derivativesMatchFiniteDifferences(tc, signal)
            % interior points away from kinks of windowed signals
            t = linspace(1e-4, 0.55e-3, 7);
            h = 1e-9;
            s = @(tt) ExcitationSignal(signal{1}, signal{2}, tt);
            d1 = ExcitationSignalDerivative(signal{1}, signal{2}, t, 1);
            d2 = ExcitationSignalDerivative(signal{1}, signal{2}, t, 2);
            fd1 = (s(t+h) - s(t-h))/(2*h);
            fd2 = (s(t+1e3*h) - 2*s(t) + s(t-1e3*h))/(1e3*h)^2;
            scale1 = max(abs(fd1)) + 1;
            scale2 = max(abs(fd2)) + 1;
            tc.verifyEqual(d1(:), fd1(:), 'AbsTol', 1e-5*scale1);
            tc.verifyEqual(d2(:), fd2(:), 'AbsTol', 1e-3*scale2);
        end

        function sineBurstVanishesOutsideWindow(tc)
            f = 5e3; n = 3;
            t = [0, n/f, 1.01*n/f, 2*n/f];
            x = ExcitationSignal(1, [f n], t);
            tc.verifyEqual(x, zeros(1, 4), 'AbsTol', 1e-12);
            % x(t) = sin(2*pi*f*t) * sin(pi*f*t/n)^2 (Hann window)
            tq = n/(2*f) + 0.25/f;
            x = ExcitationSignal(1, [f n], tq);
            tc.verifyEqual(x, sin(2*pi*f*tq)*sin(pi*f*tq/n)^2, 'AbsTol', 1e-12);
            tc.verifyEqual(abs(x), cos(pi/(4*n))^2, 'AbsTol', 1e-12);
        end

        function piecewiseLinearInterpolatesTable(tc)
            TS = [0 0; 1 1; 3 1; 4 0];
            t = [0 0.5 1 2 3 3.5 4];
            x = ExcitationSignal(5, TS, t);
            tc.verifyEqual(x(:), [0 0.5 1 1 1 0.5 0]', 'AbsTol', 1e-14);
        end

        function rickerWaveletPeaksAtDelay(tc)
            f = 5e3;
            tc.verifyEqual(ExcitationSignal(6, f, 1/f), 1, 'AbsTol', 1e-14);
        end

        function rayleighParametersReproduceTargetDamping(tc)
            w = [2*pi*100, 2*pi*1000];
            z = [0.02, 0.05];
            [a, b] = RayleighDamping(w, z, false);
            zeta = a./(2*w) + b*w/2;
            tc.verifyEqual(zeta, z, 'RelTol', 1e-5);
        end

        function rayleighWithEqualRatios(tc)
            w = [10, 50];
            [a, b] = RayleighDamping(w, 0.03, false);
            tc.verifyEqual(a./(2*w) + b*w/2, [0.03 0.03], 'RelTol', 1e-5);
        end

        function noDampingForZeroFrequencies(tc)
            [a, b] = RayleighDamping([0 0], [0 0], false);
            tc.verifyEqual([a b], [0 0]);
        end
    end
end
