classdef QuadratureTest < matlab.unittest.TestCase
    %QUADRATURETEST Gauss-Legendre (lgwt) and Gauss-Lobatto-Legendre
    %(lglnodes) rules.

    properties (TestParameter)
        n = {1, 2, 3, 4, 6, 9, 12};
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
        function gaussLegendreIsExactForDegree2nMinus1(tc, n)
            [x, w] = lgwt(n, -1, 1);
            tc.verifySize(x, [n 1]);
            for k = 0:2*n-1
                exact = (1 - (-1)^(k+1))/(k+1);
                tc.verifyEqual(w'*x.^k, exact, 'AbsTol', 1e-13, ...
                    sprintf('n = %d, monomial degree %d', n, k));
            end
        end

        function gaussLegendreMapsToArbitraryInterval(tc, n)
            a = 0.5; b = 3.25;
            [x, w] = lgwt(n, a, b);
            tc.verifyTrue(all(x > a & x < b));
            tc.verifyEqual(sum(w), b - a, 'AbsTol', 1e-13);
            k = 2*n - 1;
            tc.verifyEqual(w'*x.^k, (b^(k+1) - a^(k+1))/(k+1), 'RelTol', 1e-12);
        end

        function gaussLegendrePointsAreSymmetric(tc, n)
            [x, w] = lgwt(n, -1, 1);
            tc.verifyEqual(sort(x), sort(-x), 'AbsTol', 1e-14);
            tc.verifyTrue(all(w > 0));
        end

        function gaussLobattoIsExactForDegree2nMinus1(tc, n)
            % lglnodes(N) returns N+1 points that integrate degree 2N-1
            [x, w] = lglnodes(n);
            tc.verifySize(x, [n+1 1]);
            tc.verifyEqual(sort(abs(x([1 end]))), [1; 1], 'AbsTol', 1e-15);
            for k = 0:2*n-1
                exact = (1 - (-1)^(k+1))/(k+1);
                tc.verifyEqual(w'*x.^k, exact, 'AbsTol', 1e-13, ...
                    sprintf('N = %d, monomial degree %d', n, k));
            end
        end
    end
end
