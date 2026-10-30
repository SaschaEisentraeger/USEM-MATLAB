classdef NodalDistributionTest < matlab.unittest.TestCase
    %NODALDISTRIBUTIONTEST Supporting points of the Lagrange polynomials.

    properties (TestParameter)
        p = {1, 2, 3, 4, 5, 8};
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
        function nodesAreAscendingAndSpanReferenceInterval(tc, p)
            x = nodal_distribution(p);
            tc.verifySize(x, [p+1 1]);
            tc.verifyEqual(x(1), -1, 'AbsTol', 1e-14);
            tc.verifyEqual(x(end), 1, 'AbsTol', 1e-14);
            tc.verifyTrue(all(diff(x) > 0), 'Nodes must be strictly ascending.');
        end

        function nodesAreSymmetric(tc, p)
            x = nodal_distribution(p);
            tc.verifyEqual(x, -flipud(x), 'AbsTol', 1e-14);
        end

        function gllNodesAreLegendreDerivativeRoots(tc, p)
            % Interior GLL nodes are the roots of P_p'(x)
            tc.assumeGreaterThan(p, 1, 'No interior nodes for p = 1.');
            x = nodal_distribution(p);
            dP = polyder(legendreCoefficients(p));
            tc.verifyEqual(polyval(dP, x(2:end-1)), zeros(p-1,1), ...
                'AbsTol', 1e-10*max(abs(dP)));
        end

    end
end

function P = legendreCoefficients(n)
% Coefficients (highest degree first) of the Legendre polynomial P_n
P0 = 1; P1 = [1 0];
if n == 0, P = P0; return, end
for k = 1:n-1
    P2 = ((2*k+1)*[P1 0] - k*[0 0 P0])/(k+1);
    P0 = P1; P1 = P2;
end
P = P1;
end
