classdef NodalDistributionTest < matlab.unittest.TestCase
    %NODALDISTRIBUTIONTEST Supporting points of the Lagrange polynomials.

    properties (TestParameter)
        p = {1, 2, 3, 4, 5, 8};
        ndist = struct('equidistant', 1, 'GLL', 2, 'GLC', 3);
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
        function nodesAreAscendingAndSpanReferenceInterval(tc, p, ndist)
            x = nodal_distribution(p, ndist);
            tc.verifySize(x, [p+1 1]);
            tc.verifyEqual(x(1), -1, 'AbsTol', 1e-14);
            tc.verifyEqual(x(end), 1, 'AbsTol', 1e-14);
            tc.verifyTrue(all(diff(x) > 0), 'Nodes must be strictly ascending.');
        end

        function nodesAreSymmetric(tc, p, ndist)
            x = nodal_distribution(p, ndist);
            tc.verifyEqual(x, -flipud(x), 'AbsTol', 1e-14);
        end

        function equidistantNodesAreEquallySpaced(tc, p)
            x = nodal_distribution(p, 1);
            tc.verifyEqual(diff(x), 2/p*ones(p,1), 'AbsTol', 1e-14);
        end

        function gllNodesAreLegendreDerivativeRoots(tc, p)
            % Interior GLL nodes are the roots of P_p'(x)
            tc.assumeGreaterThan(p, 1, 'No interior nodes for p = 1.');
            x = nodal_distribution(p, 2);
            dP = polyder(legendreCoefficients(p));
            tc.verifyEqual(polyval(dP, x(2:end-1)), zeros(p-1,1), ...
                'AbsTol', 1e-10*max(abs(dP)));
        end

        function glcNodesAreChebyshevExtrema(tc, p)
            x = nodal_distribution(p, 3);
            tc.verifyEqual(x, -cos(pi*(0:p)'/p), 'AbsTol', 1e-14);
        end
    end

        function unknownDistributionIsRejected(tc)
            tc.verifyError(@() nodal_distribution(3, 4), 'nodal_distribution:ndist');
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
