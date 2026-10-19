classdef ShapeFunctionsTest < matlab.unittest.TestCase
    %SHAPEFUNCTIONSTEST Parametric (Lagrange) shape functions in 1D and 2D.

    properties (TestParameter)
        p = {1, 2, 3, 4, 6};
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
        %% 1D
        function lagrange1DHasKroneckerProperty(tc, p, ndist)
            xn = nodal_distribution(p, ndist);
            N = ShapeFunctions1D_SEM(p, xn, ndist);
            % Order of the 1D shape functions: [xi=-1, xi=+1, interior]
            order = [1, p+1, 2:p];
            tc.verifyEqual(N(:, order), eye(p+1), 'AbsTol', 1e-12);
        end

        function lagrange1DReproducesPolynomials(tc, p, ndist)
            xn = nodal_distribution(p, ndist);
            xq = linspace(-1, 1, 23)';
            [N, dN] = ShapeFunctions1D_SEM(p, xq, ndist);
            order = [1, p+1, 2:p];
            for k = 0:p
                f  = xn(order).^k;
                tc.verifyEqual(N'*f, xq.^k, 'AbsTol', 1e-11);
                df = k*xq.^max(k-1, 0);
                tc.verifyEqual(dN'*f, df, 'AbsTol', 1e-9);
            end
        end

        %% 2D spectral element
        function sem2DHasKroneckerProperty(tc, p, ndist)
            P.pxi = p; P.peta = p; P.ndist = ndist;
            [xi, eta] = referenceNodes(P);
            N = ShapeFunctions2D_SEM(xi, eta, P);
            tc.verifyEqual(N, eye((p+1)^2), 'AbsTol', 1e-12);
        end

        function sem2DIsPartitionOfUnity(tc, p, ndist)
            P.pxi = p; P.peta = p; P.ndist = ndist;
            [xq, eq] = randomPoints(40);
            [N, dNxi, dNeta] = ShapeFunctions2D_SEM(xq, eq, P);
            tc.verifyEqual(sum(N, 1), ones(1, 40), 'AbsTol', 1e-12);
            tc.verifyEqual(sum(dNxi, 1), zeros(1, 40), 'AbsTol', 1e-10);
            tc.verifyEqual(sum(dNeta, 1), zeros(1, 40), 'AbsTol', 1e-10);
        end

        function sem2DReproducesTensorProductPolynomials(tc, p)
            P.pxi = p; P.peta = p; P.ndist = 2;
            [xi, eta] = referenceNodes(P);
            [xq, eq] = randomPoints(30);
            [N, dNxi, dNeta] = ShapeFunctions2D_SEM(xq, eq, P);
            for a = 0:p
                for b = 0:p
                    f = xi.^a .* eta.^b;
                    tc.verifyEqual(N'*f, xq.^a .* eq.^b, 'AbsTol', 1e-10);
                    tc.verifyEqual(dNxi'*f, a*xq.^max(a-1,0) .* eq.^b, 'AbsTol', 1e-8);
                    tc.verifyEqual(dNeta'*f, xq.^a .* b.*eq.^max(b-1,0), 'AbsTol', 1e-8);
                end
            end
        end

        function sem2DSupportsAnisotropicDegrees(tc)
            P.pxi = 3; P.peta = 1; P.ndist = 2;
            [xi, eta] = referenceNodes(P);
            N = ShapeFunctions2D_SEM(xi, eta, P);
            tc.verifyEqual(N, eye(8), 'AbsTol', 1e-12);
        end

        %% Classical 4-node and 8-node shape functions
        function bilinearAndSerendipityAreInterpolatory(tc)
            xc = [-1 1 1 -1]'; ec = [-1 -1 1 1]';
            N4 = ShapeFunctions2D_4nodes(xc, ec);
            tc.verifyEqual(N4, eye(4), 'AbsTol', 1e-15);
            x8 = [xc; 0; 1; 0; -1]; e8 = [ec; -1; 0; 1; 0];
            N8 = ShapeFunctions2D_8nodes(x8, e8);
            tc.verifyEqual(N8, eye(8), 'AbsTol', 1e-15);
        end

        function serendipityDerivativesMatchFiniteDifferences(tc)
            [xq, eq] = randomPoints(10);
            [~, dNxi, dNeta] = ShapeFunctions2D_8nodes(xq, eq);
            h = 1e-6;
            Np = ShapeFunctions2D_8nodes(xq+h, eq); Nm = ShapeFunctions2D_8nodes(xq-h, eq);
            tc.verifyEqual(dNxi, (Np-Nm)/(2*h), 'AbsTol', 1e-8);
            Np = ShapeFunctions2D_8nodes(xq, eq+h); Nm = ShapeFunctions2D_8nodes(xq, eq-h);
            tc.verifyEqual(dNeta, (Np-Nm)/(2*h), 'AbsTol', 1e-8);
        end

        function bilinearDerivativesMatchFiniteDifferences(tc)
            [xq, eq] = randomPoints(10);
            [~, dNxi, dNeta] = ShapeFunctions2D_4nodes(xq, eq);
            h = 1e-6;
            Np = ShapeFunctions2D_4nodes(xq+h, eq); Nm = ShapeFunctions2D_4nodes(xq-h, eq);
            tc.verifyEqual(dNxi, (Np-Nm)/(2*h), 'AbsTol', 1e-8);
            Np = ShapeFunctions2D_4nodes(xq, eq+h); Nm = ShapeFunctions2D_4nodes(xq, eq-h);
            tc.verifyEqual(dNeta, (Np-Nm)/(2*h), 'AbsTol', 1e-8);
        end
    end
end

function [xi, eta] = referenceNodes(P)
% Local coordinates of the element nodes in the order used by
% ShapeFunctions2D_SEM: corners, edges 1-4, interior (row by row)
x = nodal_distribution(P.pxi, P.ndist);
e = nodal_distribution(P.peta, P.ndist);
xi  = [-1; 1; 1; -1];
eta = [-1; -1; 1; 1];
xi  = [xi;  x(2:end-1);             ones(P.peta-1,1); x(2:end-1);             -ones(P.peta-1,1)];
eta = [eta; -ones(P.pxi-1,1);        e(2:end-1);       ones(P.pxi-1,1);        e(2:end-1)];
[XI, ETA] = meshgrid(x(2:end-1), e(2:end-1));
XI = XI'; ETA = ETA';
xi  = [xi;  XI(:)];
eta = [eta; ETA(:)];
end

function [xq, eq] = randomPoints(n)
s = RandStream('mt19937ar', 'Seed', 42);
xq = 2*rand(s, n, 1) - 1;
eq = 2*rand(s, n, 1) - 1;
end
