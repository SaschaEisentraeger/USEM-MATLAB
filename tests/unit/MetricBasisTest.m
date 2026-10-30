classdef MetricBasisTest < matlab.unittest.TestCase
    %METRICBASISTEST Polynomial basis of the generalized Vandermonde matrix
    %of the metric shape functions (PolynomialBasis1D, p.metricBasis in
    %ShapeFunctions2D_Metric).
    %
    %   All bases span Q_p in skew coordinates, so the metric shape
    %   functions must not depend on the basis; the bases differ only in
    %   the conditioning of the Vandermonde matrix
    %   (docs/chapters/metric-basis.tex).

    properties (TestParameter)
        basis = {'monomial', 'legendre', 'chebyshev', 'lagrange'};
        p = {1, 2, 3, 4, 6};
        pHigh = {8, 12, 16};
        shape = {'parallelogram', 'trapezoid', 'curved', 'general'};
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
        %% One-dimensional bases
        function basis1DDerivativesMatchFiniteDifferences(tc, basis, p)
            s = linspace(-1.2, 1.2, 17);
            h = 1e-6;
            [~, dPhi] = PolynomialBasis1D(s, p, basis);
            Pp = PolynomialBasis1D(s + h, p, basis);
            Pm = PolynomialBasis1D(s - h, p, basis);
            tc.verifyEqual(dPhi, (Pp - Pm)/(2*h), 'AbsTol', 1e-6);
        end

        function basis1DSpansPolynomialsOfDegreeP(tc, basis, p)
            % Interpolation of s^k, k <= p, in the basis is exact
            sn = nodal_distribution(p);
            V = PolynomialBasis1D(sn, p, basis);
            sq = linspace(-1, 1, 11);
            Q = PolynomialBasis1D(sq, p, basis);
            for k = 0:p
                c = V'\(sn.^k);
                tc.verifyEqual(Q'*c, sq'.^k, 'AbsTol', 1e-12);
            end
        end

        function legendreAndChebyshevValues(tc)
            s = linspace(-1, 1, 9);
            L = PolynomialBasis1D(s, 3, 'legendre');
            tc.verifyEqual(L(3,:), (3*s.^2 - 1)/2, 'AbsTol', 1e-15);
            tc.verifyEqual(L(4,:), (5*s.^3 - 3*s)/2, 'AbsTol', 1e-15);
            T = PolynomialBasis1D(s, 8, 'chebyshev');
            tc.verifyEqual(T, cos((0:8)'*acos(s)), 'AbsTol', 1e-13);
        end

        function lagrange1DMatchesParametricShapeFunctions(tc, p)
            xq = [linspace(-1, 1, 13)'; nodal_distribution(p)];
            [L, dL] = PolynomialBasis1D(xq, p, 'lagrange');
            [N, dN] = ShapeFunctions1D_SEM(p, xq);
            order = [1, p+1, 2:p];
            tc.verifyEqual(L(order,:), N, 'AbsTol', 1e-12);
            tc.verifyEqual(dL(order,:), dN, 'AbsTol', 1e-10);
        end

        function unknownBasisIsRejected(tc)
            tc.verifyError(@() PolynomialBasis1D(0, 2, 'bernstein'), ...
                'PolynomialBasis1D:basis');
        end

        %% Metric shape functions
        function metricFunctionsDoNotDependOnBasis(tc, basis, p, shape)
            [x, y, P] = usemTestElement(shape, p);
            [xq, yq] = interiorPoints(P, x, y);
            P.metricBasis = 'legendre';
            [N0, dNx0, dNy0] = ShapeFunctions2D_Metric(xq, yq, x, y, P);
            P.metricBasis = basis;
            [N, dNx, dNy] = ShapeFunctions2D_Metric(xq, yq, x, y, P);
            tc.verifyEqual(N, N0, 'AbsTol', 1e-11);
            tc.verifyEqual(dNx, dNx0, 'AbsTol', 1e-10);
            tc.verifyEqual(dNy, dNy0, 'AbsTol', 1e-10);
        end

        function metricFunctionsAreFrameInvariantForAllBases(tc, basis, shape)
            [x, y, P] = usemTestElement(shape, 3);
            P.metricBasis = basis;
            [xq, yq] = interiorPoints(P, x, y);
            [N0, dNx0, dNy0] = ShapeFunctions2D_Metric(xq, yq, x, y, P);
            th = -1.1; s = 0.2; t = [-3, 7];
            R = [cos(th) -sin(th); sin(th) cos(th)];
            XY  = s*([x y]*R') + t;
            XYq = s*([xq yq]*R') + t;
            [N1, dNx1, dNy1] = ShapeFunctions2D_Metric(XYq(:,1), XYq(:,2), ...
                                                       XY(:,1), XY(:,2), P);
            tc.verifyEqual(N1, N0, 'AbsTol', 1e-10);
            G0 = R*[dNx0(:)'; dNy0(:)']/s;
            tc.verifyEqual(dNx1(:)', G0(1,:), 'AbsTol', 1e-8);
            tc.verifyEqual(dNy1(:)', G0(2,:), 'AbsTol', 1e-8);
        end

        function defaultBasisIsLagrange(tc)
            [x, y, P] = usemTestElement('general', 4);
            [~, ~, ~, V0] = ShapeFunctions2D_Metric(x, y, x, y, P);
            P.metricBasis = 'lagrange';
            [~, ~, ~, V1] = ShapeFunctions2D_Metric(x, y, x, y, P);
            tc.verifyEqual(V0, V1);
        end

        %% Conditioning at high polynomial degree
        function lagrangeVandermondeIsPermutationOnAffineElements(tc, pHigh)
            % The basis is ordered by reference nodes, the columns by the
            % element node numbering: V is a permutation matrix, cond(V) = 1
            [x, y, P] = usemTestElement('parallelogram', pHigh);
            [~, ~, ~, V] = ShapeFunctions2D_Metric(x, y, x, y, P);
            n = numel(x);
            tc.verifyEqual(V*V', eye(n), 'AbsTol', 1e-11);
            tc.verifyEqual(sum(abs(V) > 0.5, 1), ones(1, n));
        end

        function monomialVandermondeGrowsExponentially(tc)
            % The reason for the change of the default basis: already on
            % the affine element cond(V) > 1e8 for p = 12
            [x, y, P] = usemTestElement('parallelogram', 12);
            P.metricBasis = 'monomial';
            [~, ~, ~, V] = ShapeFunctions2D_Metric(x, y, x, y, P);
            tc.verifyGreaterThan(cond(V), 1e8);
        end

        function highOrderMetricEqualsParametricOnParallelogram(tc, pHigh)
            [x, y, P] = usemTestElement('parallelogram', pHigh);
            [XI, ETA] = meshgrid(linspace(-0.95, 0.95, 7));
            Xg = Local2GlobalCoordinatesQuad(P, [XI(:), ETA(:)], [x, y]);
            Nm = ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, P);
            Np = ShapeFunctions2D_SEM(XI(:), ETA(:), P);
            tc.verifyEqual(Nm, Np, 'AbsTol', 1e-12);
        end

        function highOrderStiffnessIsSymmetricOnParallelogram(tc)
            [x, y, P] = usemTestElement('parallelogram', 10);
            K = Ke_2D(P, HookeMatrix2d(1, 0.3), x, y);
            tc.verifyLessThan(norm(K - K', 'fro')/norm(K, 'fro'), 1e-13);
        end

        function highOrderKroneckerPropertyOnCurvedElement(tc)
            % On distorted elements part of the ill-conditioning is
            % intrinsic (Lebesgue constant of the metric interpolation),
            % but the default basis removes the monomial part of it
            [x, y, P] = usemTestElement('curved', 12);
            P.metricBasis = 'monomial';
            e0 = max(max(abs(ShapeFunctions2D_Metric(x, y, x, y, P) - eye(numel(x)))));
            P.metricBasis = 'lagrange';
            e1 = max(max(abs(ShapeFunctions2D_Metric(x, y, x, y, P) - eye(numel(x)))));
            tc.verifyLessThan(e1, 1e-10);
            tc.verifyLessThan(e1, e0/100);
        end
    end
end

function [xq, yq] = interiorPoints(P, x, y)
[XI, ETA] = meshgrid(linspace(-0.9, 0.9, 4));
Xg = Local2GlobalCoordinatesQuad(P, [XI(:), ETA(:)], [x, y]);
xq = Xg(:,1); yq = Xg(:,2);
end
