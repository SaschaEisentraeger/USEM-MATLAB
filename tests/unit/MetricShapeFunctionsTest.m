classdef MetricShapeFunctionsTest < matlab.unittest.TestCase
    %METRICSHAPEFUNCTIONSTEST Metric (trial) shape functions written in
    %corotational skew coordinates (ShapeFunctions2D_Metric,
    %CorotationalSkewCoordinates).

    properties (TestParameter)
        p = {1, 2, 3, 4};
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
        function metricFunctionsHaveKroneckerProperty(tc, p, shape)
            [x, y, P] = usemTestElement(shape, p);
            N = ShapeFunctions2D_Metric(x, y, x, y, P);
            tc.verifyEqual(N, eye(numel(x)), 'AbsTol', 1e-9);
        end

        function metricFunctionsAreCompleteInCartesianCoordinates(tc, p, shape)
            % Any polynomial of total degree <= p in (x,y) is reproduced
            % exactly, together with its gradient, also on distorted
            % elements. This is the property that makes the element
            % immune to mesh distortion.
            [x, y, P] = usemTestElement(shape, p);
            [xq, yq] = interiorPoints(P, x, y);
            [N, dNx, dNy] = ShapeFunctions2D_Metric(xq, yq, x, y, P);
            for a = 0:p
                for b = 0:p-a
                    f = x.^a .* y.^b;
                    tc.verifyEqual(N'*f, xq.^a .* yq.^b, 'AbsTol', 1e-9);
                    tc.verifyEqual(dNx'*f, a*xq.^max(a-1,0) .* yq.^b, 'AbsTol', 1e-8);
                    tc.verifyEqual(dNy'*f, xq.^a .* b.*yq.^max(b-1,0), 'AbsTol', 1e-8);
                end
            end
        end

        function derivativesMatchFiniteDifferences(tc, shape)
            [x, y, P] = usemTestElement(shape, 3);
            [xq, yq] = interiorPoints(P, x, y);
            [~, dNx, dNy] = ShapeFunctions2D_Metric(xq, yq, x, y, P);
            h = 1e-6;
            Np = ShapeFunctions2D_Metric(xq+h, yq, x, y, P);
            Nm = ShapeFunctions2D_Metric(xq-h, yq, x, y, P);
            tc.verifyEqual(dNx, (Np-Nm)/(2*h), 'AbsTol', 1e-6);
            Np = ShapeFunctions2D_Metric(xq, yq+h, x, y, P);
            Nm = ShapeFunctions2D_Metric(xq, yq-h, x, y, P);
            tc.verifyEqual(dNy, (Np-Nm)/(2*h), 'AbsTol', 1e-6);
        end

        function metricEqualsParametricOnParallelogram(tc, p)
            % For an affine element the skew coordinates coincide with
            % (xi, eta), so trial and test spaces are identical.
            [x, y, P] = usemTestElement('parallelogram', p);
            [xi, eta] = refGrid(5);
            Xg = Local2GlobalCoordinatesQuad(P, [xi, eta], [x, y]);
            Nm = ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, P);
            Np = ShapeFunctions2D_SEM(xi, eta, P);
            tc.verifyEqual(Nm, Np, 'AbsTol', 1e-10);
        end

        function metricFunctionsAreFrameInvariant(tc, p, shape)
            % Rigid rotation, translation and uniform scaling of the
            % element leave the shape function values unchanged; the
            % gradients rotate with the frame and scale with 1/s.
            [x, y, P] = usemTestElement(shape, p);
            [xq, yq] = interiorPoints(P, x, y);
            [N0, dNx0, dNy0] = ShapeFunctions2D_Metric(xq, yq, x, y, P);
            th = 0.7; s = 3.5; t = [10, -4];
            R = [cos(th) -sin(th); sin(th) cos(th)];
            XY  = s*([x y]*R') + t;
            XYq = s*([xq yq]*R') + t;
            [N1, dNx1, dNy1] = ShapeFunctions2D_Metric(XYq(:,1), XYq(:,2), ...
                                                       XY(:,1), XY(:,2), P);
            tc.verifyEqual(N1, N0, 'AbsTol', 1e-9);
            G0 = R*[dNx0(:)'; dNy0(:)']/s;      % rotated reference gradient
            tc.verifyEqual(dNx1(:)', G0(1,:), 'AbsTol', 1e-8);
            tc.verifyEqual(dNy1(:)', G0(2,:), 'AbsTol', 1e-8);
        end

        function globalCoordinateVariantIsFrameDependent(tc)
            % method = 0 (centroid-based Cartesian monomials, Rajendran &
            % Liew 2003) is NOT rotation invariant for p >= 2. This is the
            % reason why the code uses skew coordinates by default.
            [x, y, P] = usemTestElement('general', 2);
            [xq, yq] = interiorPoints(P, x, y);
            N0 = ShapeFunctions2D_Metric(xq, yq, x, y, P, 0);
            th = 0.7; R = [cos(th) -sin(th); sin(th) cos(th)];
            XY = [x y]*R'; XYq = [xq yq]*R';
            N1 = ShapeFunctions2D_Metric(XYq(:,1), XYq(:,2), XY(:,1), XY(:,2), P, 0);
            tc.verifyGreaterThan(max(abs(N1(:) - N0(:))), 1e-4);
            % ... but it is still complete in (x, y)
            f = x.^2 + 3*x.*y - y.^2;
            tc.verifyEqual(N0'*f, xq.^2 + 3*xq.*yq - yq.^2, 'AbsTol', 1e-9);
        end

        function singularInterpolationIsRejected(tc)
            % interior node moved onto the mid-side node of edge 3: the
            % metric shape functions are undefined
            [x, y, P] = usemTestElement('square', 2);
            x(9) = x(7); y(9) = y(7);
            tc.verifyError(@() ShapeFunctions2D_Metric(1, 0.5, x, y, P), ...
                'ShapeFunctions2D_Metric:singular');
        end

        function skewCoordinatesOfAffineElementAreReferenceCoordinates(tc, p)
            [x, y, P] = usemTestElement('parallelogram', p);
            [xi, eta] = refGrid(4);
            Xg = Local2GlobalCoordinatesQuad(P, [xi, eta], [x, y]);
            [s1, s2, J] = CorotationalSkewCoordinates(Xg(:,1), Xg(:,2), x, y, P);
            tc.verifyEqual([s1, s2], [xi, eta], 'AbsTol', 1e-12);
            tc.verifyEqual(J, Jacobi_Matrix(P, [0 0], x, y), 'AbsTol', 1e-14);
        end

        function skewOriginIsElementCentre(tc, shape)
            [x, y, P] = usemTestElement(shape, 2);
            Xc = Local2GlobalCoordinatesQuad(P, [0 0], [x, y]);
            [s1, s2] = CorotationalSkewCoordinates(Xc(1), Xc(2), x, y, P);
            tc.verifyEqual([s1, s2], [0, 0], 'AbsTol', 1e-14);
        end
    end
end

function [xi, eta] = refGrid(n)
[XI, ETA] = meshgrid(linspace(-0.9, 0.9, n));
xi = XI(:); eta = ETA(:);
end

function [xq, yq] = interiorPoints(P, x, y)
[xi, eta] = refGrid(4);
Xg = Local2GlobalCoordinatesQuad(P, [xi, eta], [x, y]);
xq = Xg(:,1); yq = Xg(:,2);
end
