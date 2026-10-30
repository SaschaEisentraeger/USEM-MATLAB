classdef MetricHybridTest < matlab.unittest.TestCase
    %METRICHYBRIDTEST Hybrid trial space (ShapeFunctions2D_Metric, method 2)
    %
    %   T_p = P_p(x,y) + span{L_a(xi) L_b(eta) : a, b <= p, a + b > p}:
    %   complete in Cartesian coordinates, frame invariant, identical to
    %   the parametric space on straight-edged (bilinearly mapped) elements,
    %   Lebesgue constant close to that of the parametric interpolation
    %   (docs/chapters/metric-basis.tex). Also covers the helpers
    %   ReferenceNodes2D.

    properties (TestParameter)
        p = {1, 2, 3, 4, 6};
        shape = {'parallelogram', 'trapezoid', 'curved', 'general'};
        straight = {'square', 'parallelogram', 'trapezoid'};
        curved = {'curved', 'general'};
    end

    properties
        C       % plane stress elasticity matrix
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'code', 'subroutines'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'tests', 'helpers')));
        end

        function setMaterial(tc)
            tc.C = HookeMatrix2d(1, 0.3);
        end
    end

    methods (Test)
        %% Helpers
        function referenceNodesAreTheSEMNodes(tc, p)
            P.pxi = p; P.peta = p;
            [xi, eta] = ReferenceNodes2D(P);
            N = ShapeFunctions2D_SEM(xi, eta, P);
            tc.verifyEqual(N, eye((p+1)^2), 'AbsTol', 1e-12);
        end

        function referenceNodesForAnisotropicDegrees(tc)
            P.pxi = 4; P.peta = 2;
            [xi, eta] = ReferenceNodes2D(P);
            tc.verifyEqual(ShapeFunctions2D_SEM(xi, eta, P), eye(15), 'AbsTol', 1e-12);
        end

        %% Shape functions
        function hybridFunctionsHaveKroneckerProperty(tc, p, shape)
            [x, y, P] = usemTestElement(shape, p);
            [xiN, etaN] = ReferenceNodes2D(P);
            N = ShapeFunctions2D_Metric(x, y, x, y, P, 2, [xiN(:), etaN(:)]);
            tc.verifyEqual(N, eye(numel(x)), 'AbsTol', 1e-12);
        end

        function hybridFunctionsAreCompleteInCartesianCoordinates(tc, p, shape)
            [x, y, P] = usemTestElement(shape, p);
            [xq, yq, IP] = interiorPoints(P, x, y);
            [N, dNx, dNy] = ShapeFunctions2D_Metric(xq, yq, x, y, P, 2, IP);
            xc = mean(x); yc = mean(y);
            for a = 0:p
                for b = 0:p-a
                    f = (x-xc).^a.*(y-yc).^b;
                    tc.verifyEqual(N'*f, (xq-xc).^a.*(yq-yc).^b, 'AbsTol', 1e-12);
                    tc.verifyEqual(dNx'*f, a*(xq-xc).^max(a-1,0).*(yq-yc).^b, 'AbsTol', 1e-11);
                    tc.verifyEqual(dNy'*f, b*(xq-xc).^a.*(yq-yc).^max(b-1,0), 'AbsTol', 1e-11);
                end
            end
        end

        function derivativesMatchFiniteDifferences(tc, shape)
            [x, y, P] = usemTestElement(shape, 3);
            [xq, yq, IP] = interiorPoints(P, x, y);
            [~, dNx, dNy] = ShapeFunctions2D_Metric(xq, yq, x, y, P, 2, IP);
            % Central differences along xi and eta, chain rule in (x,y)
            h = 1e-6;
            for d = [1 0; 0 1]'
                Xp = Local2GlobalCoordinatesQuad(P, IP + h*d', [x, y]);
                Xm = Local2GlobalCoordinatesQuad(P, IP - h*d', [x, y]);
                Np = ShapeFunctions2D_Metric(Xp(:,1), Xp(:,2), x, y, P, 2, IP + h*d');
                Nm = ShapeFunctions2D_Metric(Xm(:,1), Xm(:,2), x, y, P, 2, IP - h*d');
                dX = (Xp - Xm)/(2*h);
                tc.verifyEqual(dNx.*dX(:,1)' + dNy.*dX(:,2)', (Np-Nm)/(2*h), 'AbsTol', 1e-6);
            end
        end

        function methodIsSelectedByTheDegreeStructure(tc)
            [x, y, P] = usemTestElement('general', 3);
            [xq, yq, IP] = interiorPoints(P, x, y);
            N0 = ShapeFunctions2D_Metric(xq, yq, x, y, P, 2, IP);
            P.metricMethod = 2;
            N1 = ShapeFunctions2D_Metric(xq, yq, x, y, P, [], IP);
            tc.verifyEqual(N1, N0);
        end

        function negativeJacobianIsAllowed(tc, shape)
            % Mirrored element (det J < 0 everywhere): the functions are
            % evaluated through (xi,eta), completeness is unaffected
            [x, y, P] = usemTestElement(shape, 3);
            y = -y;
            [xq, yq, IP] = interiorPoints(P, x, y);
            [N, dNx, dNy] = ShapeFunctions2D_Metric(xq, yq, x, y, P, 2, IP);
            xc = mean(x); yc = mean(y);
            for a = 0:3
                for b = 0:3-a
                    f = (x-xc).^a.*(y-yc).^b;
                    tc.verifyEqual(N'*f, (xq-xc).^a.*(yq-yc).^b, 'AbsTol', 1e-12);
                    tc.verifyEqual(dNx'*f, a*(xq-xc).^max(a-1,0).*(yq-yc).^b, 'AbsTol', 1e-11);
                    tc.verifyEqual(dNy'*f, b*(xq-xc).^a.*(yq-yc).^max(b-1,0), 'AbsTol', 1e-11);
                end
            end
        end

        function singularJacobianAtQueryPointIsRejected(tc)
            % Edge 2-3 collapsed (node 3 onto node 2): x_eta = 0 on xi = 1,
            % so det J = 0 there and the gradients are undefined; the
            % Vandermonde matrix stays regular (xi*eta separates the nodes)
            [x, y, P] = usemTestElement('square', 1);
            x(3) = x(2); y(3) = y(2);
            tc.verifyError(@() ShapeFunctions2D_Metric(x(2), y(2), x, y, P, 2, ...
                [1, 0]), 'ShapeFunctions2D_Metric:singularJacobian');
        end

        function hybridEqualsParametricOnStraightEdgedElements(tc, straight, p)
            % Bilinear map: x, y in Q_1(xi,eta), hence P_p(x,y) is in
            % Q_p(xi,eta) and T_p = Q_p(xi,eta) (dimension count)
            [x, y, P] = usemTestElement(straight, p);
            [xq, yq, IP] = interiorPoints(P, x, y);
            Nm = ShapeFunctions2D_Metric(xq, yq, x, y, P, 2, IP);
            Np = ShapeFunctions2D_SEM(IP(:,1), IP(:,2), P);
            tc.verifyEqual(Nm, Np, 'AbsTol', 1e-12);
        end

        function hybridDiffersFromParametricOnCurvedElements(tc, curved)
            [x, y, P] = usemTestElement(curved, 3);
            [xq, yq, IP] = interiorPoints(P, x, y);
            Nm = ShapeFunctions2D_Metric(xq, yq, x, y, P, 2, IP);
            Np = ShapeFunctions2D_SEM(IP(:,1), IP(:,2), P);
            tc.verifyGreaterThan(max(abs(Nm(:) - Np(:))), 1e-3);
        end

        function hybridFunctionsAreFrameInvariant(tc, shape)
            [x, y, P] = usemTestElement(shape, 3);
            [xq, yq, IP] = interiorPoints(P, x, y);
            [N0, dNx0, dNy0] = ShapeFunctions2D_Metric(xq, yq, x, y, P, 2, IP);
            th = 0.7; s = 3.5; t = [10, -4];
            R = [cos(th) -sin(th); sin(th) cos(th)];
            XY  = s*([x y]*R') + t;
            XYq = s*([xq yq]*R') + t;
            [N1, dNx1, dNy1] = ShapeFunctions2D_Metric(XYq(:,1), XYq(:,2), ...
                                   XY(:,1), XY(:,2), P, 2, IP);
            tc.verifyEqual(N1, N0, 'AbsTol', 1e-11);
            G0 = R*[dNx0(:)'; dNy0(:)']/s;
            tc.verifyEqual(dNx1(:)', G0(1,:), 'AbsTol', 1e-10);
            tc.verifyEqual(dNy1(:)', G0(2,:), 'AbsTol', 1e-10);
        end

        function lebesgueConstantStaysSmallAtHighDegree(tc, curved)
            % Q_p(s,t) (method 1): > 1e3 at p = 12; hybrid: close to the
            % parametric value (5.1 at p = 12)
            [x, y, P] = usemTestElement(curved, 12);
            [xi, eta] = refGrid(41, 1);
            Xg = Local2GlobalCoordinatesQuad(P, [xi, eta], [x, y]);
            M1 = ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, P, 1);
            M2 = ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, P, 2, [xi, eta]);
            tc.verifyGreaterThan(max(sum(abs(M1), 1)), 1e3);
            tc.verifyLessThan(max(sum(abs(M2), 1)), 7);
        end

        %% Element matrices
        function rigidBodyModesAreInRightAndLeftNullSpace(tc, shape, p)
            [x, y, P] = usemTestElement(shape, p);
            P.metricMethod = 2;
            K = Ke_2D(P, tc.C, x, y);
            R = rigidBodyModes(x, y);
            tc.verifyLessThan(norm(K*R)/norm(K), 1e-12);
            tc.verifyLessThan(norm(R'*K)/norm(K), 1e-12);
        end

        function stiffnessHasExactlyThreeZeroEigenvalues(tc, shape, p)
            [x, y, P] = usemTestElement(shape, p);
            P.metricMethod = 2;
            s = svd(Ke_2D(P, tc.C, x, y));
            tc.verifyEqual(nnz(s < 1e-10*s(1)), 3);
        end

        function stiffnessIsClassicalOnStraightEdgedElements(tc, straight)
            [x, y, P] = usemTestElement(straight, 4);
            P.metricMethod = 2;
            K = Ke_2D(P, tc.C, x, y);
            Ksym = usemKeSymmetric(P, tc.C, x, y);
            tc.verifyEqual(K, Ksym, 'AbsTol', 1e-12*norm(K, 'fro'));
        end

        function stiffnessIsUnsymmetricOnCurvedElements(tc, curved)
            [x, y, P] = usemTestElement(curved, 3);
            P.metricMethod = 2;
            K = Ke_2D(P, tc.C, x, y);
            tc.verifyGreaterThan(norm(K - K', 'fro')/norm(K, 'fro'), 1e-4);
        end

    end
end

function [xi, eta] = refGrid(n, a)
[XI, ETA] = meshgrid(linspace(-a, a, n));
xi = XI(:); eta = ETA(:);
end

function [xq, yq, IP] = interiorPoints(P, x, y)
[xi, eta] = refGrid(4, 0.9);
IP = [xi, eta];
Xg = Local2GlobalCoordinatesQuad(P, IP, [x, y]);
xq = Xg(:,1); yq = Xg(:,2);
end

function R = rigidBodyModes(x, y)
n = numel(x);
R = zeros(2*n, 3);
R(1:2:end, 1) = 1;
R(2:2:end, 2) = 1;
R(1:2:end, 3) = -y;
R(2:2:end, 3) = x;
end

function A = elementArea(P, x, y)
[g, w] = lgwt(2*P.pxi + 2, -1, 1);
A = 0;
for i = 1:numel(g)
    for j = 1:numel(g)
        A = A + det(Jacobi_Matrix(P, [g(i) g(j)], x, y))*w(i)*w(j);
    end
end
end
