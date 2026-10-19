classdef ElementMatricesTest < matlab.unittest.TestCase
    %ELEMENTMATRICESTEST Unsymmetric element stiffness and mass matrices
    %(Ke_2D, Me_2D, KMe_2D) and the lumping schemes.

    properties (TestParameter)
        p = {1, 2, 3, 4};
        shape = {'square', 'parallelogram', 'trapezoid', 'curved', 'general'};
        distorted = {'trapezoid', 'curved', 'general'};
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
            tc.C = HookeMatrix2d(1, 0.3, 1);
        end
    end

    methods (Test)
        %% Stiffness matrix
        function rigidBodyModesAreInRightAndLeftNullSpace(tc, shape, p)
            % Trial side: the metric functions reproduce linear fields, so
            % B_R * u_rigid = 0. Test side: the parametric functions are a
            % partition of unity and isoparametric, so B_L * u_rigid = 0.
            [x, y, P] = usemTestElement(shape, p);
            K = Ke_2D(P, tc.C, x, y);
            R = rigidBodyModes(x, y);
            tc.verifyLessThan(norm(K*R)/norm(K), 1e-12);
            tc.verifyLessThan(norm(R'*K)/norm(K), 1e-12);
        end

        function stiffnessHasExactlyThreeZeroEigenvalues(tc, shape, p)
            [x, y, P] = usemTestElement(shape, p);
            K = Ke_2D(P, tc.C, x, y);
            s = svd(K);
            tc.verifyEqual(nnz(s < 1e-10*s(1)), 3, ...
                'Expected rank(K) = ndof - 3 (no spurious zero-energy modes).');
        end

        function stiffnessIsSymmetricOnAffineElements(tc, p)
            % For parallelograms trial and test spaces coincide, so the
            % unsymmetric element reduces to the classical spectral element.
            [x, y, P] = usemTestElement('parallelogram', p);
            K = Ke_2D(P, tc.C, x, y);
            tc.verifyLessThan(norm(K - K', 'fro')/norm(K, 'fro'), 1e-12);
            Ksym = usemKeSymmetric(P, tc.C, x, y);
            tc.verifyEqual(K, Ksym, 'AbsTol', 1e-12*norm(K, 'fro'));
        end

        function stiffnessIsUnsymmetricOnDistortedElements(tc, distorted)
            [x, y, P] = usemTestElement(distorted, 2);
            K = Ke_2D(P, tc.C, x, y);
            tc.verifyGreaterThan(norm(K - K', 'fro')/norm(K, 'fro'), 1e-4);
        end

        function constantStrainEnergyIsExact(tc, distorted, p)
            % For u = eps0 * x (constant strain) and any test function,
            % K*u equals the consistent nodal force of the constant stress.
            % In particular u'*K*u = eps0'*C*eps0 * area.
            [x, y, P] = usemTestElement(distorted, p);
            K = Ke_2D(P, tc.C, x, y);
            ex = 1e-3; ey = -4e-4; gxy = 2e-4;
            u = zeros(2*numel(x), 1);
            u(1:2:end) = ex*x + gxy/2*y;
            u(2:2:end) = ey*y + gxy/2*x;
            eps0 = [ex; ey; gxy];
            tc.verifyEqual(u'*K*u, eps0'*tc.C*eps0*elementArea(P, x, y), ...
                'RelTol', 1e-10);
        end

        function stiffnessScalesWithMaterial(tc)
            [x, y, P] = usemTestElement('general', 3);
            K1 = Ke_2D(P, tc.C, x, y);
            K2 = Ke_2D(P, 7.5*tc.C, x, y);
            tc.verifyEqual(K2, 7.5*K1, 'AbsTol', 1e-12*norm(K2, 'fro'));
        end

        function stiffnessIsInvariantUnderRigidRotation(tc, distorted)
            % Rotating the element by theta must give T*K*T' with the
            % nodal rotation matrix T (frame invariance of the formulation).
            [x, y, P] = usemTestElement(distorted, 2);
            K = Ke_2D(P, tc.C, x, y);
            th = 0.9; R = [cos(th) -sin(th); sin(th) cos(th)];
            XY = [x y]*R';
            Kr = Ke_2D(P, tc.C, XY(:,1), XY(:,2));
            T = kron(eye(numel(x)), R);
            tc.verifyEqual(Kr, T*K*T', 'AbsTol', 1e-10*norm(K, 'fro'));
        end

        function combinedRoutineMatchesSeparateStiffness(tc, shape)
            [x, y, P] = usemTestElement(shape, 3);
            K = Ke_2D(P, tc.C, x, y);
            [K2, M2] = KMe_2D(P, tc.C, 2.5, x, y, 1);
            tc.verifyEqual(K2, K, 'AbsTol', 1e-12*norm(K, 'fro'));
            % total mass is integrated exactly by both routines
            ex = repmat([1; 0], numel(x), 1);
            tc.verifyEqual(ex'*M2*ex, 2.5*elementArea(P, x, y), 'RelTol', 1e-10);
        end

        function combinedRoutineMatchesSeparateMassOnAffineElements(tc, p)
            [x, y, P] = usemTestElement('parallelogram', p);
            [~, M2] = KMe_2D(P, tc.C, 2.5, x, y, 1);
            M = Me_2D(P, 2.5, x, y, 1);
            tc.verifyEqual(M2, M, 'AbsTol', 1e-12*norm(M, 'fro'));
        end

        %% Mass matrix
        function consistentMassReproducesTotalMass(tc, shape, p)
            rho = 7850;
            [x, y, P] = usemTestElement(shape, p);
            M = Me_2D(P, rho, x, y, 1);
            ex = repmat([1; 0], numel(x), 1);
            ey = repmat([0; 1], numel(x), 1);
            m = rho*elementArea(P, x, y);
            tc.verifyEqual(ex'*M*ex, m, 'RelTol', 1e-10);
            tc.verifyEqual(ey'*M*ey, m, 'RelTol', 1e-10);
            tc.verifyEqual(ex'*M*ey, 0, 'AbsTol', 1e-12*m);
        end

        function massIsSymmetricOnAffineElements(tc, p)
            [x, y, P] = usemTestElement('parallelogram', p);
            M = Me_2D(P, 1, x, y, 1);
            tc.verifyLessThan(norm(M - M', 'fro')/norm(M, 'fro'), 1e-12);
        end

        function nodalQuadratureGivesDiagonalMass(tc, distorted, p)
            % GLL nodes + GLL quadrature: the parametric AND the metric
            % functions are interpolatory at the quadrature points, so the
            % unsymmetric mass matrix is diagonal as well.
            [x, y, P] = usemTestElement(distorted, p, 2);
            M = Me_2D(P, 1, x, y, 2);
            tc.verifyLessThan(norm(M - diag(diag(M)), 'fro')/norm(M, 'fro'), 1e-10);
            tc.verifyGreaterThan(min(diag(M)), 0);
            tc.verifyEqual(sum(diag(M))/2, elementArea(P, x, y), 'RelTol', 1e-10);
        end

        function rowSumAndHRZLumpingPreserveMass(tc, distorted)
            [x, y, P] = usemTestElement(distorted, 3);
            Mc = Me_2D(P, 1, x, y, 1);
            for scheme = [3 4]
                M = Me_2D(P, 1, x, y, scheme);
                tc.verifyTrue(isdiag(M));
                tc.verifyEqual(sum(M(:)), sum(Mc(:)), 'RelTol', 1e-12);
            end
            tc.verifyEqual(RowSum_Method(Mc), diag(sum(Mc, 2)));
            Mh = HRZ_Method(Mc);
            tc.verifyEqual(diag(Mh)/sum(diag(Mh)), diag(Mc)/sum(diag(Mc)), 'AbsTol', 1e-14);
        end
    end

    methods (Test, TestTags = {'KnownIssue'})
        function nodalQuadratureFallsBackForNonGLLNodes(tc, p)
            % lumping = 2 with ndist = 1 or 3 is documented to fall back to
            % HRZ / row-sum lumping, but the quadrature points are never
            % defined in that branch (see docs/known-issues.md).
            [x, y, P] = usemTestElement('trapezoid', p, 1);
            M = Me_2D(P, 1, x, y, 2);
            tc.verifyTrue(isdiag(M));
            [~, M2] = KMe_2D(P, tc.C, 1, x, y, 2);
            tc.verifyTrue(isdiag(M2));
        end

        function combinedRoutineMatchesSeparateMassOnDistortedElements(tc, distorted)
            % KMe_2D (used by the dynamic, modal and frequency-domain
            % analyses) integrates the unsymmetric mass matrix with p+1
            % Gauss points per direction, Me_2D with 2p. On distorted
            % elements the metric functions are not polynomials of degree
            % p in (xi, eta), so p+1 points under-integrate N'*M*detJ
            % (see docs/known-issues.md).
            [x, y, P] = usemTestElement(distorted, 3);
            [~, M2] = KMe_2D(P, tc.C, 2.5, x, y, 1);
            M = Me_2D(P, 2.5, x, y, 1);
            tc.verifyEqual(M2, M, 'AbsTol', 1e-10*norm(M, 'fro'));
        end
    end
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
