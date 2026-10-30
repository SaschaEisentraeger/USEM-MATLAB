classdef GeometryMappingTest < matlab.unittest.TestCase
    %GEOMETRYMAPPINGTEST Isoparametric mapping, Jacobian and its inverse.

    properties (TestParameter)
        shape = {'square', 'parallelogram', 'trapezoid', 'curved', 'general'};
        p = {1, 2, 4};
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
        function jacobianOfAffineMapIsConstant(tc, p)
            % x = A*[xi; eta] + c  =>  J = [dx/dxi dy/dxi; dx/deta dy/deta] = A'
            [x, y, P] = usemTestElement('parallelogram', p);
            A = [1.0 0.3; 0.15 0.6];               % from the corner nodes
            for pt = [-0.7 0.2; 0.1 -0.5; 0.9 0.9]'
                J = Jacobi_Matrix(P, pt', x, y);
                tc.verifyEqual(J, A', 'AbsTol', 1e-12);
            end
        end

        function jacobianMatchesFiniteDifferences(tc, shape, p)
            [x, y, P] = usemTestElement(shape, p);
            pt = [0.31, -0.42]; h = 1e-6;
            J = Jacobi_Matrix(P, pt, x, y);
            X = @(q) Local2GlobalCoordinatesQuad(P, q, [x y]);
            Jfd = [(X(pt + [h 0]) - X(pt - [h 0]))/(2*h); ...
                   (X(pt + [0 h]) - X(pt - [0 h]))/(2*h)];
            tc.verifyEqual(J, Jfd, 'AbsTol', 1e-8);
            tc.verifyGreaterThan(det(J), 0, 'Element must not be inverted.');
        end

        function localToGlobalInterpolatesNodes(tc, shape, p)
            [x, y, P] = usemTestElement(shape, p);
            xi = nodal_distribution(P.pxi, P.ndist);
            Xc = Local2GlobalCoordinatesQuad(P, [-1 -1; 1 -1; 1 1; -1 1], [x y]);
            tc.verifyEqual(Xc, [x(1:4) y(1:4)], 'AbsTol', 1e-13);
            if p > 1
                % first node on edge 1 (eta = -1)
                X5 = Local2GlobalCoordinatesQuad(P, [xi(2) -1], [x y]);
                tc.verifyEqual(X5, [x(5) y(5)], 'AbsTol', 1e-13);
            end
        end

        function areaFromJacobianIsExact(tc)
            % Trapezoid with straight edges: area by the shoelace formula
            [x, y, P] = usemTestElement('trapezoid', 2);
            [g, w] = lgwt(4, -1, 1);
            A = 0;
            for i = 1:4
                for j = 1:4
                    A = A + det(Jacobi_Matrix(P, [g(i) g(j)], x, y))*w(i)*w(j);
                end
            end
            xc = x(1:4); yc = y(1:4);
            Ashoelace = 0.5*abs(sum(xc.*circshift(yc,-1) - circshift(xc,-1).*yc));
            tc.verifyEqual(A, Ashoelace, 'RelTol', 1e-13);
        end
    end

    methods (Test)
    end
end
