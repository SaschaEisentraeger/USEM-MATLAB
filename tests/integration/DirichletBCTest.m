classdef DirichletBCTest < matlab.unittest.TestCase
    %DIRICHLETBCTEST Homogeneous and prescribed displacements
    %(DirichletBCs_homogeneous, DirichletBCs_nonhomogeneous).
    %
    %   Uniaxial stretch of the strip [0,3] x [0,2] on a distorted mesh:
    %   u_x = 0 on the left edge, u_x = d on the right edge, u_y = 0 at the
    %   origin. Plane stress, exact solution u_x = d*x/L,
    %   u_y = -nu*d*y/L.

    properties (TestParameter)
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
        function prescribedStretchIsExact(tc, p)
            nu = 0.3; L = 3; d = 1e-3;
            [nodes, topo] = usemTestMesh(3, 2, L, 2, 2);
            nodes = usemDistortMesh(nodes, topo, 'combined', 0.3);
            model = usemBuildModel(nodes, topo, p);
            K = usemAssemble(model, HookeMatrix2d(1, nu, 1));
            x = model.xCoor; y = model.yCoor;
            left   = find(abs(x) < 1e-12);
            right  = find(abs(x - L) < 1e-12);
            origin = find(abs(x) < 1e-12 & abs(y) < 1e-12);

            bc = struct('set_name', {'left', 'origin', 'right'}, ...
                        'amplitude', {[], [], d}, ...
                        'type', {'Ux 0', 'Uy 0', 'Ux'}, ...
                        'nset', {left, origin, right}, ...
                        'elset', {[], [], []}, ...
                        'step_name', {'Step-1', 'Step-1', 'Step-1'});
            bc = bc(:);         % the readers return column struct arrays
            [~, fixed] = evalc('DirichletBCs_homogeneous(bc)');
            tc.verifyEqual(sort(fixed(:)), sort([2*left-1; 2*origin]));

            F = zeros(model.num_DoF, 1);
            [~, K, F] = evalc('DirichletBCs_nonhomogeneous(K, F, model.num_DoF, bc, x, y)');
            act = setdiff((1:model.num_DoF)', unique(fixed));
            u = zeros(model.num_DoF, 1);
            u(act) = K(act, act)\F(act);

            tc.verifyEqual(u(2*right-1), d*ones(numel(right), 1), 'AbsTol', 1e-14);
            tc.verifyEqual(u(1:2:end), d*x/L, 'AbsTol', 1e-12);
            tc.verifyEqual(u(2:2:end), -nu*d*y/L, 'AbsTol', 1e-12);
        end
    end

    methods (Test)
        function threeNodeSetIsNotMistakenForGenerateTriplet(tc)
            % Node sets are used as they are; "generate" sets are expanded
            % by the reader. A set {2, 8, 3} (e.g. an edge of a p = 2
            % element) used to be turned into the range 2:3:8.
            bc = struct('set_name', 'edge', 'amplitude', [], 'type', 'Ux 0', ...
                        'nset', [2; 8; 3], 'elset', [], 'step_name', 'Step-1');
            [~, fixed] = evalc('DirichletBCs_homogeneous(bc)');
            tc.verifyEqual(sort(fixed(:)), [3; 5; 15]);
        end
    end
end
