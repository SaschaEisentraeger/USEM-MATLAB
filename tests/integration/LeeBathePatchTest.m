classdef LeeBathePatchTest < matlab.unittest.TestCase
    %LEEBATHEPATCHTEST Load-driven patch tests of Lee & Bathe (1993).
    %
    %   Plate [0,10] x [-1,1] (E = 70 GPa, nu = 0.3, plane stress), u_x = 0
    %   on the left edge, u_y = 0 at (0,-1), edge loads on the right edge.
    %   Meshes: undistorted, angular and mid-edge (curved) distortion of a
    %   two-element mesh, and a five-element mesh (Rajendran & Liew 2003).
    %
    %   The input data (mesh, loads, boundary conditions after they were
    %   adapted to the spectral elements) is stored in
    %   tests/fixtures/reference/*.mat; see make_reference_fixtures.m. The
    %   tests repeat the static analysis of call_Static.m with the current
    %   code and check
    %     - the quadratic test (pure bending) against the analytical
    %       solution for every p >= 2 and every mesh,
    %     - the cubic test for p = 3: identical results on all meshes and
    %       agreement with the stored reference solution.

    properties (TestParameter)
        quadraticCase = {'PatchTestLeeQuadratic_undistorted_ndof30', ...
                         'PatchTestLeeQuadratic_angular_ndof30', ...
                         'PatchTestLeeQuadratic_angular_ndof56', ...
                         'PatchTestLeeQuadratic_angular_ndof90', ...
                         'PatchTestLeeQuadratic_angular_ndof182', ...
                         'PatchTestLeeQuadratic_midedge_ndof30'};
        cubicCase = {'PatchTestLeeCubic_undistorted_ndof56', ...
                     'PatchTestLeeCubic_angular_ndof56', ...
                     'PatchTestLeeCubic_midedge_ndof56', ...
                     'PatchTestLeeCubic_Rajendran_adj_ndof104'};
        anyCase = {'PatchTestLeeQuadratic_undistorted_ndof12', ...
                   'PatchTestLeeQuadratic_angular_ndof12', ...
                   'PatchTestLeeQuadratic_midedge_ndof12', ...
                   'PatchTestLeeCubic_undistorted_ndof30', ...
                   'PatchTestLeeCubic_midedge_ndof30'};
    end

    properties (Constant)
        E  = 7e10;
        nu = 0.3;
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
        function topologyIsRegeneratedIdentically(tc, quadraticCase)
            S = loadCase(quadraticCase);
            [TN, x, y] = pSpectralElement_NodesTopology(S.topology, ...
                S.nodes(:,1), S.nodes(:,2), S.p, S.p_org);
            tc.verifyEqual(TN, S.Topology_Node);
            tc.verifyEqual([x y], [S.xCoor S.yCoor], 'AbsTol', 1e-12);
        end

        function quadraticPatchTestMatchesPureBending(tc, quadraticCase)
            % sigma_x = a*y, u_x = a*x*y/E,
            % u_y = -a/(2E)*(x^2 + nu*y^2) + a*nu/(2E)   (u_y(0,-1) = 0).
            % The load intensity a follows from the consistent nodal
            % forces on the right edge: sum(F_x,i * y_i) = 2a/3.
            S = loadCase(quadraticCase);
            [u, F] = solveStatic(S);
            right = abs(S.xCoor - 10) < 1e-9;
            a = 1.5*sum(F(2*find(right)-1).*S.yCoor(right));
            tc.assertGreaterThan(abs(a), 0);
            x = S.xCoor; y = S.yCoor;
            ux = a*x.*y/tc.E;
            uy = -a/(2*tc.E)*(x.^2 + tc.nu*y.^2) + a*tc.nu/(2*tc.E);
            tc.verifyEqual(u(1:2:end), ux, 'AbsTol', 1e-9*max(abs(uy)));
            tc.verifyEqual(u(2:2:end), uy, 'AbsTol', 1e-9*max(abs(uy)));
        end

        function cubicPatchTestIsMeshIndependent(tc, cubicCase)
            % All meshes must give the tip displacement of the
            % undistorted mesh (cubic field, p = 3).
            ref = loadCase('PatchTestLeeCubic_undistorted_ndof56');
            uRef = solveStatic(ref);
            S = loadCase(cubicCase);
            u = solveStatic(S);
            for yTip = [-1 1]
                iR = find(abs(ref.xCoor - 10) < 1e-9 & abs(ref.yCoor - yTip) < 1e-9);
                iS = find(abs(S.xCoor - 10) < 1e-9 & abs(S.yCoor - yTip) < 1e-9);
                tc.verifyEqual(u([2*iS-1, 2*iS]), uRef([2*iR-1, 2*iR]), ...
                    'RelTol', 1e-8, sprintf('tip node y = %g', yTip));
            end
        end

        function cubicPatchTestMatchesStoredReference(tc, cubicCase)
            S = loadCase(cubicCase);
            u = solveStatic(S);
            tc.verifyEqual(u, S.u, 'AbsTol', 1e-8*max(abs(S.u)));
        end

        function lowOrderCasesRun(tc, anyCase)
            % p below the degree of the exact field: the analysis must run
            % and give a finite, non-trivial solution.
            S = loadCase(anyCase);
            u = solveStatic(S);
            tc.verifyTrue(all(isfinite(u)));
            tc.verifyGreaterThan(max(abs(u)), 0);
        end
    end
end

function S = loadCase(name)
here = fileparts(mfilename('fullpath'));
S = load(fullfile(fileparts(here), 'fixtures', 'reference', [name '.mat']));
end

function [u, F] = solveStatic(S)
% Same steps as call_Static.m (assembly, Neumann and Dirichlet boundary
% conditions, solution of the reduced system).
[~, K] = evalc(['AssemblyFEM(S.Topology_Node, S.p, S.xCoor, S.yCoor, ' ...
    'S.num_Elem, S.num_DoF, S.C, S.rho, numel(S.C), S.mass_matrix, ' ...
    'S.material_elset, 0, 0, S.zeta, 0)']);
MAT.elset = S.material_elset;
MAT.prop  = S.material_properties;
[~, F] = evalc(['NeumannBCs(S.p, S.cload, S.dsload, S.body_force, ' ...
    'S.gravity, MAT, S.num_DoF, S.Topology_Node, S.xCoor, S.yCoor, ' ...
    'S.num_Nodes)']);
[~, BC_DOFs] = evalc('DirichletBCs_homogeneous(S.boundary_conditions)');
[~, K, F] = evalc(['DirichletBCs_nonhomogeneous(K, F, S.num_DoF, ' ...
    'S.boundary_conditions, S.xCoor, S.yCoor)']);
act = setdiff((1:S.num_DoF)', unique(BC_DOFs));
u = zeros(S.num_DoF, 1);
u(act) = K(act, act)\F(act);
end
