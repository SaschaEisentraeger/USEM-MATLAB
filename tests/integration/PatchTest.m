classdef PatchTest < matlab.unittest.TestCase
    %PATCHTEST Displacement patch tests on distorted meshes.
    %
    %   The defining property of the unsymmetric (Petrov-Galerkin) element
    %   is that it reproduces every displacement field that lies in the
    %   trial space exactly, regardless of the element shape (Rajendran
    %   2010). The metric trial space of degree p contains all polynomials
    %   of total degree <= p, so equilibrium fields of degree <= p must be
    %   reproduced to machine precision on angular and curved-edge
    %   distorted meshes. The exact field is prescribed on the boundary
    %   and the interior solution is compared with it.

    properties (TestParameter)
        p = {1, 2, 3, 4};
        mesh = struct('regular',  {{1, 'none',     0}}, ...
                      'angular',  {{1, 'angular',  0.3}}, ...
                      'curved',   {{2, 'curved',   0.15}}, ...
                      'combined', {{2, 'combined', 0.3}});
        analysis = struct('planeStress', 1, 'planeStrain', 2);
    end

    properties (Constant)
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
        function linearFieldIsExact(tc, mesh, p)
            % Classical patch test (constant strain)
            model = buildModel(mesh, p);
            K = usemAssemble(model, HookeMatrix2d(1, tc.nu, 1));
            field = @(x, y) deal(1e-3 + 2e-3*x - 1e-3*y, -5e-4 + 1e-3*x + 3e-3*y);
            [~, ~, err] = usemDirichletSolve(K, model, field);
            tc.verifyLessThan(err, 1e-10);
        end

        function harmonicFieldsUpToDegreePAreExact(tc, mesh, p, analysis)
            % u = (Re z^k, -Im z^k) is divergence-free and harmonic, hence
            % in equilibrium for any isotropic material without body force.
            model = buildModel(mesh, p);
            K = usemAssemble(model, HookeMatrix2d(1, tc.nu, analysis));
            for k = 1:p
                field = @(x, y) harmonicField(x - 1.3, y - 0.7, k);
                [~, ~, err] = usemDirichletSolve(K, model, field);
                tc.verifyLessThan(err, 1e-9, sprintf('degree k = %d', k));
            end
        end

        function pureBendingIsExactForQuadraticElements(tc, mesh, p)
            % Pure bending in plane stress: u = kappa*(-x*y, (x^2 + nu*y^2)/2)
            % (quadratic, non-zero divergence, depends on nu).
            tc.assumeGreaterThanOrEqual(p, 2);
            model = buildModel(mesh, p);
            K = usemAssemble(model, HookeMatrix2d(1, tc.nu, 1));
            kappa = 1e-3; nu0 = tc.nu;
            field = @(x, y) deal(-kappa*x.*y, kappa/2*(x.^2 + nu0*y.^2));
            [~, ~, err] = usemDirichletSolve(K, model, field);
            tc.verifyLessThan(err, 1e-9);
        end

        function classicalElementIsNotExactOnCurvedMesh(tc)
            % Control experiment: with the same mesh the conventional
            % (symmetric, isoparametric) spectral element loses quadratic
            % completeness, so the mesh is distorted enough to be a
            % meaningful test of the unsymmetric element.
            model = buildModel({2, 'curved', 0.15}, 2);
            C = HookeMatrix2d(1, tc.nu, 1);
            Ksym = assembleSymmetric(model, C);
            field = @(x, y) harmonicField(x - 1.3, y - 0.7, 2);
            [~, ~, errSym] = usemDirichletSolve(Ksym, model, field);
            K = usemAssemble(model, C);
            [~, ~, errUnsym] = usemDirichletSolve(K, model, field);
            tc.verifyGreaterThan(errSym, 1e-4);
            tc.verifyLessThan(errUnsym, 1e-10);
        end

        function fieldsOfHigherDegreeAreNotExact(tc)
            % Sanity check of the test itself: a cubic field is outside
            % the trial space of p = 2.
            model = buildModel({2, 'combined', 0.3}, 2);
            K = usemAssemble(model, HookeMatrix2d(1, tc.nu, 1));
            [~, ~, err] = usemDirichletSolve(K, model, ...
                @(x, y) harmonicField(x, y, 3));
            tc.verifyGreaterThan(err, 1e-5);
        end

        function pRefinementConvergesExponentially(tc)
            % Smooth non-polynomial equilibrium field
            % u = (e^x cos y, -e^x sin y) on a distorted mesh.
            pvec = 2:7;
            err = zeros(size(pvec));
            for k = 1:numel(pvec)
                model = buildModel({2, 'combined', 0.3}, pvec(k));
                K = usemAssemble(model, HookeMatrix2d(1, tc.nu, 1));
                [~, ~, err(k)] = usemDirichletSolve(K, model, ...
                    @(x, y) deal(exp(x).*cos(y), -exp(x).*sin(y)));
            end
            tc.verifyTrue(all(diff(log(err)) < 0), ...
                sprintf('errors not decreasing: %s', mat2str(err, 3)));
            tc.verifyLessThan(err(end), 1e-6);
            tc.verifyLessThan(err(end)/err(1), 1e-4);
        end
    end
end

function model = buildModel(mesh, p)
[order, kind, amount] = mesh{:};
[nodes, topo] = usemTestMesh(3, 2, 3, 2, order);
if ~strcmp(kind, 'none')
    nodes = usemDistortMesh(nodes, topo, kind, amount);
end
model = usemBuildModel(nodes, topo, p);
end

function [ux, uy] = harmonicField(x, y, k)
z = complex(x, y).^k;
ux = real(z);
uy = -imag(z);
end

function K = assembleSymmetric(model, C)
n = model.num_DoF;
K = sparse(n, n);
T = Topology_DoF_2D(model.Topology_Node);
for e = 1:model.num_Elem
    nd = model.Topology_Node(e,:);
    Ke = usemKeSymmetric(model.p, C, model.xCoor(nd), model.yCoor(nd));
    K(T(e,:), T(e,:)) = K(T(e,:), T(e,:)) + Ke; %#ok<SPRIX>
end
end
