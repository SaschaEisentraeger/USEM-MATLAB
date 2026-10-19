classdef AbaqusInputTest < matlab.unittest.TestCase
    %ABAQUSINPUTTEST Reading the example ABAQUS input files end to end.
    %
    %   The input files in code/InputFilesAbaqus are read with the same
    %   steps as in Start_Sim.m (usemReadAbaqus). The resulting model -
    %   high-order mesh, edge loads and boundary conditions adapted to the
    %   spectral elements - must equal the data stored in the reference
    %   workspaces (tests/fixtures/reference), and the static analysis of
    %   the quadratic patch test must reproduce pure bending exactly.

    properties (TestParameter)
        refCase = {'PatchTestLeeQuadratic_undistorted_ndof12', ...
                   'PatchTestLeeQuadratic_undistorted_ndof30', ...
                   'PatchTestLeeQuadratic_angular_ndof30', ...
                   'PatchTestLeeQuadratic_angular_ndof90', ...
                   'PatchTestLeeQuadratic_midedge_ndof30', ...
                   'PatchTestLeeCubic_undistorted_ndof56', ...
                   'PatchTestLeeCubic_angular_ndof56', ...
                   'PatchTestLeeCubic_midedge_ndof30', ...
                   'PatchTestLeeCubic_midedge_ndof56', ...
                   'PatchTestLeeCubic_Rajendran_adj_ndof104'};
    end

    properties
        rootDir
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            tc.rootDir = root;
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'code', 'subroutines'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'tests', 'helpers')));
        end
    end

    methods (TestMethodSetup)
        function useTemporaryFolder(tc)
            % Adjust_Set_names_input_file_ABAQUS writes to the current folder
            tc.applyFixture(matlab.unittest.fixtures.WorkingFolderFixture);
        end
    end

    methods (Test)
        function inputFileReproducesReferenceModel(tc, refCase)
            R = loadRef(tc, refCase);
            S = readCase(tc, refCase, R.p);
            tc.verifyEqual(S.nodes, R.nodes);
            tc.verifyEqual(double(S.topology), double(R.topology));
            tc.verifyEqual(S.Topology_Node, R.Topology_Node);
            tc.verifyEqual([S.xCoor S.yCoor], [R.xCoor R.yCoor], 'AbsTol', 1e-12);
            tc.verifyEqual(S.material_properties, R.material_properties);
            tc.verifyEqual(S.step_type, {'Static'});
            tc.verifyEmpty(S.cload);
            tc.verifyNumElements(S.dsload, numel(R.dsload));
            for k = 1:numel(R.dsload)
                tc.verifyEqual(S.dsload(k).type, R.dsload(k).type);
                tc.verifyEqual(S.dsload(k).amplitude, R.dsload(k).amplitude);
                tc.verifyEqual(S.dsload(k).DirectionVector, R.dsload(k).DirectionVector);
                tc.verifyEqual(sort(double(S.dsload(k).nset(:))), ...
                               sort(double(R.dsload(k).nset(:))));
            end
            tc.verifyNumElements(S.boundary_conditions, numel(R.boundary_conditions));
            for k = 1:numel(R.boundary_conditions)
                tc.verifyEqual(S.boundary_conditions(k).type, R.boundary_conditions(k).type);
                tc.verifyEqual(sort(double(S.boundary_conditions(k).nset(:))), ...
                               sort(double(R.boundary_conditions(k).nset(:))));
            end
        end

        function quadraticPatchTestFromInputFileIsExact(tc)
            % p = 3 on the curved-edge mesh: read, assemble, solve.
            p.pxi = 3; p.peta = 3; p.ndist = 2;
            S = readCase(tc, 'PatchTestLeeQuadratic_midedge', p);
            C = PlaneStressConstitutiveMatrices(HookeMatrix3d(7e10, 0.3));
            model = struct('p', p, 'Topology_Node', S.Topology_Node, ...
                'xCoor', S.xCoor, 'yCoor', S.yCoor, 'num_Elem', S.num_Elem, ...
                'num_DoF', S.num_DoF);
            K = usemAssemble(model, C);
            MAT.elset = S.material_elset; MAT.prop = S.material_properties;
            [~, F] = evalc(['NeumannBCs(p, S.cload, S.dsload, S.body_force, ' ...
                'S.gravity, MAT, S.num_DoF, S.Topology_Node, S.xCoor, ' ...
                'S.yCoor, S.num_Nodes)']);
            [~, fixed] = evalc('DirichletBCs_homogeneous(S.boundary_conditions)');
            act = setdiff((1:S.num_DoF)', unique(fixed));
            u = zeros(S.num_DoF, 1);
            u(act) = K(act, act)\F(act);
            right = abs(S.xCoor - 10) < 1e-9;
            a = 1.5*sum(F(2*find(right)-1).*S.yCoor(right));
            x = S.xCoor; y = S.yCoor; E = 7e10; nu = 0.3;
            uy = -a/(2*E)*(x.^2 + nu*y.^2) + a*nu/(2*E);
            tc.verifyEqual(u(1:2:end), a*x.*y/E, 'AbsTol', 1e-9*max(abs(uy)));
            tc.verifyEqual(u(2:2:end), uy, 'AbsTol', 1e-9*max(abs(uy)));
        end
    end
end

function R = loadRef(tc, name)
R = load(fullfile(tc.rootDir, 'tests', 'fixtures', 'reference', [name '.mat']));
end

function S = readCase(tc, name, p)
base = regexprep(name, '_ndof\d+$', '');
inp = fullfile(tc.rootDir, 'code', 'InputFilesAbaqus', [base '.inp']);
[~, S] = evalc('usemReadAbaqus(inp, p)');
end
