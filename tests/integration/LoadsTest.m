classdef LoadsTest < matlab.unittest.TestCase
    %LOADSTEST Body forces, gravity and load cases (NeumannBCs,
    %CreateBload, CreateGravity, SpatialLoadingFunction).
    %
    %   2 x 2 mesh of [0,2] x [0,2], p = 2. Elements 1, 2 (bottom row) are
    %   made of material 1, elements 3, 4 (top row) of material 2. The
    %   consistent nodal forces of a constant volume load sum up to the
    %   load times the loaded area.

    properties
        model
        MAT
    end

    methods (TestClassSetup)
        function addCodeToPath(tc)
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'code', 'subroutines'), 'IncludingSubfolders', true));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'tests', 'helpers')));
        end

        function buildModel(tc)
            [nodes, topo] = usemTestMesh(2, 2, 2, 2, 1);
            tc.model = usemBuildModel(nodes, topo, 2);
            tc.MAT.elset = {[1; 2]; [3; 4]};
            tc.MAT.prop  = {1000, [], [70e9 0.3], [], []; ...
                            2000, [], [70e9 0.3], [], []};
        end
    end

    methods (Test)
        function bodyForceOnElementSubset(tc)
            bf = struct('set_name', 'bf', 'amplitude', [3; -1], ...
                        'type', 'Body_force', 'elset', [1; 3], 'step_name', 'Step-1');
            F = neumann(tc, [], bf, [], false);
            tc.verifySize(F, [tc.model.num_DoF 1]);
            % elements 1 and 3: area 2
            tc.verifyEqual(sum(F(1:2:end)), 3*2, 'RelTol', 1e-12);
            tc.verifyEqual(sum(F(2:2:end)), -1*2, 'RelTol', 1e-12);
            % no load on nodes that only belong to elements 2 and 4
            x = tc.model.xCoor;
            tc.verifyEqual(F([2*find(x > 1.5)-1; 2*find(x > 1.5)]), ...
                           zeros(2*nnz(x > 1.5), 1));
        end

        function severalBodyForces(tc)
            bf = struct('set_name', {'a'; 'b'}, 'amplitude', {[1; 0]; [0; 2]}, ...
                        'type', 'Body_force', 'elset', {[1; 2]; [4]}, ...
                        'step_name', 'Step-1');
            F = neumann(tc, [], bf, [], false);
            tc.verifyEqual(sum(F(1:2:end)), 1*2, 'RelTol', 1e-12);
            tc.verifyEqual(sum(F(2:2:end)), 2*1, 'RelTol', 1e-12);
        end

        function gravityUsesTheDensityOfEachMaterial(tc)
            g = struct('set_name', 'whole', 'amplitude', [9.81; 0; -1], ...
                       'type', 'Body_force', 'elset', [], 'step_name', 'Step-1');
            F = neumann(tc, [], [], g, false);
            tc.verifyEqual(sum(F(1:2:end)), 0, 'AbsTol', 1e-9);
            tc.verifyEqual(sum(F(2:2:end)), -9.81*(1000*2 + 2000*2), 'RelTol', 1e-12);
        end

        function multipleRightHandSides(tc)
            bf = struct('set_name', 'bf', 'amplitude', [3; -1], ...
                        'type', 'Body_force', 'elset', [1; 3], 'step_name', 'Step-1');
            g = struct('set_name', 'whole', 'amplitude', [9.81; 0; -1], ...
                       'type', 'Body_force', 'elset', [], 'step_name', 'Step-1');
            F = neumann(tc, [], bf, g, true);
            tc.verifySize(F, [tc.model.num_DoF 2]);
            tc.verifyEqual(F(:,1), neumann(tc, [], bf, [], false), 'AbsTol', 1e-12);
            tc.verifyEqual(F(:,2), neumann(tc, [], [], g, false), 'AbsTol', 1e-9);
        end

        function spatialLoadIsScaledByAmplitude(tc)
            P = [0 0.5; 3 -1];
            tc.verifyEqual(SpatialLoadingFunction([2 1], P), 2*120*P(:,2));
            tc.verifyEqual(SpatialLoadingFunction([2 2], P), 2*(6 - 6*P(:,2).^2));
        end
    end
end

function F = neumann(tc, cload, bf, g, flag_mRHS)
m = tc.model;
[~, F] = evalc(['NeumannBCs(m.p, cload, [], bf, g, tc.MAT, m.num_DoF, ' ...
    'm.Topology_Node, m.xCoor, m.yCoor, m.num_Nodes, flag_mRHS)']);
end
