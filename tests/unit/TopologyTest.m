classdef TopologyTest < matlab.unittest.TestCase
    %TOPOLOGYTEST Connectivity data and the generation of high-order nodes
    %(Topology_*, ConnectivityEdgeElement, pSpectralElement_NodesTopology).

    properties (TestParameter)
        p = {1, 2, 3, 5};
        order = struct('linear', 1, 'quadratic', 2);
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
        function dofTopologyInterleavesComponents(tc)
            T = Topology_DoF_2D([1 2 5 4; 2 3 6 5]);
            tc.verifyEqual(T, [1 2 3 4 9 10 7 8; 3 4 5 6 11 12 9 10]);
        end

        function edgeTopologyCountsUniqueEdges(tc)
            [~, topo] = usemTestMesh(3, 2, 3, 2, 1);
            [TE, nE] = Topology_Edge_2D(topo);
            tc.verifyEqual(nE, 3*(2+1) + 2*(3+1));      % 17 edges
            tc.verifySize(TE, [6 4]);
            % interior edges are shared by exactly two elements
            [list, localEdge] = ConnectivityEdgeElement(TE);
            counts = cellfun(@numel, list);
            tc.verifyEqual(nnz(counts == 2), 7);
            tc.verifyEqual(nnz(counts == 1), 10);
            for i = 1:nE
                for k = 1:numel(list{i})
                    tc.verifyEqual(TE(list{i}(k), localEdge{i}(k)), i);
                end
            end
        end

        function highOrderMeshHasExpectedNodeCount(tc, p, order)
            nx = 3; ny = 2;
            [nodes, topo] = usemTestMesh(nx, ny, 3, 2, order);
            model = usemBuildModel(nodes, topo, p);
            tc.verifySize(model.Topology_Node, [nx*ny (p+1)^2]);
            tc.verifyEqual(model.num_Nodes, (nx*p+1)*(ny*p+1));
            tc.verifyEqual(numel(unique(model.Topology_Node)), model.num_Nodes);
            % no two nodes at the same position
            XY = round([model.xCoor model.yCoor], 10);
            tc.verifyEqual(size(unique(XY, 'rows'), 1), model.num_Nodes);
        end

        function highOrderNodesLieOnTheGeometry(tc, p, order)
            % Every element node must coincide with the image of its
            % reference position under the ORIGINAL (bilinear or
            % serendipity) geometry map, i.e. neighbouring elements share
            % edge nodes consistently.
            [nodes, topo] = usemTestMesh(2, 2, 2, 2, order);
            nodes = usemDistortMesh(nodes, topo, 'angular', 0.2);
            if order == 2
                nodes = usemDistortMesh(nodes, topo, 'curved', 0.1);
            end
            model = usemBuildModel(nodes, topo, p);
            [xi, eta] = referenceNodes(model.p);
            for e = 1:size(topo, 1)
                if order == 1
                    N = ShapeFunctions2D_4nodes(xi, eta);
                else
                    N = ShapeFunctions2D_8nodes(xi, eta);
                end
                expected = N'*nodes(topo(e,:), :);
                actual = [model.xCoor(model.Topology_Node(e,:)), ...
                          model.yCoor(model.Topology_Node(e,:))];
                if order == 2 && p == 1
                    expected = expected(1:4, :);   % only corners remain
                end
                tc.verifyEqual(actual, expected, 'AbsTol', 1e-12, ...
                    sprintf('element %d', e));
            end
        end

        function cornerNodesMustBeNumberedFirst(tc)
            % 8-node element whose mid-side nodes are numbered between
            % the corner nodes
            x = [0 1 2 2 2 1 0 0]'; y = [0 0 0 1 2 2 2 1]';
            P.pxi = 2; P.peta = 2;
            tc.verifyError(@() pSpectralElement_NodesTopology( ...
                [1 3 5 7 2 4 6 8], x, y, P, 2), ...
                'pSpectralElement_NodesTopology:cornerNodes');
        end

        function linearMeshIsUnchangedForPOne(tc)
            [nodes, topo] = usemTestMesh(2, 2, 2, 2, 1);
            model = usemBuildModel(nodes, topo, 1);
            tc.verifyEqual(model.Topology_Node, topo);
            tc.verifyEqual([model.xCoor model.yCoor], nodes);
        end

    end
end

function [xi, eta] = referenceNodes(P)
x = nodal_distribution(P.pxi);
e = nodal_distribution(P.peta);
xi  = [-1; 1; 1; -1; x(2:end-1); ones(P.peta-1,1); x(2:end-1); -ones(P.peta-1,1)];
eta = [-1; -1; 1; 1; -ones(P.pxi-1,1); e(2:end-1); ones(P.pxi-1,1); e(2:end-1)];
[XI, ETA] = meshgrid(x(2:end-1), e(2:end-1));
XI = XI'; ETA = ETA';
xi  = [xi;  XI(:)];
eta = [eta; ETA(:)];
end
