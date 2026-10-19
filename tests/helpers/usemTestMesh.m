function [nodes, topo] = usemTestMesh(nx, ny, Lx, Ly, order)
%USEMTESTMESH Structured quadrilateral mesh in ABAQUS node ordering.
%
%   [nodes, topo] = usemTestMesh(nx, ny, Lx, Ly, order) returns a mesh of
%   nx-by-ny elements on the rectangle [0,Lx] x [0,Ly].
%
%   order = 1: 4-node elements (CPS4 ordering: corners counter-clockwise)
%   order = 2: 8-node elements (CPS8 ordering: corners 1-4, then the
%              mid-side nodes 5 (edge 1-2), 6 (2-3), 7 (3-4), 8 (4-1))
%
%   nodes : [nNodes x 2] coordinates. The corner nodes come first, as
%           pSpectralElement_NodesTopology expects.
%   topo  : [nElem x 4] or [nElem x 8] connectivity.

if nargin < 5
    order = 1;
end

% Corner nodes (numbered row by row)
[X, Y] = meshgrid(linspace(0, Lx, nx+1), linspace(0, Ly, ny+1));
X = X'; Y = Y';
nodes = [X(:), Y(:)];
cid = @(i, j) (j-1)*(nx+1) + i;          % i = 1..nx+1, j = 1..ny+1

nElem = nx*ny;
topo = zeros(nElem, 4*order);
e = 0;
for j = 1:ny
    for i = 1:nx
        e = e + 1;
        topo(e, 1:4) = [cid(i,j), cid(i+1,j), cid(i+1,j+1), cid(i,j+1)];
    end
end

if order == 1
    return
end

% Mid-side nodes: one per unique edge, appended after the corner nodes
edgeMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
pairs = [1 2; 2 3; 3 4; 4 1];
for e = 1:nElem
    for k = 1:4
        a = topo(e, pairs(k,1));
        b = topo(e, pairs(k,2));
        key = sprintf('%d_%d', min(a,b), max(a,b));
        if ~isKey(edgeMap, key)
            nodes(end+1, :) = (nodes(a,:) + nodes(b,:))/2; %#ok<AGROW>
            edgeMap(key) = size(nodes, 1);
        end
        topo(e, 4+k) = edgeMap(key);
    end
end

end
