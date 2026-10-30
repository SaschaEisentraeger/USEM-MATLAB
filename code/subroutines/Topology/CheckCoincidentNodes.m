function CheckCoincidentNodes(Topology_Node,xCoor,yCoor,pxi,peta)
%CHECKCOINCIDENTNODES Stop if two nodes of a spectral element coincide.
%
% Check for coincident nodes within the elements

% Two nodes of an element are considered coincident if their distance is
% below 1e-10 times the largest distance between the corner nodes.
tol = 1e-10;
for e = 1:size(Topology_Node,1)
    x = xCoor(Topology_Node(e,:));
    y = yCoor(Topology_Node(e,:));
    x = x(:);
    y = y(:);
    h = max(max(hypot(x(1:4) - x(1:4).', y(1:4) - y(1:4).')));
    D = hypot(x - x.', y - y.');
    D(1:numel(x)+1:end) = inf;
    [dmin,k] = min(D(:));
    if dmin < tol*h
        [a,b] = ind2sub(size(D),k);
        error('pSpectralElement_NodesTopology:coincidentNodes', ...
             ['Element %d: the local nodes %d and %d coincide at ' ...
              '(%g, %g). The isoparametric map of the original element ' ...
              'places a node onto another node of the same element; ' ...
              'the element is too strongly distorted for degree ' ...
              '(%d,%d). Change the mesh (e.g. the mid-side nodes) or ' ...
              'the polynomial degree.'], e, a, b, x(a), y(a), ...
              pxi, peta);
    end
end
end
