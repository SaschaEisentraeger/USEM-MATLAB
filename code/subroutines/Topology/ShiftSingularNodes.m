function [xCoor,yCoor,shifted] = ShiftSingularNodes(Topology_Node,xCoor,yCoor,p)
%SHIFTSINGULARNODES Move single nodes of elements with a singular metric
%interpolation matrix.
%
%   [xCoor, yCoor, shifted] = ShiftSingularNodes(Topology_Node, xCoor,
%   yCoor, p) checks the interpolation matrix V of the metric (trial)
%   shape functions of every spectral element. If V is singular or nearly
%   singular (rcond(V) < p.nodeShiftTol, default 1e-12), the nodes are not
%   unisolvent for the trial space: this happens for coincident nodes or
%   for special node positions on strongly distorted elements. One node of
%   the element is then moved by p.nodeShift times the element size h
%   (largest distance between the corner nodes), which makes V regular.
%
%   Which node: the corner nodes are never moved (they define the input
%   mesh). Interior nodes are preferred, because they belong to this
%   element only and the boundary of the domain stays exact; edge nodes
%   (shared with the neighbour) are used only if the element has no
%   interior nodes. Every candidate is moved on trial, and the one that
%   gives the best-conditioned V is kept. The node moves towards the
%   centroid of the element, rotated by 30 degrees (a fixed, non-symmetric
%   direction, so that a symmetric node configuration does not stay
%   symmetric). If V is still nearly singular, the step is doubled (at most
%   six times).
%
%   The remedy is pragmatic: the geometry of the element changes by about
%   p.nodeShift*h near the moved node, and the metric shape functions
%   interpolate at the moved node. A slight perturbation of the nodal
%   coordinates is one of the remedies for singular interpolation matrices
%   of point interpolation methods (Liu and Gu) that Ooi et al. (2008)
%   discuss for unsymmetric elements; their own remedy, a rotation of the
%   coordinate system, has no effect on metric functions in corotational
%   skew coordinates, which are frame invariant.
%
%   shifted: [nShift x 3] array [element, local node, step size].

if isfield(p,'nodeShiftTol')
    tol = p.nodeShiftTol;
else
    tol = 1e-12;
end
delta = p.nodeShift;
maxDoubling = 6;

% Local node classes: corners, edge nodes, interior nodes
[xiRef,etaRef] = ReferenceNodes2D(p);
interior = find(abs(xiRef) < 1 & abs(etaRef) < 1);
edge = setdiff(5:numel(xiRef),interior);
if ~isempty(interior)
    candidates = interior(:)';
else
    candidates = edge(:)';
end

% Fixed direction: towards the centroid, rotated by 30 degrees
c = cosd(30);
s = sind(30);

shifted = zeros(0,3);
for e = 1:size(Topology_Node,1)
    nodes = Topology_Node(e,:);
    x = xCoor(nodes);
    y = yCoor(nodes);
    if MetricInterpolationCondition(x,y,p) >= tol
        continue
    end
    h = max(max(hypot(x(1:4) - x(1:4).', y(1:4) - y(1:4).')));
    xc = mean(x(1:4));
    yc = mean(y(1:4));
    step = delta*h;
    best = 0;
    for k = 0:maxDoubling
        for a = candidates
            d = [xc - x(a); yc - y(a)];
            if norm(d) == 0
                d = [1; 0];
            end
            d = [c -s; s c]*d/norm(d);
            xt = x; yt = y;
            xt(a) = x(a) + step*d(1);
            yt(a) = y(a) + step*d(2);
            rc = MetricInterpolationCondition(xt,yt,p);
            if rc > best
                best = rc;
                bestNode = a;
                bestXY = [xt(a), yt(a)];
            end
        end
        if best >= tol
            break
        end
        step = 2*step;
    end
    if best < tol
        error('ShiftSingularNodes:failed', ['Element %d: the metric ' ...
              'interpolation matrix stays singular after moving a node ' ...
              'by up to %g times the element size.'],e,step/h);
    end
    xCoor(nodes(bestNode)) = bestXY(1);
    yCoor(nodes(bestNode)) = bestXY(2);
    shifted(end+1,:) = [e, bestNode, step]; %#ok<AGROW>
    warning('ShiftSingularNodes:shifted', ['Element %d: local node %d ' ...
            'moved by %.3g (%.3g times the element size) to make the ' ...
            'metric interpolation matrix regular (rcond = %.2e).'], ...
            e,bestNode,step,step/h,best);
end

% A moved edge node also belongs to the neighbouring element: check again
if isempty(interior) && ~isempty(shifted)
    for e = 1:size(Topology_Node,1)
        nodes = Topology_Node(e,:);
        if MetricInterpolationCondition(xCoor(nodes),yCoor(nodes),p) < tol
            error('ShiftSingularNodes:neighbour', ['Element %d: moving ' ...
                  'an edge node of its neighbour made its metric ' ...
                  'interpolation matrix singular.'],e);
        end
    end
end

end
