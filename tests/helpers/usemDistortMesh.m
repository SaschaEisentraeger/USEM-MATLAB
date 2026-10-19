function nodes = usemDistortMesh(nodes, topo, kind, amount)
%USEMDISTORTMESH Deterministic distortion of interior nodes.
%
%   nodes = usemDistortMesh(nodes, topo, kind, amount) moves only nodes
%   that are not on the boundary of the rectangle, so the domain and the
%   boundary conditions stay the same.
%
%   kind = 'angular'  : interior corner nodes are shifted in x by
%                       +/- amount (alternating), so straight-sided
%                       elements become trapezoids (Lee & Bathe 1993).
%                       Mid-side nodes of 8-node elements are moved to
%                       the new edge midpoints, so the edges stay straight.
%   kind = 'curved'   : interior mid-side nodes of 8-node elements are
%                       shifted perpendicular to their edge by amount,
%                       so the edges become parabolic (curved-edge
%                       distortion).
%   kind = 'combined' : 'angular' followed by 'curved'.

xmin = min(nodes(:,1)); xmax = max(nodes(:,1));
ymin = min(nodes(:,2)); ymax = max(nodes(:,2));
tol = 1e-12*max(xmax - xmin, ymax - ymin);
onBoundary = @(P) abs(P(:,1) - xmin) < tol | abs(P(:,1) - xmax) < tol | ...
                  abs(P(:,2) - ymin) < tol | abs(P(:,2) - ymax) < tol;

corner = unique(topo(:, 1:4));
pairs = [1 2; 2 3; 3 4; 4 1];

switch kind
    case 'angular'
        interior = corner(~onBoundary(nodes(corner,:)));
        sgn = (-1).^(1:numel(interior))';
        nodes(interior, 1) = nodes(interior, 1) + amount*sgn;
        if size(topo, 2) == 8
            for e = 1:size(topo, 1)
                for k = 1:4
                    m = topo(e, 4+k);
                    nodes(m, :) = (nodes(topo(e,pairs(k,1)),:) + ...
                                   nodes(topo(e,pairs(k,2)),:))/2;
                end
            end
        end

    case 'curved'
        assert(size(topo, 2) == 8, ...
               'Curved-edge distortion requires 8-node elements.');
        done = false(size(nodes, 1), 1);
        for e = 1:size(topo, 1)
            for k = 1:4
                m = topo(e, 4+k);
                if done(m) || onBoundary(nodes(m,:))
                    continue
                end
                t = nodes(topo(e,pairs(k,2)),:) - nodes(topo(e,pairs(k,1)),:);
                n = [-t(2), t(1)]/norm(t);
                nodes(m, :) = nodes(m, :) + amount*n;
                done(m) = true;
            end
        end

    case 'combined'
        nodes = usemDistortMesh(nodes, topo, 'angular', amount);
        nodes = usemDistortMesh(nodes, topo, 'curved', amount/2);

    otherwise
        error('usemDistortMesh:kind', 'Unknown distortion type "%s".', kind);
end

end
