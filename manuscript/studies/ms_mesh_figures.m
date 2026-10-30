function ms_mesh_figures()
%MS_MESH_FIGURES Geometry data of the figures of the test elements and meshes
%(manuscript, figures fig-elements, fig-mesh32 and fig-square).
%
%   Writes to manuscript/data (pgfplots tables, columns "x y", polylines
%   separated by blank lines):
%     shape-<shape>-outline.dat : boundary of the single test elements
%                                 (parallelogram, trapezoid, curved, general)
%     shape-<shape>-p8.dat      : their GLL nodes for p = 8 (usemTestElement)
%     mesh32-<mesh>.dat         : element edges of the 3 x 2 meshes on
%                                 [0,3] x [0,2] (ms_spectrum): curved15,
%                                 curved20, combined30
%     square-d<100 delta>-N4.dat: element edges of the 4 x 4 meshes of the
%                                 h-convergence study (ms_convergence),
%                                 delta = 0.1, 0.2, 0.3
%     square-d<100 delta>-N4-negJ.dat : closed polygons bounding the part of
%                                 the elements with det(J) < 0 of the 8-node
%                                 map (only written if this part is not empty)
%   The edges are the parabolas of the 8-node geometry (25 points per edge).

[~, dataDir] = ms_setup();

% single test elements, p = 8
for s = {'parallelogram', 'trapezoid', 'curved', 'general'}
    [x, y, P] = usemTestElement(s{1}, 8);
    r = linspace(-1, 1, 41)';
    o = ones(size(r));
    B = [r, -o; o, r; flipud(r), o; -o, flipud(r)];
    Xb = Local2GlobalCoordinatesQuad(P, B, [x, y]);
    writePolylines(fullfile(dataDir, ['shape-' s{1} '-outline.dat']), {Xb});
    writePolylines(fullfile(dataDir, ['shape-' s{1} '-p8.dat']), {[x(:), y(:)]});
end

% 3 x 2 meshes of the patch test and the eigenvalue study (ms_spectrum)
[nodes0, topo] = usemTestMesh(3, 2, 3, 2, 2);
meshes = {'curved15', 'curved', 0.15; 'curved20', 'curved', 0.20;
          'combined30', 'combined', 0.30};
for m = 1:size(meshes, 1)
    nodes = usemDistortMesh(nodes0, topo, meshes{m, 2}, meshes{m, 3});
    writePolylines(fullfile(dataDir, ['mesh32-' meshes{m, 1} '.dat']), meshEdges(nodes, topo));
    fprintf('mesh32-%-10s min det J (8-node map) = %.4f\n', meshes{m, 1}, minDet(nodes, topo));
end

% unit square, N = 4, persistent curved-edge distortion (ms_convergence)
N = 4;
for delta = [0.1 0.2 0.3]
    [nodes, topo] = usemTestMesh(N, N, 1, 1, 2);
    nodes = usemDistortMesh(nodes, topo, 'curved', delta/N);
    tag = sprintf('square-d%03d-N%d', round(100*delta), N);
    writePolylines(fullfile(dataDir, [tag '.dat']), meshEdges(nodes, topo));
    neg = negativeRegion(nodes, topo);
    if ~isempty(neg)
        writePolylines(fullfile(dataDir, [tag '-negJ.dat']), neg);
    end
    fprintf('%s min det J (8-node map) = %.4f, %d polygon(s) with det J < 0\n', ...
            tag, minDet(nodes, topo), numel(neg));
end
end

function L = meshEdges(nodes, topo)
% Each edge once, as the parabola through its three nodes
pairs = [1 2; 2 3; 3 4; 4 1];
s = linspace(-1, 1, 25)';
L = {};
seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for e = 1:size(topo, 1)
    for k = 1:4
        a = topo(e, pairs(k, 1));  b = topo(e, pairs(k, 2));  c = topo(e, 4 + k);
        key = sprintf('%d_%d', min(a, b), max(a, b));
        if isKey(seen, key), continue, end
        seen(key) = true;
        L{end+1} = (s.*(s - 1)/2)*nodes(a, :) + (1 - s.^2)*nodes(c, :) + ...
                   (s.*(s + 1)/2)*nodes(b, :); %#ok<AGROW>
    end
end
end

function d = minDet(nodes, topo)
g = linspace(-1, 1, 201);
[G1, G2] = meshgrid(g, g);
d = inf;
for e = 1:size(topo, 1)
    d = min(d, min(detJ8(nodes(topo(e, :), :), G1, G2), [], 'all'));
end
end

function D = detJ8(X, G1, G2)
% det(J) of the 8-node map with nodes X (8 x 2) at the points (G1, G2)
[~, dNxi, dNeta] = ShapeFunctions2D_8nodes(G1(:), G2(:));
xx = dNxi.'*X(:, 1);  xy = dNxi.'*X(:, 2);
yx = dNeta.'*X(:, 1); yy = dNeta.'*X(:, 2);
D = reshape(xx.*yy - xy.*yx, size(G1));
end

function L = negativeRegion(nodes, topo)
% Closed polygons (physical coordinates) around {det J < 0} of every element:
% zero contour of det J on a fine reference grid, padded with positive
% values so that the contours close along the element boundary
g = linspace(-1, 1, 401);
gp = [-1 - 1e-9, g, 1 + 1e-9];
[G1, G2] = meshgrid(g, g);
L = {};
for e = 1:size(topo, 1)
    X = nodes(topo(e, :), :);
    D = detJ8(X, G1, G2);
    if all(D(:) > 0), continue, end
    Dp = ones(numel(gp));
    Dp(2:end-1, 2:end-1) = D;
    C = contourc(gp, gp, Dp, [0 0]);
    k = 1;
    while k < size(C, 2)
        n = C(2, k);
        xi = min(max(C(1, k+1:k+n), -1), 1);
        eta = min(max(C(2, k+1:k+n), -1), 1);
        N8 = ShapeFunctions2D_8nodes(xi(:), eta(:));
        L{end+1} = [N8.'*X(:, 1), N8.'*X(:, 2)]; %#ok<AGROW>
        k = k + n + 1;
    end
end
end

function writePolylines(file, L)
fid = fopen(file, 'w');
fprintf(fid, 'x y\n');
for k = 1:numel(L)
    if k > 1, fprintf(fid, '\n'); end
    fprintf(fid, '%.6f %.6f\n', L{k}.');
end
fclose(fid);
end
