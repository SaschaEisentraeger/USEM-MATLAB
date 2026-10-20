function R = ms_distortion_sweep(pvec)
%MS_DISTORTION_SWEEP Lebesgue constants and conditioning versus the amount
%of distortion of a single element (manuscript, Sec. 5.2).
%
%   Element families on the square [-1,1]^2 (8-node input, lifted to p):
%     'bend'    : mid-side nodes of the bottom and top edge moved by +delta
%                 in y; x = xi, y = eta + delta (1 - xi^2), det J = 1
%     'onecurve': mid-side node of the bottom edge moved by +delta in y
%     'trapezoid': top corners moved inwards by delta (bilinear map)
%   For p in pvec: LN, LQ, LT (Lebesgue constants), kQ, kT (2-norm
%   condition numbers of the Vandermonde matrices, Lagrange basis for
%   Q_p(s,t), Legendre basis for T_p), minDetJ.
%   Writes manuscript/data/sweep-<family>.dat.

if nargin < 1, pvec = [4 8 12]; end
[~, dataDir] = ms_setup();
fam = {'bend', 'onecurve', 'trapezoid'};
dmax = [0.5, 1.0, 0.9];
for f = 1:numel(fam)
    d = linspace(0, dmax(f), 21)';
    T = struct('delta', d);
    for p = pvec
        for q = {'LN', 'LQ', 'LT', 'kQ', 'kT', 'minDetJ'}
            T.(sprintf('%s_p%d', q{1}, p)) = nan(numel(d), 1);
        end
        for i = 1:numel(d)
            nodes = familyNodes(fam{f}, d(i));
            model = usemBuildModel(nodes, 1:size(nodes, 1), p);
            P = model.p;
            x = model.xCoor(model.Topology_Node(1,:));
            y = model.yCoor(model.Topology_Node(1,:));
            [xiN, etaN] = ReferenceNodes2D(P);
            try
                [LN, LQ, LT] = ms_lebesgue(x, y, P);
            catch
                LN = NaN; LQ = NaN; LT = NaN;
            end
            P.metricBasis = 'lagrange';
            try
                [~, ~, ~, V] = ShapeFunctions2D_Metric(x, y, x, y, P, 1);
                kQ = cond(V);
            catch
                kQ = NaN;
            end
            [~, ~, ~, V] = ShapeFunctions2D_Metric(x, y, x, y, P, 2, [xiN, etaN]);
            kT = cond(V);
            [g, ~] = lgwt(2*p, -1, 1);
            [Xi, Eta] = meshgrid(g);
            dj = zeros(numel(Xi), 1);
            for k = 1:numel(Xi)
                dj(k) = det(Jacobi_Matrix(P, [Xi(k), Eta(k)], x, y));
            end
            T.(sprintf('LN_p%d', p))(i) = LN;
            T.(sprintf('LQ_p%d', p))(i) = LQ;
            T.(sprintf('LT_p%d', p))(i) = LT;
            T.(sprintf('kQ_p%d', p))(i) = kQ;
            T.(sprintf('kT_p%d', p))(i) = kT;
            T.(sprintf('minDetJ_p%d', p))(i) = min(dj);
            fprintf('%-9s p=%2d delta=%5.3f  LN %5.2f  LQ %9.2e  LT %7.2f  kQ %9.2e  kT %9.2e\n', ...
                    fam{f}, p, d(i), LN, LQ, LT, kQ, kT);
        end
    end
    R.(fam{f}) = T;
    ms_write_table(fullfile(dataDir, ['sweep-' fam{f} '.dat']), T);
end
end

function nodes = familyNodes(name, d)
nodes = [-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0];
switch name
    case 'bend'
        nodes(5, 2) = nodes(5, 2) + d;
        nodes(7, 2) = nodes(7, 2) + d;
    case 'onecurve'
        nodes(5, 2) = nodes(5, 2) + d;
    case 'trapezoid'
        nodes(3, 1) = 1 - d;
        nodes(4, 1) = -1 + d;
        nodes(6, :) = (nodes(2, :) + nodes(3, :))/2;
        nodes(7, :) = (nodes(3, :) + nodes(4, :))/2;
        nodes(8, :) = (nodes(4, :) + nodes(1, :))/2;
end
end
