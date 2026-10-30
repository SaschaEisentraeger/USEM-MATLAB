function R = probe_hybrid_space(varargin)
%PROBE_HYBRID_SPACE Lebesgue constant of the metric interpolation for the
%trial space Q_p(s,t) and for a hybrid space (outlook, not used by the code).
%
%   R = probe_hybrid_space()                       p = 2:2:16
%   R = probe_hybrid_space('pvec', 2:2:12, 'writeData', false)
%
%   Q_p(s,t):  span{L_a(s) L_b(t) : a, b <= p}   (skew coordinates; the
%              space of ShapeFunctions2D_Metric)
%   hybrid:    span{L_a(s) L_b(t) : a + b <= p}  (= P_p(x,y))
%            + span{L_a(xi) L_b(eta) : a, b <= p, a + b > p}
%              (high-order complement in the parametric coordinates)
%
%   Both spaces have dimension (p+1)^2 and contain P_p(x,y), i.e. both
%   satisfy the completeness condition for distortion immunity. The
%   Lebesgue constant max_x sum_k |M_k(x)| is evaluated on a dense grid in
%   the element; it is independent of the basis. L_a: Legendre polynomials,
%   (s,t) scaled to the bounding box of the nodes as in
%   ShapeFunctions2D_Metric.

ip = inputParser;
ip.addParameter('pvec', 2:2:16);
ip.addParameter('writeData', true);
ip.parse(varargin{:});
opt = ip.Results;

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(fullfile(root, 'code', 'subroutines')));
addpath(fullfile(root, 'tests', 'helpers'));

shapes = {'trapezoid', 'curved', 'general'};
pvec = opt.pvec(:);
R = struct('p', pvec);
for s = 1:numel(shapes)
    LQ = zeros(numel(pvec), 1);
    LH = zeros(numel(pvec), 1);
    for i = 1:numel(pvec)
        p = pvec(i);
        [x, y, P] = usemTestElement(shapes{s}, p);

        % Reference coordinates of the element nodes
        r = nodal_distribution(p, P.ndist);
        [XI, ETA] = meshgrid(r);
        [~, k] = max(ShapeFunctions2D_SEM(XI(:), ETA(:), P), [], 2);
        xiN = XI(k); etaN = ETA(k);

        % Dense evaluation grid
        [Xd, Ed] = meshgrid(linspace(-1, 1, 4*p + 9));
        Xg = Local2GlobalCoordinatesQuad(P, [Xd(:), Ed(:)], [x, y]);

        % Skew coordinates, scaled to the bounding box of the nodes
        [sN, tN] = CorotationalSkewCoordinates(x, y, x, y, P);
        [sQ, tQ] = CorotationalSkewCoordinates(Xg(:,1), Xg(:,2), x, y, P);
        cs = (max(sN) + min(sN))/2; rs = (max(sN) - min(sN))/2;
        ct = (max(tN) + min(tN))/2; rt = (max(tN) - min(tN))/2;
        sN = (sN - cs)/rs; sQ = (sQ - cs)/rs;
        tN = (tN - ct)/rt; tQ = (tQ - ct)/rt;

        L = @(v) PolynomialBasis1D(v, p, 'legendre', P.ndist);
        Sn = L(sN); Tn = L(tN); Sq = L(sQ); Tq = L(tQ);
        Xn = L(xiN); En = L(etaN); Xq = L(Xd(:)); Eq = L(Ed(:));

        n = (p+1)^2;
        VQ = zeros(n); PQ = zeros(n, numel(sQ));
        VH = zeros(n); PH = zeros(n, numel(sQ));
        c = 0;
        for a = 0:p
            for b = 0:p
                c = c + 1;
                VQ(c,:) = Sn(a+1,:).*Tn(b+1,:);
                PQ(c,:) = Sq(a+1,:).*Tq(b+1,:);
                if a + b <= p
                    VH(c,:) = VQ(c,:);
                    PH(c,:) = PQ(c,:);
                else
                    VH(c,:) = Xn(a+1,:).*En(b+1,:);
                    PH(c,:) = Xq(a+1,:).*Eq(b+1,:);
                end
            end
        end
        LQ(i) = max(sum(abs(VQ\PQ), 1));
        LH(i) = max(sum(abs(VH\PH), 1));
        fprintf('%-10s p = %2d   Lebesgue constant  Q_p(s,t): %9.2e   hybrid: %6.2f\n', ...
                shapes{s}, p, LQ(i), LH(i));
    end
    R.(['LebQ_' shapes{s}]) = LQ;
    R.(['LebH_' shapes{s}]) = LH;
end

if opt.writeData
    out = fullfile(root, 'docs', 'figures', 'data');
    if ~exist(out, 'dir'), mkdir(out); end
    names = fieldnames(R);
    fid = fopen(fullfile(out, 'metric-basis-hybrid.dat'), 'w');
    fprintf(fid, '%s ', names{:});
    fprintf(fid, '\n');
    for i = 1:numel(pvec)
        for k = 1:numel(names)
            fprintf(fid, '%.6e ', R.(names{k})(i));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end

end
