function T = ms_scaling_split(pvec)
%MS_SCALING_SPLIT Separates the effect of the bounding-box scaling from that
%of the basis: 2-norm condition number of the Q_p(s,t) Vandermonde matrix
%for monomials and Lagrange polynomials, each without and with scaling
%(manuscript, Sec. 6.2). Writes manuscript/data/scaling-split.dat.
if nargin < 1, pvec = [4 8 12 16]; end
[~, dataDir] = ms_setup();
shapes = {'parallelogram', 'trapezoid', 'curved', 'general'};
T = struct('p', pvec(:));
for s = 1:numel(shapes)
    for q = {'mon0', 'mon1', 'lag0', 'lag1'}
        T.([q{1} '_' shapes{s}]) = nan(numel(pvec), 1);
    end
    for i = 1:numel(pvec)
        [x, y, P] = usemTestElement(shapes{s}, pvec(i));
        [sN, tN] = CorotationalSkewCoordinates(x, y, x, y, P);
        cs = (max(sN) + min(sN))/2; rs = (max(sN) - min(sN))/2;
        ct = (max(tN) + min(tN))/2; rt = (max(tN) - min(tN))/2;
        S = {sN, tN; (sN - cs)/rs, (tN - ct)/rt};
        B = {'monomial', 'lagrange'}; tag = {'mon', 'lag'};
        for b = 1:2
            for k = 1:2
                Vs = PolynomialBasis1D(S{k,1}, P.pxi, B{b});
                Vt = PolynomialBasis1D(S{k,2}, P.peta, B{b});
                V = zeros(size(Vs,1)*size(Vt,1), numel(sN)); c = 0;
                for bb = 1:size(Vt,1)
                    for aa = 1:size(Vs,1)
                        c = c + 1; V(c,:) = Vs(aa,:).*Vt(bb,:);
                    end
                end
                T.(sprintf('%s%d_%s', tag{b}, k-1, shapes{s}))(i) = cond(V);
            end
        end
        fprintf('%-13s p=%2d  mon: %9.2e -> %9.2e (scaled) | lag: %9.2e -> %9.2e (scaled)\n', shapes{s}, pvec(i), ...
            T.(['mon0_' shapes{s}])(i), T.(['mon1_' shapes{s}])(i), T.(['lag0_' shapes{s}])(i), T.(['lag1_' shapes{s}])(i));
    end
end
ms_write_table(fullfile(dataDir, 'scaling-split.dat'), T);
end
