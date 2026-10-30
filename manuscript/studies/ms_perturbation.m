function R = ms_perturbation(pvec)
%MS_PERTURBATION Deviation of the metric element from the classical one,
%and effect of the bounding-box scaling (manuscript, Secs. 5.1, 4).
%
%   For the distorted test elements of usemTestElement and p in pvec:
%     dN_PM, dN_PH  max_k max_q |M_k(x_q) - N_k(x_q)| at the (2p)^2 Gauss
%                   points (shape function deviation, Proposition 4)
%     dK_PM, dK_PH  ||K_e - K_e^PP||_F / ||K_e^PP||_F
%     kNoBox_lag, kBox_lag, kNoBox_leg, kBox_leg
%                   2-norm condition number of V for Q_p(s,t) without and
%                   with the bounding-box scaling (Lagrange, Legendre)
%   Writes manuscript/data/perturbation-<shape>.dat.

if nargin < 1, pvec = 1:12; end
[~, dataDir] = ms_setup();
shapes = {'trapezoid', 'curved', 'general'};
C = HookeMatrix2d(1, 0.3);
for s = 1:numel(shapes)
    T = struct('p', pvec(:));
    for q = {'dN_PM', 'dN_PH', 'dK_PM', 'dK_PH', 'kNoBox_lag', 'kBox_lag', ...
             'kNoBox_leg', 'kBox_leg'}
        T.(q{1}) = nan(numel(pvec), 1);
    end
    for i = 1:numel(pvec)
        p = pvec(i);
        [x, y, P] = usemTestElement(shapes{s}, p);
        [g, ~] = lgwt(2*p, -1, 1);
        [Xi, Eta] = meshgrid(g);
        IP = [Xi(:), Eta(:)];
        Xg = Local2GlobalCoordinatesQuad(P, IP, [x, y]);
        N = ShapeFunctions2D_SEM(IP(:,1), IP(:,2), P);
        MQ = ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, P, 1);
        MT = ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, P, 2, IP);
        T.dN_PM(i) = max(abs(MQ(:) - N(:)));
        T.dN_PH(i) = max(abs(MT(:) - N(:)));
        KP = usemKeSymmetric(P, C, x, y);
        P1 = P; P1.metricMethod = 1;
        P2 = P; P2.metricMethod = 2;
        KQ = Ke_2D(P1, C, x, y);
        KT = Ke_2D(P2, C, x, y);
        T.dK_PM(i) = norm(KQ - KP, 'fro')/norm(KP, 'fro');
        T.dK_PH(i) = norm(KT - KP, 'fro')/norm(KP, 'fro');
        [T.kNoBox_lag(i), T.kBox_lag(i)] = scalingEffect(x, y, P, 'lagrange');
        [T.kNoBox_leg(i), T.kBox_leg(i)] = scalingEffect(x, y, P, 'legendre');
        fprintf(['%-9s p=%2d  dN: PM %8.2e PH %8.2e | dK: PM %8.2e PH %8.2e | ' ...
                 'k lag %8.2e -> %8.2e  leg %8.2e -> %8.2e\n'], ...
                shapes{s}, p, T.dN_PM(i), T.dN_PH(i), T.dK_PM(i), T.dK_PH(i), ...
                T.kNoBox_lag(i), T.kBox_lag(i), T.kNoBox_leg(i), T.kBox_leg(i));
    end
    ms_write_table(fullfile(dataDir, ['perturbation-' shapes{s} '.dat']), T);
    R.(shapes{s}) = T;
end
end

function [k0, k1] = scalingEffect(x, y, P, basis)
% Condition number of the Q_p(s,t) Vandermonde matrix in unscaled skew
% coordinates (k0) and in bounding-box scaled ones (k1, as in the code)
[sN, tN] = CorotationalSkewCoordinates(x, y, x, y, P);
V0 = vdm(sN, tN, P, basis);
cs = (max(sN) + min(sN))/2; rs = (max(sN) - min(sN))/2;
ct = (max(tN) + min(tN))/2; rt = (max(tN) - min(tN))/2;
V1 = vdm((sN - cs)/rs, (tN - ct)/rt, P, basis);
k0 = cond(V0);
k1 = cond(V1);
end

function V = vdm(s, t, P, basis)
Vs = PolynomialBasis1D(s, P.pxi, basis);
Vt = PolynomialBasis1D(t, P.peta, basis);
V = zeros(size(Vs, 1)*size(Vt, 1), numel(s));
c = 0;
for b = 1:size(Vt, 1)
    for a = 1:size(Vs, 1)
        c = c + 1;
        V(c, :) = Vs(a, :).*Vt(b, :);
    end
end
end
