function R = ms_element_study(pvec)
%MS_ELEMENT_STUDY Element-level study of the manuscript (Secs. 5.1-5.2).
%
%   R = ms_element_study()          p = 1..16
%
%   For the four test elements of usemTestElement and every degree p:
%     kappa2_<b>, kappa1_<b>  condition numbers of the metric Vandermonde
%                             matrix V of Q_p(s,t) in the basis b
%                             (mon, leg, che, lag) and of the hybrid space
%                             (hyb, Legendre basis)
%     eKron_<b>               max_{l,k} |M_l(x_k) - delta_lk| (computed)
%     LN, LQ, LT              Lebesgue constants: parametric interpolation,
%                             Q_p(s,t), hybrid space T_p
%     idErr                   max |(V^-1)_ld - M_l(s_d)| (Proposition 2,
%                             Lagrange basis)
%     lamRef                  max_k lambda_ref(s_k) (first factor of
%                             kappa_1 in Proposition 2)
%     asym_<b>                ||K_e - K_e^T||_F/||K_e||_F (parallelogram)
%   Writes manuscript/data/element-<shape>.dat.

if nargin < 1, pvec = 1:16; end
[~, dataDir] = ms_setup();

bases  = {'monomial', 'legendre', 'chebyshev', 'lagrange'};
tags   = {'mon', 'leg', 'che', 'lag'};
shapes = {'parallelogram', 'trapezoid', 'curved', 'general'};
C = HookeMatrix2d(1, 0.3, 1);

for s = 1:numel(shapes)
    T = struct('p', pvec(:));
    for q = [tags, {'hyb'}]
        T.(['kappa2_' q{1}]) = nan(numel(pvec), 1);
        T.(['kappa1_' q{1}]) = nan(numel(pvec), 1);
        T.(['eKron_'  q{1}]) = nan(numel(pvec), 1);
        T.(['asym_'   q{1}]) = nan(numel(pvec), 1);
    end
    T.LN = nan(numel(pvec), 1); T.LQ = T.LN; T.LT = T.LN;
    T.idErr = T.LN; T.lamRef = T.LN;
    for ii = 1:numel(pvec)
        p = pvec(ii);
        [x, y, P] = usemTestElement(shapes{s}, p);
        n = numel(x);
        [xiN, etaN] = ReferenceNodes2D(P);
        for b = 1:numel(bases)
            P.metricBasis = bases{b};
            [Nk, ~, ~, V] = ShapeFunctions2D_Metric(x, y, x, y, P, 1);
            T.(['kappa2_' tags{b}])(ii) = cond(V);
            T.(['kappa1_' tags{b}])(ii) = cond(V, 1);
            T.(['eKron_'  tags{b}])(ii) = max(max(abs(Nk - eye(n))));
        end
        [Nk, ~, ~, V] = ShapeFunctions2D_Metric(x, y, x, y, P, 2, [xiN, etaN]);
        T.kappa2_hyb(ii) = cond(V);
        T.kappa1_hyb(ii) = cond(V, 1);
        T.eKron_hyb(ii)  = max(max(abs(Nk - eye(n))));
        [T.LN(ii), T.LQ(ii), T.LT(ii)] = ms_lebesgue(x, y, P);
        [T.idErr(ii), T.lamRef(ii)] = checkProposition2(x, y, P);
        if strcmp(shapes{s}, 'parallelogram')
            for b = [1 4]
                P.metricBasis = bases{b};
                K = Ke_2D(P, C, x, y);
                T.(['asym_' tags{b}])(ii) = norm(K - K', 'fro')/norm(K, 'fro');
            end
            P.metricMethod = 2;
            K = Ke_2D(P, C, x, y);
            T.asym_hyb(ii) = norm(K - K', 'fro')/norm(K, 'fro');
        end
        fprintf(['%-13s p=%2d  k2: mon %8.2e leg %8.2e che %8.2e lag %8.2e hyb %8.2e' ...
                 ' | eK: mon %8.1e lag %8.1e hyb %8.1e | Leb N %5.2f Q %8.2e T %6.2f' ...
                 ' | id %7.1e\n'], shapes{s}, p, T.kappa2_mon(ii), T.kappa2_leg(ii), ...
                T.kappa2_che(ii), T.kappa2_lag(ii), T.kappa2_hyb(ii), T.eKron_mon(ii), ...
                T.eKron_lag(ii), T.eKron_hyb(ii), T.LN(ii), T.LQ(ii), T.LT(ii), T.idErr(ii));
    end
    R.(shapes{s}) = T;
    ms_write_table(fullfile(dataDir, ['element-' shapes{s} '.dat']), T);
end
end

function [err, lamRef] = checkProposition2(x, y, P)
% (V^-1)_ld = M_l(s_d): evaluate the metric functions (Lagrange basis) at
% the points whose scaled skew coordinates are the reference tensor nodes
P.metricBasis = 'lagrange';
[~, ~, ~, V] = ShapeFunctions2D_Metric(x, y, x, y, P, 1);
[sN, tN, J0] = CorotationalSkewCoordinates(x, y, x, y, P);
cs = (max(sN) + min(sN))/2; rs = (max(sN) - min(sN))/2;
ct = (max(tN) + min(tN))/2; rt = (max(tN) - min(tN))/2;
r = nodal_distribution(P.pxi, P.ndist); r = r(:);
q = nodal_distribution(P.peta, P.ndist); q = q(:);
[S, Tt] = meshgrid(r, q);              % c = i*(px+1) + j + 1: s fastest
S = S'; Tt = Tt';
st = [cs + rs*S(:), ct + rt*Tt(:)]';
N0 = ShapeFunctions2D_SEM(0, 0, P);
x0 = [N0'*x; N0'*y];
X = J0'*st + x0;                       % (s,t) = J0^-T (x - x0)
M = ShapeFunctions2D_Metric(X(1,:)', X(2,:)', x, y, P, 1);
err = max(max(abs(M - inv(V))));
% first factor: max_k sum_c |phi_c(s_k)| = ||V||_1
lamRef = norm(V, 1);
end
