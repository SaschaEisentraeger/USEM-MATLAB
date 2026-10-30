function ms_deviation_check(pvec)
%MS_DEVIATION_CHECK Deviation d_K of the unsymmetric stiffness matrices from
%the classical one (Table "deviation", ms_perturbation) with calibrated
%Gauss rules.
%   For the distorted test elements and p in pvec, the rule of every
%   formulation is calibrated on the single element (ms_calibrate), and
%   d_K = ||K_e - K_e^PP||_F/||K_e^PP||_F is computed with these rules and
%   with the former default rule (2p points, as in ms_perturbation).
%   Writes manuscript/data/quadrature/deviation-check.dat.
if nargin < 1, pvec = [4 8 12]; end
[~, dataDir] = ms_setup();
C = HookeMatrix2d(1, 0.3, 1);
shapes = {'trapezoid', 'curved', 'general'};
rows = [];
for s = 1:numel(shapes)
    for p = pvec
        [x, y, P] = usemTestElement(shapes{s}, p);
        nq = zeros(1, 3);  K2 = cell(1, 3);  Kc = cell(1, 3);
        forms = {'PP', 'PM', 'PH'};
        for f = 1:3
            R = ms_calibrate({x, y}, P, C, forms{f});
            nq(f) = R.nG;
            Kc{f} = ms_Ke(P, C, x, y, nq(f), forms{f});
            K2{f} = ms_Ke(P, C, x, y, 2*p, forms{f});
        end
        d = @(K, f) norm(K{f} - K{1}, 'fro')/norm(K{1}, 'fro');
        row = [s, p, nq, d(K2, 2), d(K2, 3), d(Kc, 2), d(Kc, 3)];
        rows = [rows; row]; %#ok<AGROW>
        fprintf('%-9s p=%2d  nG %2d %2d %2d | dK(2p): PM %8.2e PH %8.2e | dK(cal): PM %8.2e PH %8.2e\n', ...
                shapes{s}, p, nq, row(6:9));
    end
end
names = {'shape', 'p', 'nG_PP', 'nG_PM', 'nG_PH', 'dK2p_PM', 'dK2p_PH', 'dK_PM', 'dK_PH'};
T = struct();
for k = 1:numel(names), T.(names{k}) = rows(:, k); end
ms_write_table(fullfile(dataDir, 'quadrature', 'deviation-check.dat'), T);
end
