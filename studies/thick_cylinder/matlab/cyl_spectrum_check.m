function cyl_spectrum_check(cases, outfile)
%CYL_SPECTRUM_CHECK Singular value spectra and error decomposition.
%   For each case [alpha N p] (smooth MMS, both displacement components
%   prescribed on x = 0 and y = 0): full SVD of K_aa = U S V' for PM, PP
%   and PH (hybrid trial space).
%   The nodal error vector e = u_h - u_I = -V S^{-1} U' r is decomposed into
%   the contributions of the singular triplets. Reported: number of singular
%   values of PM below the smallest one of PP (would indicate spurious
%   modes), smallest singular values, and the index k (sorted ascending) and
%   share of the triplet that dominates ||e||.
here = fileparts(mfilename('fullpath'));
cyl_setup_paths();
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(1000, 0.25));
prob = cyl_problem('smooth', 1000, 0.25);
L = {};
for i = 1:size(cases, 1)
    a = cases(i, 1);  N = cases(i, 2);  pdeg = cases(i, 3);
    m = cyl_model(fullfile(here, '..', 'mesh', 'inp', sprintf('ThickCylinder_N%d_a%03d.inp', N, a)), pdeg);
    F = cyl_loads(m, prob);
    forms = {'PM', 'PP', 'PH'};  s = cell(1, 3);  info = cell(1, 3);
    for f = 1:3
        try
            K = cyl_assemble(m, C, forms{f});
        catch err
            if ~startsWith(err.identifier, 'ShapeFunctions2D_Metric:'), rethrow(err); end
            s{f} = NaN(60, 1);  info{f} = NaN(1, 4);
            continue
        end
        [~, uEx, act] = cyl_solve(m, K, F, prob);
        fix = setdiff((1:m.num_DoF).', act);
        r = K(act, act)*uEx(act) + K(act, fix)*uEx(fix) - F(act);
        [U, S, V] = svd(full(K(act, act)));
        sv = diag(S);  [sv, o] = sort(sv);  U = U(:, o);  V = V(:, o);
        c = (U.'*r)./sv;                  % coefficients of e in the basis V
        [cmax, k] = max(abs(c));
        info{f} = [sv(1), k, cmax^2/sum(c.^2), sv(k)/sv(end)];
        s{f} = sv;
    end
    nsp = sum(s{1} < s{2}(1));
    L{end+1} = sprintf(['alpha=%2d N=%d p=%d: #sv(PM) < min sv(PP): %d | PM: dominant k=%d ' ...
        '(share %.2f, s_k/s_max=%.2e) | PP: dominant k=%d (share %.2f, s_k/s_max=%.2e)'], ...
        a, N, pdeg, nsp, info{1}(2), info{1}(3), info{1}(4), info{2}(2), info{2}(3), info{2}(4)); %#ok<AGROW>
    fprintf('%s\n', L{end});
    q = [1 2 5 10 20 50];
    L{end+1} = ['   s_k(PM)/s_k(PP), k = 1 2 5 10 20 50: ' sprintf('%.3f ', s{1}(q)./s{2}(q))]; %#ok<AGROW>
    fprintf('%s\n', L{end});
    L{end+1} = sprintf(['   PH: #sv(PH) < min sv(PP): %d, dominant k=%d (share %.2f, s_k/s_max=%.2e); ' ...
        's_k(PH)/s_k(PP), k = 1 2 5 10 20 50: %s'], sum(s{3} < s{2}(1)), info{3}(2), info{3}(3), ...
        info{3}(4), sprintf('%.3f ', s{3}(q)./s{2}(q))); %#ok<AGROW>
    fprintf('%s\n', L{end});
end
fid = fopen(outfile, 'w');  fprintf(fid, '%s\n', L{:});  fclose(fid);
end
