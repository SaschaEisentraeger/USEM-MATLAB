function nq = cyl_rule(nqfun, alpha, N, form, p)
%CYL_RULE Gauss points per direction of the stiffness matrices: 2p (default
%of Ke_2D) or nqfun(alpha, N, form, p), e.g. the calibrated rules of the
%manuscript (manuscript/studies/ms_rule). A missing rule (trial space
%rejected during the calibration) falls back to 2p, so that the rejection
%is reproduced by the assembly.
nq = 2*p;
if ~isempty(nqfun)
    n = nqfun(alpha, N, form, p);
    if ~isnan(n), nq = n; end
end
end
