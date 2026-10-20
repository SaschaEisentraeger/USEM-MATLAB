function cyl_conditioning_check(outfile)
%CYL_CONDITIONING_CHECK Consistency vs. stability on selected meshes.
%   For the polynomial field (in the PM trial space for p >= 5) the PM
%   stiffness is integrated exactly (nq = p^2 + p/2 rounded up, see
%   run_quadrature_study). Reported per case:
%     res   = ||K_aa u_ex,a + K_ab u_b - F_a|| / ||F_a||  (consistency residual
%             of the exact nodal values; zero in exact arithmetic)
%     eN    = nodal error of the computed solution
%     smin  = smallest singular value of K_aa (scaled by the largest),
%             PM and PP
here = fileparts(mfilename('fullpath'));
if nargin < 1, outfile = fullfile(here, '..', 'results', 'conditioning.txt'); end
cyl_setup_paths();
inpDir = fullfile(here, '..', 'mesh', 'inp');
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(1000, 0.25));
prob = cyl_problem('poly', 1000, 0.25);
L = {sprintf('%-28s %10s %10s %12s %12s', 'case', 'res', 'eN', 'smin/smax PM', 'smin/smax PP')};
for c = [0 1 5; 30 1 5; 60 1 5; 60 2 5; 90 1 5; 0 1 8; 30 1 8; 60 1 8; 60 2 8; 90 1 8; 90 2 8].'
    alpha = c(1);  N = c(2);  pdeg = c(3);
    m = cyl_model(fullfile(inpDir, sprintf('ThickCylinder_N%d_a%03d.inp', N, alpha)), pdeg);
    nq = ceil(pdeg^2 + pdeg/2);
    K = cyl_assemble(m, C, 'PM', nq);
    F = cyl_loads(m, prob, nq);
    [u, uEx, act] = cyl_solve(m, K, F, prob);
    fix = setdiff((1:m.num_DoF).', act);
    r = K(act, act)*uEx(act) + K(act, fix)*uEx(fix) - F(act);
    eN = max(abs(u - uEx))/max(abs(uEx));
    s1 = svd(full(K(act, act)));
    K2 = cyl_assemble(m, C, 'PP');
    s2 = svd(full(K2(act, act)));
    L{end+1} = sprintf('alpha=%2d N=%d p=%d nq=%2d      %10.2e %10.2e %12.2e %12.2e', ...
        alpha, N, pdeg, nq, norm(r)/norm(F(act)), eN, s1(end)/s1(1), s2(end)/s2(1)); %#ok<AGROW>
    fprintf('%s\n', L{end});
end
fid = fopen(outfile, 'w');  fprintf(fid, '%s\n', L{:});  fclose(fid);
end
