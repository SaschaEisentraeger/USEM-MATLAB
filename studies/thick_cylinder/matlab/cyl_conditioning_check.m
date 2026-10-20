function cyl_conditioning_check(outfile)
%CYL_CONDITIONING_CHECK Consistency vs. stability on selected meshes.
%   For the polynomial field (in both trial spaces of the unsymmetric
%   element for p >= 5) the stiffness is integrated with nq = p^2 + p/2
%   (rounded up) points, which is exact for PM (see run_quadrature_study).
%   Reported per case and trial space (PM: skew coordinates, PH: hybrid):
%     res   = ||K_aa u_ex,a + K_ab u_b - F_a|| / ||F_a||  (consistency residual
%             of the exact nodal values; zero in exact arithmetic for PM)
%     eN    = nodal error of the computed solution
%     smin  = smallest singular value of K_aa (scaled by the largest);
%             also for the conventional element (PP)
%   "rejected": the code base stops (ShapeFunctions2D_Metric).
here = fileparts(mfilename('fullpath'));
if nargin < 1, outfile = fullfile(here, '..', 'results', 'conditioning.txt'); end
cyl_setup_paths();
inpDir = fullfile(here, '..', 'mesh', 'inp');
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(1000, 0.25));
prob = cyl_problem('poly', 1000, 0.25);
L = {sprintf('%-26s | %-32s | %-32s | %10s', 'case', 'PM: res  eN  smin/smax', ...
             'PH: res  eN  smin/smax', 'PP smin/smax')};
for c = [0 1 5; 30 1 5; 60 1 5; 60 2 5; 90 1 5; 0 1 8; 30 1 8; 60 1 8; 60 2 8; 90 1 8; 90 2 8].'
    alpha = c(1);  N = c(2);  pdeg = c(3);
    m = cyl_model(fullfile(inpDir, sprintf('ThickCylinder_N%d_a%03d.inp', N, alpha)), pdeg);
    nq = ceil(pdeg^2 + pdeg/2);
    F = cyl_loads(m, prob, nq);
    cells = {};
    for form = {'PM', 'PH'}
        try
            K = cyl_assemble(m, C, form{1}, nq);
        catch err
            if ~startsWith(err.identifier, 'ShapeFunctions2D_Metric:'), rethrow(err); end
            cells{end+1} = sprintf('%-32s', 'rejected'); %#ok<AGROW>
            continue
        end
        [u, uEx, act] = cyl_solve(m, K, F, prob);
        fix = setdiff((1:m.num_DoF).', act);
        r = K(act, act)*uEx(act) + K(act, fix)*uEx(fix) - F(act);
        eN = max(abs(u - uEx))/max(abs(uEx));
        s = svd(full(K(act, act)));
        cells{end+1} = sprintf('%9.2e %9.2e %10.2e  ', norm(r)/norm(F(act)), eN, s(end)/s(1)); %#ok<AGROW>
    end
    K2 = cyl_assemble(m, C, 'PP');
    act = setdiff((1:m.num_DoF).', unique([2*m.set.xsym-1; 2*m.set.ysym; 2*m.set.xsym; 2*m.set.ysym-1]));
    s2 = svd(full(K2(act, act)));
    L{end+1} = sprintf('alpha=%2d N=%d p=%d nq=%2d   | %s | %s | %10.2e', alpha, N, pdeg, nq, ...
                       cells{1}, cells{2}, s2(end)/s2(1)); %#ok<AGROW>
    fprintf('%s\n', L{end});
end
fid = fopen(outfile, 'w');  fprintf(fid, '%s\n', L{:});  fclose(fid);
end
