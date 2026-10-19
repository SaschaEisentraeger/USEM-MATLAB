function cyl_residual_check(cases, kind, outfile)
%CYL_RESIDUAL_CHECK Consistency residual of the exact nodal values.
%   For each case [alpha N p] and both elements: r = K_aa u_I,a + K_ab u_b - F_a
%   with the exact nodal values u_I (default quadrature). Since
%   u_h - u_I = K_aa^{-1} (-r), the nodal error is controlled by
%   ||r|| and the smallest singular value of K_aa. Reported: ||r||/||F_a||
%   and the nodal error e_n.
here = fileparts(mfilename('fullpath'));
cyl_setup_paths();
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(1000, 0.25));
prob = cyl_problem(kind, 1000, 0.25);
L = {sprintf('%-18s %12s %12s %12s %12s', ['case (' kind ')'], 'res PM', 'res PP', 'e_n PM', 'e_n PP')};
for i = 1:size(cases, 1)
    a = cases(i, 1);  N = cases(i, 2);  pdeg = cases(i, 3);
    m = cyl_model(fullfile(here, '..', 'mesh', 'inp', sprintf('ThickCylinder_N%d_a%03d.inp', N, a)), pdeg);
    F = cyl_loads(m, prob);
    res = [0 0];  en = [0 0];  forms = {'PM', 'PP'};
    for f = 1:2
        K = cyl_assemble(m, C, forms{f});
        [u, uEx, act] = cyl_solve(m, K, F, prob);
        fix = setdiff((1:m.num_DoF).', act);
        r = K(act, act)*uEx(act) + K(act, fix)*uEx(fix) - F(act);
        res(f) = norm(r)/norm(F(act));
        en(f) = max(hypot(u(1:2:end) - uEx(1:2:end), u(2:2:end) - uEx(2:2:end))) / ...
                max(hypot(uEx(1:2:end), uEx(2:2:end)));
    end
    L{end+1} = sprintf('alpha=%2d N=%d p=%d %12.3e %12.3e %12.3e %12.3e', a, N, pdeg, res, en); %#ok<AGROW>
    fprintf('%s\n', L{end});
end
fid = fopen(outfile, 'w');  fprintf(fid, '%s\n', L{:});  fclose(fid);
end
