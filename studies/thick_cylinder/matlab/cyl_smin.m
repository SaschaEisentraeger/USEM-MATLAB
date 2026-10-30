function S = cyl_smin(cases, outfile)
%CYL_SMIN Smallest singular values of the global stiffness matrices.
%   cases: rows [alpha N p]. For each case and the three elements (PM, PP,
%   PH) the smallest singular value of K_aa (symmetry BCs of the Lame
%   problem) is computed (svds, shift 0) together with the largest one; the
%   ratios PM/PP and PH/PP of the smallest singular values measure the loss
%   of stability of the unsymmetric discretizations relative to the
%   conventional one on the same mesh. The nodal location of the largest
%   entry of the PM singular vector is reported (r, theta).
here = fileparts(mfilename('fullpath'));
cyl_setup_paths();
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(1000, 0.25));
S = zeros(size(cases, 1), 10);
for i = 1:size(cases, 1)
    a = cases(i, 1);  N = cases(i, 2);  pdeg = cases(i, 3);
    m = cyl_model(fullfile(here, '..', 'mesh', 'inp', sprintf('ThickCylinder_N%d_a%03d.inp', N, a)), pdeg);
    fix = unique([2*m.set.xsym - 1; 2*m.set.ysym]);
    act = setdiff((1:m.num_DoF).', fix);
    out = zeros(1, 6);  loc = [NaN NaN];
    for f = 1:3
        forms = {'PM', 'PP', 'PH'};
        K = cyl_assemble(m, C, forms{f});
        A = K(act, act);
        [~, smin, V] = svds(A, 1, 0);
        smax = svds(A, 1);
        out(2*f-1:2*f) = [smin, smax];
        if f == 1
            v = zeros(m.num_DoF, 1);  v(act) = V(:, 1);
            [~, j] = max(hypot(v(1:2:end), v(2:2:end)));
            loc = [hypot(m.xCoor(j), m.yCoor(j)), atan2(m.yCoor(j), m.xCoor(j))*180/pi];
        end
    end
    S(i, :) = [a, N, pdeg, out(1), out(3), out(1)/out(3), loc, out(5), out(5)/out(3)];
    fprintf('alpha=%2d N=%d p=%d: smin PM %.3e  PP %.3e  ratio %.3e | smax PM %.3e PP %.3e | PM mode max at r=%.3f theta=%.1f | PH smin %.3e ratio %.3e smax %.3e\n', ...
            a, N, pdeg, out(1), out(3), out(1)/out(3), out(2), out(4), loc, out(5), out(5)/out(3), out(6));
end
if nargin >= 2
    fid = fopen(outfile, 'w');
    fprintf(fid, 'alpha N p sminPM sminPP ratio r_mode theta_mode sminPH ratioPH\n');
    fprintf(fid, '%d %d %d %.4e %.4e %.4e %.3f %.2f %.4e %.4e\n', S.');
    fclose(fid);
end
end
