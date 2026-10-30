function ms_cylinder(alphas)
%MS_CYLINDER Thick-cylinder study of the manuscript (Sec. 5.8) with the
%calibrated Gauss rules.
%
%   ms_cylinder(alphas) runs studies/thick_cylinder/matlab/run_study and
%   run_interpolation_study for the distortion angles alphas (default: all)
%   with the stiffness rules of ms_quadrature_calibration (families
%   cyl_aAAA_NN); loads use max(3p+2, nG) and the error integrals
%   max(3p+2, nG+4) points. The results are written to
%   studies/thick_cylinder/results/raw_calibrated (study_aXXX.csv,
%   interp_aXXX.csv); the original results of the study report (2p rule)
%   in results/raw are not touched.
if nargin < 1 || isempty(alphas), alphas = [0 15 30 45 60 90]; end
root = ms_setup();
cyl = fullfile(root, 'studies', 'thick_cylinder', 'matlab');
addpath(cyl);
cyl_setup_paths();
outdir = fullfile(root, 'studies', 'thick_cylinder', 'results', 'raw_calibrated');
if ~exist(outdir, 'dir'), mkdir(outdir); end
nqfun = @(alpha, N, form, p) ms_rule(sprintf('cyl_a%03d_N%d', alpha, N), form, p);
for a = alphas
    run_study(a, [1 2 4 8], 2:8, outdir, nqfun);
    run_interpolation_study(a, [1 2 4 8], 2:8, outdir, nqfun);
end
end
