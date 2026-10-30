function run_study(alphas, Ns, ps, outdir, nqfun)
%RUN_STUDY Distortion and MMS convergence study, thick cylinder.
%
%   run_study(alphas, Ns, ps, outdir) solves, for every distortion angle in
%   alphas (degree), every mesh N (N x N super-elements, 2N x 2N elements)
%   and every polynomial degree in ps,
%     - the thick cylinder under internal pressure (Lame solution),
%     - the manufactured polynomial solution (total degree 5),
%     - the manufactured smooth non-polynomial solution,
%   with the unsymmetric spectral element of the code base with its two
%   trial spaces (PM: skew coordinates, ShapeFunctions2D_Metric method 1,
%   the default; PH: hybrid space, method 2) and the conventional spectral
%   element (PP). The stiffness matrices use the default quadrature of
%   Ke_2D (2p x 2p Gauss points). Cases the code base rejects (singular
%   metric interpolation, det J = 0 for the hybrid space) are written as
%   NaN. The columns lebesgue, condV, condKe and condK refer to the trial
%   space of the row (PP: lebesgue = GLL value, condV = NaN).
%   One CSV file per angle is written to outdir:  study_aXXX.csv
%
%   nqfun (optional): handle nq = nqfun(alpha, N, form, p) with the number of
%   Gauss points per direction of the stiffness matrices (e.g. the
%   calibrated rules of manuscript/studies/ms_rule); the loads then use
%   max(3p+2, nq) and the error integrals max(3p+2, nq+4) points. Without
%   nqfun, the stiffness uses 2p and loads and errors 3p+2 points.
%
%   Defaults: alphas = [0 15 30 45 60 90], Ns = [1 2 4 8], ps = 2:8,
%   outdir = ../results/raw

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(alphas), alphas = [0 15 30 45 60 90]; end
if nargin < 2 || isempty(Ns), Ns = [1 2 4 8]; end
if nargin < 3 || isempty(ps), ps = 2:8; end
if nargin < 4 || isempty(outdir), outdir = fullfile(here, '..', 'results', 'raw'); end
if nargin < 5, nqfun = []; end
cyl_setup_paths();
if ~exist(outdir, 'dir'), mkdir(outdir); end
inpDir = fullfile(here, '..', 'mesh', 'inp');

E = 1000;  nu = 0.25;
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(E, nu));
kinds = {'lame', 'poly', 'smooth'};
forms = {'PM', 'PH', 'PP'};
for k = 1:numel(kinds)
    prob{k} = cyl_problem(kinds{k}, E, nu);          %#ok<AGROW>
    ref{k} = cyl_reference_norms(prob{k}, C);       %#ok<AGROW>
end

for alpha = alphas
    fname = fullfile(outdir, sprintf('study_a%03d.csv', alpha));
    fid = fopen(fname, 'w');
    fprintf(fid, ['problem,form,alpha,N,p,ndof,energy,L2,nodal,energySigned,' ...
                  'minDetJ,negDetJ,area,lebesgue,lebesgueGLL,condV\n']);
    fclose(fid);
    for N = Ns
        inp = fullfile(inpDir, sprintf('ThickCylinder_N%d_a%03d.inp', N, alpha));
        for pdeg = ps
            t0 = tic;
            model = cyl_model(inp, pdeg);
            rows = {};
            dg = cyl_diagnostics(model);
            leb = struct('PM', dg.lebesgue, 'PH', dg.lebesguePH, 'PP', dg.lebesgueGLL);
            cV = struct('PM', dg.condV, 'PH', dg.condVPH, 'PP', NaN);
            for f = 1:numel(forms)
                try
                    nqK = cyl_rule(nqfun, alpha, N, forms{f}, pdeg);
                    K = cyl_assemble(model, C, forms{f}, nqK);
                catch err
                    % trial space rejected by the code base on this mesh
                    if ~startsWith(err.identifier, 'ShapeFunctions2D_Metric:'), rethrow(err); end
                    for k = 1:numel(kinds)
                        rows{end+1} = sprintf('%s,%s,%d,%d,%d,%d,NaN,NaN,NaN,NaN,%.6e,%d,%.15e,NaN,NaN,NaN\n', ...
                            kinds{k}, forms{f}, alpha, N, pdeg, model.num_DoF, ...
                            model.minDetJ, model.negDetJ, model.area); %#ok<AGROW>
                    end
                    continue
                end
                for k = 1:numel(kinds)
                    F = cyl_loads(model, prob{k}, max(3*pdeg + 2, nqK));
                    [u, uEx] = cyl_solve(model, K, F, prob{k});
                    er = cyl_errors(model, u, uEx, prob{k}, C, forms{f}, max(3*pdeg + 2, nqK + 4), ref{k});
                    rows{end+1} = sprintf('%s,%s,%d,%d,%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%d,%.15e,%.6e,%.6e,%.6e\n', ...
                        kinds{k}, forms{f}, alpha, N, pdeg, model.num_DoF, ...
                        er.energy, er.L2, er.nodal, er.energySigned, ...
                        model.minDetJ, model.negDetJ, model.area, ...
                        leb.(forms{f}), dg.lebesgueGLL, cV.(forms{f}));
                end
            end
            fid = fopen(fname, 'a');
            fprintf(fid, '%s', rows{:});
            fclose(fid);
            fprintf('alpha = %3d, N = %d, p = %d: %6d dof, %.1f s\n', ...
                    alpha, N, pdeg, model.num_DoF, toc(t0));
        end
    end
end
end
