function run_study(alphas, Ns, ps, outdir)
%RUN_STUDY Distortion and MMS convergence study, thick cylinder.
%
%   run_study(alphas, Ns, ps, outdir) solves, for every distortion angle in
%   alphas (degree), every mesh N (N x N super-elements, 2N x 2N elements)
%   and every polynomial degree in ps,
%     - the thick cylinder under internal pressure (Lame solution),
%     - the manufactured polynomial solution (total degree 5),
%     - the manufactured smooth non-polynomial solution,
%   with the unsymmetric spectral element (PM, code base) and the
%   conventional spectral element (PP). The stiffness matrices use the
%   default quadrature of Ke_2D (2p x 2p Gauss points).
%   One CSV file per angle is written to outdir:  study_aXXX.csv
%
%   Defaults: alphas = [0 15 30 45 60 90], Ns = [1 2 4 8], ps = 2:8,
%   outdir = ../results/raw

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(alphas), alphas = [0 15 30 45 60 90]; end
if nargin < 2 || isempty(Ns), Ns = [1 2 4 8]; end
if nargin < 3 || isempty(ps), ps = 2:8; end
if nargin < 4 || isempty(outdir), outdir = fullfile(here, '..', 'results', 'raw'); end
cyl_setup_paths();
if ~exist(outdir, 'dir'), mkdir(outdir); end
inpDir = fullfile(here, '..', 'mesh', 'inp');

E = 1000;  nu = 0.25;
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(E, nu));
kinds = {'lame', 'poly', 'smooth'};
forms = {'PM', 'PP'};
for k = 1:numel(kinds)
    prob{k} = cyl_problem(kinds{k}, E, nu);          %#ok<AGROW>
    ref{k} = cyl_reference_norms(prob{k}, C);       %#ok<AGROW>
end

for alpha = alphas
    fname = fullfile(outdir, sprintf('study_a%03d.csv', alpha));
    % The ABAQUS readers close all files, so the rows are appended per case
    fid = fopen(fname, 'w');
    fprintf(fid, ['problem,form,alpha,N,p,ndof,energy,L2,nodal,energySigned,' ...
                  'minDetJ,negDetJ,area,lebesgueM,lebesgueGLL,condV,condKe,condK\n']);
    fclose(fid);
    for N = Ns
        inp = fullfile(inpDir, sprintf('ThickCylinder_N%d_a%03d.inp', N, alpha));
        for pdeg = ps
            t0 = tic;
            model = cyl_model(inp, pdeg);
            dg = cyl_diagnostics(model, C);
            rows = {};
            for f = 1:numel(forms)
                K = cyl_assemble(model, C, forms{f});
                for k = 1:numel(kinds)
                    F = cyl_loads(model, prob{k});
                    [u, uEx, act] = cyl_solve(model, K, F, prob{k});
                    if k == 1
                        % 1-norm condition estimate (symmetry BCs); costly,
                        % therefore only for the coarse meshes
                        if N <= 2, condK = condest(K(act, act)); else, condK = NaN; end
                    end
                    er = cyl_errors(model, u, uEx, prob{k}, C, forms{f}, [], ref{k});
                    rows{end+1} = sprintf('%s,%s,%d,%d,%d,%d,%.6e,%.6e,%.6e,%.6e,%.6e,%d,%.15e,%.6e,%.6e,%.6e,%.6e,%.6e\n', ...
                        kinds{k}, forms{f}, alpha, N, pdeg, model.num_DoF, ...
                        er.energy, er.L2, er.nodal, er.energySigned, ...
                        model.minDetJ, model.negDetJ, model.area, ...
                        dg.lebesgue, dg.lebesgueGLL, dg.condV, ...
                        dg.(['condKe' forms{f}]), condK);
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
