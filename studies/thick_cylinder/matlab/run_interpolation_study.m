function run_interpolation_study(alphas, Ns, ps, outdir, nqfun)
%RUN_INTERPOLATION_STUDY Interpolation errors in the trial spaces.
%   For the same meshes as run_study, the exact nodal values are inserted
%   into the trial functions of the unsymmetric element (metric functions,
%   element by element; PM: skew coordinates, PH: hybrid space) and of the
%   conventional element (parametric functions, PP), and the relative errors are evaluated with cyl_errors. The
%   ratio (discretization error)/(interpolation error) separates the
%   approximation property of the trial space from the stability of the
%   discrete problem. Writes outdir/interp_aXXX.csv. The optional nqfun is
%   that of run_study, so that both studies use the same error integrals.
here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(alphas), alphas = [0 15 30 45 60 90]; end
if nargin < 2 || isempty(Ns), Ns = [1 2 4 8]; end
if nargin < 3 || isempty(ps), ps = 2:8; end
if nargin < 4 || isempty(outdir), outdir = fullfile(here, '..', 'results', 'raw'); end
if nargin < 5, nqfun = []; end
cyl_setup_paths();
inpDir = fullfile(here, '..', 'mesh', 'inp');
E = 1000;  nu = 0.25;
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(E, nu));
kinds = {'lame', 'poly', 'smooth'};
for k = 1:3
    prob{k} = cyl_problem(kinds{k}, E, nu);      %#ok<AGROW>
    ref{k} = cyl_reference_norms(prob{k}, C);   %#ok<AGROW>
end
for alpha = alphas
    fname = fullfile(outdir, sprintf('interp_a%03d.csv', alpha));
    fid = fopen(fname, 'w');
    fprintf(fid, 'problem,form,alpha,N,p,energy,L2\n');
    fclose(fid);
    for N = Ns
        for pdeg = ps
            model = cyl_model(fullfile(inpDir, sprintf('ThickCylinder_N%d_a%03d.inp', N, alpha)), pdeg);
            rows = {};
            for form = {'PM', 'PH', 'PP'}
                for k = 1:3
                    uEx = zeros(model.num_DoF, 1);
                    uEx(1:2:end) = prob{k}.ux(model.xCoor, model.yCoor);
                    uEx(2:2:end) = prob{k}.uy(model.xCoor, model.yCoor);
                    try
                        nqK = cyl_rule(nqfun, alpha, N, form{1}, pdeg);
                        er = cyl_errors(model, uEx, uEx, prob{k}, C, form{1}, max(3*pdeg + 2, nqK + 4), ref{k});
                    catch err
                        % trial space rejected by the code base
                        if ~startsWith(err.identifier, 'ShapeFunctions2D_Metric:'), rethrow(err); end
                        er = struct('energy', NaN, 'L2', NaN);
                    end
                    rows{end+1} = sprintf('%s,%s,%d,%d,%d,%.6e,%.6e\n', kinds{k}, ...
                        form{1}, alpha, N, pdeg, er.energy, er.L2); %#ok<AGROW>
                end
            end
            fid = fopen(fname, 'a');
            fprintf(fid, '%s', rows{:});
            fclose(fid);
        end
        fprintf('alpha = %d, N = %d done\n', alpha, N);
    end
end
end
