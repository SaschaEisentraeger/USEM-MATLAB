function run_quadrature_study(outfile)
%RUN_QUADRATURE_STUDY Sensitivity of the results to the stiffness quadrature.
%   Solves the three problems with nq = 2p (default of Ke_2D), 3p, 4p and 6p
%   Gauss points per direction for selected meshes and degrees, for the
%   unsymmetric element (PM: skew coordinates, PH: hybrid trial space) and
%   the conventional (PP) element; cases rejected by the code base are NaN.
%   The PM integrand is a polynomial in (xi,eta) of degree (2p-1)(q+1) per direction for a
%   geometry of degree q, i.e. nq >= p^2 + p/2 points are needed for exact
%   integration on the elements with exact arcs (q = p; nq = 6p suffices
%   for p <= 5) and nq >= 3p - 1 on the interior elements (q = 2). The PH
%   integrand is not polynomial (its high-order part is defined in the
%   reference coordinates, its gradient contains J^-1).
here = fileparts(mfilename('fullpath'));
if nargin < 1, outfile = fullfile(here, '..', 'results', 'raw', 'quadrature.csv'); end
cyl_setup_paths();
inpDir = fullfile(here, '..', 'mesh', 'inp');
E = 1000;  nu = 0.25;
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(E, nu));
kinds = {'lame', 'poly', 'smooth'};
for k = 1:3
    prob{k} = cyl_problem(kinds{k}, E, nu);      %#ok<AGROW>
    ref{k} = cyl_reference_norms(prob{k}, C);   %#ok<AGROW>
end
fid = fopen(outfile, 'w');
fprintf(fid, 'problem,form,alpha,N,p,nq,energy,nodal\n');
fclose(fid);
for alpha = [0 30 60 90]
    for N = [1 2]
        for pdeg = [5 8]
            model = cyl_model(fullfile(inpDir, sprintf('ThickCylinder_N%d_a%03d.inp', N, alpha)), pdeg);
            rows = {};
            for form = {'PM', 'PH', 'PP'}
                for nq = [2 3 4 6]*pdeg
                    try
                        K = cyl_assemble(model, C, form{1}, nq);
                    catch err
                        % trial space rejected by the code base
                        if ~startsWith(err.identifier, 'ShapeFunctions2D_Metric:'), rethrow(err); end
                        for k = 1:3
                            rows{end+1} = sprintf('%s,%s,%d,%d,%d,%d,NaN,NaN\n', kinds{k}, ...
                                form{1}, alpha, N, pdeg, nq); %#ok<AGROW>
                        end
                        continue
                    end
                    for k = 1:3
                        F = cyl_loads(model, prob{k});
                        [u, uEx] = cyl_solve(model, K, F, prob{k});
                        er = cyl_errors(model, u, uEx, prob{k}, C, form{1}, [], ref{k});
                        rows{end+1} = sprintf('%s,%s,%d,%d,%d,%d,%.6e,%.6e\n', kinds{k}, ...
                            form{1}, alpha, N, pdeg, nq, er.energy, er.nodal); %#ok<AGROW>
                    end
                end
            end
            fid = fopen(outfile, 'a');
            fprintf(fid, '%s', rows{:});
            fclose(fid);
            fprintf('alpha = %d, N = %d, p = %d done\n', alpha, N, pdeg);
        end
    end
end
end
