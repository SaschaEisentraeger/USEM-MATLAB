function frame_cylinder_check(outfile)
%FRAME_CYLINDER_CHECK Thick-cylinder study of PH with the complete polynomials
%written in the principal axes of the nodes (Appendix B).
%   Solves the three problems of run_study for PH with the calibrated rules on
%   all meshes without folded elements and p = 2..8 and lists the errors next
%   to those with the skew coordinates (results/raw_calibrated, run ms_cylinder
%   first). Writes outfile (default manuscript/data/frames-cylinder.txt).
[root, dataDir] = ms_setup();
if nargin < 1, outfile = fullfile(dataDir, 'frames-cylinder.txt'); end
cyl = fullfile(root, 'studies', 'thick_cylinder');
addpath(fullfile(cyl, 'matlab'));
cyl_setup_paths();
axesDir = fullfile(fileparts(mfilename('fullpath')), 'principal_axes');

E = 1000;  nu = 0.25;
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(E, nu));
kinds = {'lame', 'poly', 'smooth'};
for k = 1:numel(kinds)
    prob{k} = cyl_problem(kinds{k}, E, nu);     %#ok<AGROW>
    ref{k} = cyl_reference_norms(prob{k}, C);   %#ok<AGROW>
end
rows = {sprintf('problem alpha N p energySkew energyPrincipal nodalSkew nodalPrincipal\n')};
for alpha = [0 15 30 45]
    T = readtable(fullfile(cyl, 'results', 'raw_calibrated', sprintf('study_a%03d.csv', alpha)));
    T = T(strcmp(T.form, 'PH'), :);
    for N = [1 2 4 8]
        for pdeg = 2:8
            model = cyl_model(fullfile(cyl, 'mesh', 'inp', sprintf('ThickCylinder_N%d_a%03d.inp', N, alpha)), pdeg);
            if model.negDetJ > 0, continue, end
            nq = ms_rule(sprintf('cyl_a%03d_N%d', alpha, N), 'PH', pdeg);
            % ms_rule calls ms_setup, which puts the code base first on the path
            addpath(axesDir, '-begin');
            assert(startsWith(which('CorotationalSkewCoordinates'), axesDir));
            K = cyl_assemble(model, C, 'PH', nq);
            for k = 1:numel(kinds)
                F = cyl_loads(model, prob{k}, max(3*pdeg + 2, nq));
                [u, uEx] = cyl_solve(model, K, F, prob{k});
                er = cyl_errors(model, u, uEx, prob{k}, C, 'PH', max(3*pdeg + 2, nq + 4), ref{k});
                i = strcmp(T.problem, kinds{k}) & T.N == N & T.p == pdeg;
                rows{end+1} = sprintf('%s %d %d %d %.6e %.6e %.6e %.6e\n', kinds{k}, alpha, N, pdeg, ...
                                      T.energy(i), er.energy, T.nodal(i), er.nodal); %#ok<AGROW>
            end
            fprintf('alpha = %d, N = %d, p = %d done\n', alpha, N, pdeg);
        end
    end
end
fid = fopen(outfile, 'w');
fprintf(fid, '%s', rows{:});
fclose(fid);
end
