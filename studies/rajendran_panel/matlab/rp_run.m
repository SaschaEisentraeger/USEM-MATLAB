function rp_run(outdir)
%RP_RUN Nodal rotation patch test and convergence test of Rajendran (2010),
%   Secs. 5.5-5.7, with the spectral elements of the code base.
%
%   Variants (columns of the result files):
%     PP9, MM9, PM9 : p = 2 replicas of PP-QUAD9, MM-QUAD9 and PM-QUAD9*
%                     (9th node at the corner centroid, metric functions in
%                     Cartesian coordinates, 5 x 5 Gauss points)
%     PM9skew       : as PM9, metric functions in skew coordinates (Ke_2D)
%     PMusem        : p = 2 unsymmetric element of the code base as is
%                     (Ke_2D: skew coordinates, 4 x 4 points, 9th node from
%                     the 8-node map)
%     PPp, PMp      : p = 3..8, interior nodes from the bilinear corner map,
%                     skew coordinates, 2p x 2p points (Ke_2D)
%     PMpusem       : p = 3..8, code base as is (interior nodes from the
%                     8-node map)
%   Writes outdir/patch.csv (Meshes a-f, theta = 0..270 deg, three loads)
%   and outdir/convergence.csv (sequences A, B, C; shear load; N = 1, 2, 4).
here = fileparts(mfilename('fullpath'));
if nargin < 1, outdir = fullfile(here, '..', 'results'); end
rp_setup_paths();
inp = @(N, t) fullfile(here, '..', 'mesh', 'inp', sprintf('Panel_N%d_t%03d.inp', N, t));
E = 2000;  nu = 0.25;
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(E, nu));
loads = {'tension', 'moment', 'shear'};
for k = 1:3
    prob{k} = rp_exact(loads{k}, E, nu);   %#ok<AGROW>
    ref{k} = rp_reference(prob{k}, C);    %#ok<AGROW>
end
% variant: name, p, interior nodes, metric method, nq, form
V = {'PP9', 2, 'bilinear', 0, 5, 'PP'; 'MM9', 2, 'bilinear', 0, 5, 'MM';
     'PM9', 2, 'bilinear', 0, 5, 'PM'; 'PM9skew', 2, 'bilinear', 1, 5, 'PM';
     'PMusem', 2, 'serendipity', 1, 4, 'PM'};
for pd = 3:8
    V(end+1, :) = {sprintf('PP%d', pd), pd, 'bilinear', 1, 2*pd, 'PP'};      %#ok<AGROW>
    V(end+1, :) = {sprintf('PM%d', pd), pd, 'bilinear', 1, 2*pd, 'PM'};      %#ok<AGROW>
    V(end+1, :) = {sprintf('PM%dusem', pd), pd, 'serendipity', 1, 2*pd, 'PM'}; %#ok<AGROW>
end
fid = fopen(fullfile(outdir, 'patch.csv'), 'w');
fprintf(fid, 'variant,theta,load,uxA,uyA,energy,nodal,negDetJ\n');
fclose(fid);
for t = [0 90 135 175 225 270]
    rows = solveCases(inp(1, t), 1, t, 1:3, false);
    fid = fopen(fullfile(outdir, 'patch.csv'), 'a');  fprintf(fid, '%s', rows{:});  fclose(fid);
end
fid = fopen(fullfile(outdir, 'convergence.csv'), 'w');
fprintf(fid, 'variant,theta,N,load,uxA,uyA,energy,nodal,negDetJ\n');
fclose(fid);
for t = [0 90 135]
    for N = [1 2 4]
        rows = solveCases(inp(N, t), N, t, 3, true);
        fid = fopen(fullfile(outdir, 'convergence.csv'), 'a');  fprintf(fid, '%s', rows{:});  fclose(fid);
    end
end

    function rows = solveCases(file, N, t, kk, withN)
        % kk: load cases; withN: write the mesh parameter N (convergence file)
        rows = {};
        for v = 1:size(V, 1)
            try
                m = rp_model(file, V{v, 2}, V{v, 3});
                m.metricMethod = V{v, 4};
                K = rp_assemble(m, C, V{v, 5}, V{v, 6}, V{v, 4});
                neg = negJ(m);
            catch err
                % coincident nodes (pSpectralElement_NodesTopology) or
                % singular metric interpolation (ShapeFunctions2D_Metric):
                % the element cannot be built, the case is recorded as NaN
                if ~any(strcmp(err.identifier, {'pSpectralElement_NodesTopology:coincidentNodes', ...
                                                 'ShapeFunctions2D_Metric:singular'}))
                    rethrow(err);
                end
                for k = kk
                    if ~withN
                        rows{end+1} = sprintf('%s,%d,%s,NaN,NaN,NaN,NaN,NaN\n', V{v, 1}, t, loads{k}); %#ok<AGROW>
                    else
                        rows{end+1} = sprintf('%s,%d,%d,%s,NaN,NaN,NaN,NaN,NaN\n', V{v, 1}, t, N, loads{k}); %#ok<AGROW>
                    end
                end
                continue
            end
            for k = kk
                F = rp_loads(m, prob{k});
                u = rp_solve(m, K, F);
                uEx = zeros(m.num_DoF, 1);
                uEx(1:2:end) = prob{k}.ux(m.xCoor, m.yCoor);
                uEx(2:2:end) = prob{k}.uy(m.xCoor, m.yCoor);
                errForm = V{v, 6};
                er = cyl_errors(m, u, uEx, prob{k}, C, errForm, [], ref{k});
                if ~withN
                    rows{end+1} = sprintf('%s,%d,%s,%.10f,%.10f,%.6e,%.6e,%d\n', V{v, 1}, t, ...
                        loads{k}, u(2*m.A-1), u(2*m.A), er.energySigned, er.nodal, neg); %#ok<AGROW>
                else
                    rows{end+1} = sprintf('%s,%d,%d,%s,%.10f,%.10f,%.6e,%.6e,%d\n', V{v, 1}, t, N, ...
                        loads{k}, u(2*m.A-1), u(2*m.A), er.energySigned, er.nodal, neg); %#ok<AGROW>
                end
            end
        end
        fprintf('theta = %3d, N = %d done\n', t, N);
    end
end

function n = negJ(m)
% number of elements with det J <= 0 somewhere (41 x 41 grid)
g = linspace(-1, 1, 41);  [G1, G2] = meshgrid(g, g);
[~, dNxi, dNeta] = ShapeFunctions2D_SEM(G1(:), G2(:), m.p);
n = 0;
for e = 1:m.num_Elem
    ids = m.Topology_Node(e, :).';
    [~, ~, ~, ~, dJ] = cyl_jacobian(dNxi, dNeta, m.xCoor(ids), m.yCoor(ids));
    n = n + any(dJ <= 0);
end
end
