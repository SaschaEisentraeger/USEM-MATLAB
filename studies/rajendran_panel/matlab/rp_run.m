function rp_run(outdir)
%RP_RUN Nodal rotation patch test of Rajendran (2010), Secs. 5.5-5.7, with
%   the spectral elements of the code base.
%
%   Variants (first column of the result file), p = 3..8, interior nodes
%   from the bilinear corner map, 2p x 2p Gauss points:
%     PPp : conventional spectral element
%     PMp : unsymmetric element, skew coordinates
%     PHp : unsymmetric element, hybrid trial space
%   Writes outdir/patch.csv (Meshes a-f, theta = 0..270 deg, three loads).
here = fileparts(mfilename('fullpath'));
rp_setup_paths();
E = 2000;  nu = 0.25;
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(E, nu));
loads = {'tension', 'moment', 'shear'};
for k = 1:3
    prob{k} = rp_exact(loads{k}, E, nu);   %#ok<AGROW>
    ref{k} = rp_reference(prob{k}, C);    %#ok<AGROW>
end
forms = {'PP', 'PM', 'PH'};
fid = fopen(fullfile(outdir, 'patch.csv'), 'w');
fprintf(fid, 'variant,theta,load,uxA,uyA,energy,nodal,negDetJ\n');
for t = [0 90 135 175 225 270]
    file = fullfile(here, '..', 'mesh', 'inp', sprintf('Panel_N1_t%03d.inp', t));
    for pd = 3:8
        m = rp_model(file, pd);
        neg = negJ(m);
        for f = 1:3
            K = cyl_assemble(m, C, forms{f}, 2*pd);
            for k = 1:3
                F = rp_loads(m, prob{k});
                u = rp_solve(m, K, F);
                uEx = zeros(m.num_DoF, 1);
                uEx(1:2:end) = prob{k}.ux(m.xCoor, m.yCoor);
                uEx(2:2:end) = prob{k}.uy(m.xCoor, m.yCoor);
                er = cyl_errors(m, u, uEx, prob{k}, C, forms{f}, [], ref{k});
                fprintf(fid, '%s%d,%d,%s,%.10f,%.10f,%.6e,%.6e,%d\n', forms{f}, pd, t, ...
                    loads{k}, u(2*m.A-1), u(2*m.A), er.energySigned, er.nodal, neg);
            end
        end
    end
    fprintf('theta = %3d done\n', t);
end
fclose(fid);
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
