function T = cyl_element_analysis(N, alpha, pdeg, kind, outfile)
%CYL_ELEMENT_ANALYSIS Element-wise Lebesgue constants and error shares.
%   For one mesh and degree, computes per element: the Lebesgue constant
%   of the metric interpolation, cond of the element Vandermonde matrix,
%   sigma_max/sigma_(2n-3) of the PM element stiffness, and the share of
%   the squared energy error (PM solution) and of the squared interpolation
%   error contributed by the element. Writes a table (sorted by the
%   Lebesgue constant) to outfile if given.
here = fileparts(mfilename('fullpath'));
cyl_setup_paths();
inp = fullfile(here, '..', 'mesh', 'inp', sprintf('ThickCylinder_N%d_a%03d.inp', N, alpha));
C = PlaneStressConstitutiveMatrices(HookeMatrix3d(1000, 0.25));
prob = cyl_problem(kind, 1000, 0.25);
ref = cyl_reference_norms(prob, C);
m = cyl_model(inp, pdeg);
K = cyl_assemble(m, C, 'PM');
F = cyl_loads(m, prob);
[u, uEx] = cyl_solve(m, K, F, prob);
nE = m.num_Elem;
T = zeros(nE, 7);
g = linspace(-1, 1, 31);  [G1, G2] = meshgrid(g, g);
Ng = ShapeFunctions2D_SEM(G1(:), G2(:), m.p);
full = m;
for e = 1:nE
    ids = m.Topology_Node(e, :).';
    x = m.xCoor(ids);  y = m.yCoor(ids);
    M = ShapeFunctions2D_Metric(Ng.'*x, Ng.'*y, x, y, m.p);
    s = svd(cyl_Ke_PM(m.p, C, x, y, 2*pdeg));
    sub = full;  sub.Topology_Node = m.Topology_Node(e, :);  sub.num_Elem = 1;
    eF = cyl_errors(sub, u, uEx, prob, C, 'PM', [], ref);
    eI = cyl_errors(sub, uEx, uEx, prob, C, 'PM', [], ref);
    c = mean([x y]);
    T(e, :) = [e, max(sum(abs(M), 1)), s(1)/s(end-3), eF.energySigned, eI.energySigned, hypot(c(1), c(2)), atan2(c(2), c(1))*180/pi];
end
T(:, 4) = T(:, 4)/sum(T(:, 4));  T(:, 5) = T(:, 5)/sum(T(:, 5));
[~, o] = sort(T(:, 2), 'descend');  T = T(o, :);
if nargin >= 5
    fid = fopen(outfile, 'w');
    fprintf(fid, '%% N=%d alpha=%d p=%d problem=%s\n', N, alpha, pdeg, kind);
    fprintf(fid, 'elem lebesgueM condKe shareErrFE shareErrInterp r theta\n');
    fprintf(fid, '%d %.4e %.4e %.4e %.4e %.4f %.2f\n', T.');
    fclose(fid);
end
end
