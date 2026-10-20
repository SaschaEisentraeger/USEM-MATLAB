function T = unisolvence_check(cases, outfile)
%UNISOLVENCE_CHECK Factorization of the hybrid interpolation matrix.
%   For the hybrid trial space T_p = P_m(x) + span{L_a(xi)L_b(eta): a+b>m}
%   the interpolation matrix factorizes as det(V_T) = +-det(V_ref)*det(A),
%   where V_ref is the Legendre Vandermonde matrix of the reference GLL grid
%   and A (dim P_m x dim P_m) maps the Legendre coefficients of q in P_m(x)
%   onto the low-order (a+b <= m) Legendre coefficients of the Q_p(xi)
%   interpolant of q o x. V_T is singular iff A is.
%   cases: rows [alpha N p] of the thick-cylinder meshes. For every element
%   the script computes rcond(V_T) (as ShapeFunctions2D_Metric, method 2),
%   rcond(A), |det V_T|/|det V_ref det A| - 1 and rcond(V_Q) of Q_p(s,t)
%   (Lagrange basis, method 1), and reports the minimum over the elements.
here = fileparts(mfilename('fullpath'));
cyl = fullfile(here, '..', '..', 'studies', 'thick_cylinder', 'matlab');
addpath(cyl);  cyl_setup_paths();
T = zeros(size(cases, 1), 8);
for i = 1:size(cases, 1)
    a = cases(i, 1);  N = cases(i, 2);  pdeg = cases(i, 3);
    m = cyl_model(fullfile(cyl, '..', 'mesh', 'inp', sprintf('ThickCylinder_N%d_a%03d.inp', N, a)), pdeg);
    p = m.p;  n = (p.pxi+1)*(p.peta+1);  mm = min(p.pxi, p.peta);
    [xiN, etaN] = ReferenceNodes2D(p);
    [Lxi, ~] = PolynomialBasis1D(xiN, p.pxi, 'legendre');
    [Leta, ~] = PolynomialBasis1D(etaN, p.peta, 'legendre');
    rcT = inf;  rcA = inf;  rcQ = inf;  dmax = 0;  nbad = 0;
    for e = 1:m.num_Elem
        ids = m.Topology_Node(e, :).';
        x = m.xCoor(ids);  y = m.yCoor(ids);
        [sN, tN] = CorotationalSkewCoordinates(x, y, x, y, p);
        sN = (sN - (max(sN)+min(sN))/2)/((max(sN)-min(sN))/2);
        tN = (tN - (max(tN)+min(tN))/2)/((max(tN)-min(tN))/2);
        [Ls, ~] = PolynomialBasis1D(sN, p.pxi, 'legendre');
        [Lt, ~] = PolynomialBasis1D(tN, p.peta, 'legendre');
        VT = zeros(n);  Vref = zeros(n);  low = false(n, 1);  c = 0;
        for b = 0:p.peta
            for aa = 0:p.pxi
                c = c + 1;
                Vref(c, :) = Lxi(aa+1, :).*Leta(b+1, :);
                if aa + b <= mm
                    VT(c, :) = Ls(aa+1, :).*Lt(b+1, :);  low(c) = true;
                else
                    VT(c, :) = Vref(c, :);
                end
            end
        end
        % Legendre coefficients (in xi) of the interpolants of the P_m(x) rows
        Cq = (Vref.') \ (VT(low, :).');      % n x dim P_m
        A = Cq(low, :);
        rcT = min(rcT, rcond(VT));  rcA = min(rcA, rcond(A));
        dmax = max(dmax, abs(abs(det(VT))/abs(det(Vref)*det(A)) - 1));
        [~, ~, ~, VQ] = ShapeFunctions2D_Metric_noCheck(x, y, p);
        rcQ = min(rcQ, rcond(VQ));
        nbad = nbad + (rcond(VT) < eps);
    end
    T(i, :) = [a, N, pdeg, rcT, rcA, dmax, rcQ, nbad];
    fprintf('alpha=%2d N=%d p=%d: min rcond(V_T) %.2e  min rcond(A) %.2e  max|det ratio-1| %.1e  min rcond(V_Q) %.2e  elements with rcond(V_T)<eps: %d\n', ...
            a, N, pdeg, rcT, rcA, dmax, rcQ, nbad);
end
if nargin >= 2
    fid = fopen(outfile, 'w');
    fprintf(fid, 'alpha N p rcondVT rcondA detrel rcondVQ nsingular\n');
    fprintf(fid, '%d %d %d %.3e %.3e %.1e %.3e %d\n', T.');
    fclose(fid);
end
end

function [N, dNx, dNy, V] = ShapeFunctions2D_Metric_noCheck(x, y, p)
% Interpolation matrix of Q_p(s,t) (method 1, Lagrange basis) without the
% unisolvence check of the code base
[sN, tN] = CorotationalSkewCoordinates(x, y, x, y, p);
sN = (sN - (max(sN)+min(sN))/2)/((max(sN)-min(sN))/2);
tN = (tN - (max(tN)+min(tN))/2)/((max(tN)-min(tN))/2);
[Ps, ~] = PolynomialBasis1D(sN, p.pxi, 'lagrange', p.ndist);
[Pt, ~] = PolynomialBasis1D(tN, p.peta, 'lagrange', p.ndist);
n = (p.pxi+1)*(p.peta+1);  V = zeros(n);  c = 0;
for b = 0:p.peta
    for a = 0:p.pxi
        c = c + 1;  V(c, :) = Ps(a+1, :).*Pt(b+1, :);
    end
end
N = [];  dNx = [];  dNy = [];
end
