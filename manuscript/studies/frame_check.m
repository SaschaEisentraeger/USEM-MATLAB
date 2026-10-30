function frame_check(outfile)
%FRAME_CHECK Coordinates of the complete polynomials of the hybrid space.
%   The hybrid shape functions do not depend on the affine coordinates in
%   which P_p(x,y) is written, but the Vandermonde matrix V_T does. For
%     1 skew coordinates (code base),
%     2 principal axes of the nodes,
%     3 affine least-squares fit of the reference coordinates to the nodes,
%   all with the bounding-box scaling, the script computes on the test
%   elements and on all thick-cylinder meshes (p = 2..8)
%     kappa   largest cond(V_T) over the elements,
%     rc      smallest rcond(V_T) over the elements,
%     nbad    number of elements with rcond(V_T) < eps,
%     eK      largest max|V_T\V_T - I| over the elements.
%   Writes outfile (default manuscript/data/frames.txt), Appendix B.
[~, dataDir] = ms_setup();
if nargin < 1, outfile = fullfile(dataDir, 'frames.txt'); end
here = fileparts(mfilename('fullpath'));
cyl = fullfile(here, '..', '..', 'studies', 'thick_cylinder', 'matlab');
addpath(cyl);  cyl_setup_paths();
% The ABAQUS readers close all files, so the rows are written at the end
rows = {sprintf(['case alpha N p kappaSkew kappaPrincipal kappaFit rcSkew rcPrincipal rcFit ' ...
                 'nbadSkew nbadPrincipal nbadFit eKSkew eKPrincipal eKFit\n'])};
for shape = {'parallelogram', 'trapezoid', 'curved', 'general'}
    for pdeg = [4 8 12 16]
        [x, y, p] = usemTestElement(shape{1}, pdeg);
        S = accumulate(initStats(), x, y, p);
        rows{end+1} = formatRow(shape{1}, 0, 0, pdeg, S); %#ok<AGROW>
    end
end
for alpha = [0 15 30 45 60 90]
    for N = [1 2 4 8]
        for pdeg = 2:8
            m = cyl_model(fullfile(cyl, '..', 'mesh', 'inp', sprintf('ThickCylinder_N%d_a%03d.inp', N, alpha)), pdeg);
            S = initStats();
            for e = 1:m.num_Elem
                ids = m.Topology_Node(e, :).';
                S = accumulate(S, m.xCoor(ids), m.yCoor(ids), m.p);
            end
            rows{end+1} = formatRow('cylinder', alpha, N, pdeg, S); %#ok<AGROW>
        end
    end
end
fid = fopen(outfile, 'w');
fprintf(fid, '%s', rows{:});
fclose(fid);
end

function S = initStats()
S = struct('kappa', zeros(1, 3), 'rc', inf(1, 3), 'nbad', zeros(1, 3), 'eK', zeros(1, 3));
end

function S = accumulate(S, x, y, p)
for f = 1:3
    V = hybridV(x, y, p, f);
    r = rcond(V);
    S.kappa(f) = max(S.kappa(f), cond(V));
    S.rc(f) = min(S.rc(f), r);
    S.nbad(f) = S.nbad(f) + (r < eps);
    S.eK(f) = max(S.eK(f), max(max(abs(V\V - eye(size(V))))));
end
end

function row = formatRow(name, alpha, N, pdeg, S)
row = sprintf('%s %d %d %d %.3e %.3e %.3e %.3e %.3e %.3e %d %d %d %.3e %.3e %.3e\n', ...
        name, alpha, N, pdeg, S.kappa, S.rc, S.nbad, S.eK);
end

function V = hybridV(x, y, p, frame)
[xiN, etaN] = ReferenceNodes2D(p);
X = [x - mean(x), y - mean(y)];
switch frame
    case 1
        [s, t] = CorotationalSkewCoordinates(x, y, x, y, p);
    case 2
        [~, ~, W] = svd(X, 0);
        Z = X*W;  s = Z(:, 1);  t = Z(:, 2);
    case 3
        Z = X*(X\[xiN(:), etaN(:)]);  s = Z(:, 1);  t = Z(:, 2);
end
s = (s - (max(s)+min(s))/2)/((max(s)-min(s))/2);
t = (t - (max(t)+min(t))/2)/((max(t)-min(t))/2);
Ls = PolynomialBasis1D(s, p.pxi, 'legendre');
Lt = PolynomialBasis1D(t, p.peta, 'legendre');
Lxi = PolynomialBasis1D(xiN, p.pxi, 'legendre');
Leta = PolynomialBasis1D(etaN, p.peta, 'legendre');
V = zeros((p.pxi+1)*(p.peta+1));  c = 0;
for b = 0:p.peta
    for a = 0:p.pxi
        c = c + 1;
        if a + b <= min(p.pxi, p.peta)
            V(c, :) = Ls(a+1, :).*Lt(b+1, :);
        else
            V(c, :) = Lxi(a+1, :).*Leta(b+1, :);
        end
    end
end
end
