function R = ms_timing(pvec, shape)
%MS_TIMING Run times (median of 7 runs after one warm-up run, otherwise idle
%machine) of the metric shape functions at the (2p)^2 Gauss points and of
%the element stiffness matrix (manuscript, Sec. 6). Writes
%manuscript/data/timing.dat. Run with matlab -singleCompThread (one core).
if nargin < 1, pvec = 2:2:12; end
if nargin < 2, shape = 'general'; end
[~, dataDir] = ms_setup();
C = HookeMatrix2d(1, 0.3);
T = struct('p', pvec(:));
for q = {'sfMon', 'sfLag', 'sfHyb', 'kePP', 'keMon', 'keLag', 'keHyb'}
    T.(q{1}) = nan(numel(pvec), 1);
end
for i = 1:numel(pvec)
    p = pvec(i);
    [x, y, P] = usemTestElement(shape, p);
    [g, ~] = lgwt(2*p, -1, 1);
    [Xi, Eta] = meshgrid(g);
    IP = [Xi(:), Eta(:)];
    Xg = Local2GlobalCoordinatesQuad(P, IP, [x, y]);
    Pm = P; Pm.metricBasis = 'monomial'; Pm.metricMethod = 1;
    Pl = P; Pl.metricBasis = 'lagrange'; Pl.metricMethod = 1;
    Ph = P; Ph.metricMethod = 2;
    T.sfMon(i) = runTime(@() ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, Pm, 1, IP));
    T.sfLag(i) = runTime(@() ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, Pl, 1, IP));
    T.sfHyb(i) = runTime(@() ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, Ph, 2, IP));
    T.kePP(i)  = runTime(@() usemKeSymmetric(P, C, x, y));
    T.keMon(i) = runTime(@() Ke_2D(Pm, C, x, y));
    T.keLag(i) = runTime(@() Ke_2D(Pl, C, x, y));
    T.keHyb(i) = runTime(@() Ke_2D(Ph, C, x, y));
    fprintf('p=%2d  SF: mon %7.4f lag %7.4f hyb %7.4f | Ke: PP %7.3f mon %7.3f lag %7.3f hyb %7.3f s\n', ...
            p, T.sfMon(i), T.sfLag(i), T.sfHyb(i), T.kePP(i), T.keMon(i), T.keLag(i), T.keHyb(i));
end
ms_write_table(fullfile(dataDir, 'timing.dat'), T);
R = T;
end

function t = runTime(fun)
fun();   % untimed warm-up run (JIT compilation, file caching)
tt = zeros(7, 1);
for k = 1:7
    t0 = tic;
    fun();
    tt(k) = toc(t0);
end
t = median(tt);
end
