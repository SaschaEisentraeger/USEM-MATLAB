function R = ms_transient(p, meshName, nPeriods, stepsPerPeriod)
%MS_TRANSIENT Free vibration with the trapezoidal rule (Newmark,
%beta = 1/4, gamma = 1/2) on a curved 3 x 2 mesh (manuscript, Sec. 5.5).
%
%   R = ms_transient(p, meshName, nPeriods, stepsPerPeriod)
%
%   [0,3] x [0,2], u = 0 on the boundary, consistent mass, rho = 1,
%   E = 1, nu = 0.3 (plane stress). Initial displacement: nodal
%   interpolant of u0 = (sin(pi x/3) sin(pi y/2), sin(2 pi x/3) sin(pi y/2))
%   plus a nodal perturbation of relative size 1e-10 (round-off level), zero
%   initial velocity. dt = T1/stepsPerPeriod with the fundamental period T1
%   of the classical element. For the classical element (PP) and the
%   unsymmetric element with Q_p(s,t) (PM) and with T_p (PH), the history
%   of max|u| / max|u0| is written to manuscript/data/transient-<mesh>-p<p>.dat
%   together with the spectral data (number of eigenvalues with negative
%   real part, max |Im|/|lambda|).

if nargin < 1 || isempty(p), p = 6; end
if nargin < 2 || isempty(meshName), meshName = 'curved20'; end
if nargin < 3 || isempty(nPeriods), nPeriods = 50; end
if nargin < 4 || isempty(stepsPerPeriod), stepsPerPeriod = 20; end
[~, dataDir] = ms_setup();

C = HookeMatrix2d(1, 0.3, 1);
[nodes, topo] = usemTestMesh(3, 2, 3, 2, 2);
switch meshName
    case 'curved15',   nodes = usemDistortMesh(nodes, topo, 'curved', 0.15);
    case 'curved20',   nodes = usemDistortMesh(nodes, topo, 'curved', 0.20);
    case 'combined30', nodes = usemDistortMesh(nodes, topo, 'combined', 0.30);
end
model = usemBuildModel(nodes, topo, p);
x = model.xCoor; y = model.yCoor;
bnd = abs(x) < 1e-10 | abs(x - 3) < 1e-10 | abs(y) < 1e-10 | abs(y - 2) < 1e-10;
act = setdiff((1:model.num_DoF)', sort([2*find(bnd) - 1; 2*find(bnd)]));

u0 = zeros(model.num_DoF, 1);
u0(1:2:end) = sin(pi*x/3).*sin(pi*y/2);
u0(2:2:end) = sin(2*pi*x/3).*sin(pi*y/2);
u0 = u0(act);
pert = 1e-10*cos(17*(1:numel(u0))');       % deterministic, round-off level
u0 = u0 + pert*max(abs(u0));

forms = {'PP', 'PM', 'PH'};
for f = 1:numel(forms)
    [K, M] = ms_assemble_KM(model, C, forms{f});
    K = K(act, act); M = M(act, act);
    lam = eig(full(K), full(M));
    S.(forms{f}) = struct('neg', nnz(real(lam) < 0), ...
        'imag', max(abs(imag(lam))./abs(lam)), 'lamMin', min(real(lam)));
    KM.(forms{f}) = {K, M};
    if f == 1
        T1 = 2*pi/sqrt(min(real(lam)));
    end
end
dt = T1/stepsPerPeriod;
nSteps = nPeriods*stepsPerPeriod;
R = struct('t', (0:nSteps)'*dt/T1);
for f = 1:numel(forms)
    K = KM.(forms{f}){1}; M = KM.(forms{f}){2};
    a0 = 4/dt^2; a1 = 4/dt;
    Keff = K + a0*M;
    [L, U, Pp] = lu(Keff);
    u = u0; v = zeros(size(u)); a = M\(-K*u);
    h = zeros(nSteps + 1, 1); h(1) = 1;
    for n = 1:nSteps
        rhs = M*(a0*u + a1*v + a);
        un = U\(L\(Pp*rhs));
        an = a0*(un - u) - a1*v - a;
        v = v + dt/2*(a + an);
        u = un; a = an;
        h(n + 1) = max(abs(u))/max(abs(u0));
        if ~isfinite(h(n + 1)) || h(n + 1) > 1e30
            h(n + 1:end) = NaN;
            break
        end
    end
    R.(['amp_' forms{f}]) = h;
    fprintf('%s p=%d %s: neg %3d  max|Im|/|lam| %8.2e  max amp %9.3e  final amp %9.3e\n', ...
            meshName, p, forms{f}, S.(forms{f}).neg, S.(forms{f}).imag, ...
            max(h), h(find(isfinite(h), 1, 'last')));
end
ms_write_table(fullfile(dataDir, sprintf('transient-%s-p%d.dat', meshName, p)), R);
R.spectral = S;
end
