function R = study_hybrid_space(varargin)
%STUDY_HYBRID_SPACE Mesh-level comparison of the trial spaces of the
%unsymmetric spectral element: Q_p(s,t) (method 1), hybrid space (method 2)
%and the classical (symmetric, isoparametric) spectral element.
%
%   R = study_hybrid_space()                                  all meshes
%   R = study_hybrid_space('mesh', 'curved30', 'pvec', 1:10)  one mesh
%
%   Domain [0,3] x [0,2], 3 x 2 elements from an 8-node mesh with
%   curved-edge distortion of the interior mid-side nodes:
%     'curved15'  : curved, amount 0.15   (min det J = 0.16)
%     'curved20'  : curved, amount 0.20   (min det J = 0.09)
%     'combined30': angular 0.3 + curved 0.15 (min det J = 0.014)
%   For each p:
%     patch   max relative nodal error, harmonic field of degree p
%             (exact for methods 1 and 2)
%     smooth  max relative nodal error, u = (e^x cos y, -e^x sin y)
%     eig     max relative error of the 8 eigenvalues of smallest modulus
%             (u = 0 on the boundary, consistent mass, rho = 1); reference:
%             classical element on the undistorted mesh, p = pref
%     neg     number of eigenvalues lambda = omega^2 with Re(lambda) < 0
%             (whole spectrum; exponentially growing modes in dynamics)
%     imag    max |Im(lambda)|/|lambda| over the whole spectrum
%     asym    ||K - K^T||_F/||K||_F of the assembled stiffness
%
%   Writes docs/figures/data/metric-hybrid-<mesh>.dat ('writeData').
%   Runs in MATLAB and in GNU Octave.

ip = inputParser;
ip.addParameter('mesh', {'curved15', 'curved20', 'combined30'});
ip.addParameter('pvec', 1:10);
ip.addParameter('pref', 12);
ip.addParameter('nEig', 8);
ip.addParameter('writeData', true);
ip.parse(varargin{:});
opt = ip.Results;
if ischar(opt.mesh), opt.mesh = {opt.mesh}; end

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(fullfile(root, 'code', 'subroutines')));
addpath(fullfile(root, 'tests', 'helpers'));
if exist('OCTAVE_VERSION', 'builtin')
    addpath(fullfile(here, 'octave'));
    warning('off', 'all');
end

nu = 0.3;
C = HookeMatrix2d(1, nu, 1);
[nodes0, topo] = usemTestMesh(3, 2, 3, 2, 2);

% Reference eigenvalues: classical element, undistorted mesh
model = usemBuildModel(nodes0, topo, opt.pref);
[K, M] = assembleClassical(model, C);
lamRef = real(lowestEigenvalues(K, M, model, opt.nEig));

tags = {'cla', 'skew', 'hyb'};
for im = 1:numel(opt.mesh)
    nodes = distort(nodes0, topo, opt.mesh{im});
    pvec = opt.pvec(:);
    T = struct('p', pvec);
    for t = 1:3
        for q = {'patch', 'smooth', 'eig', 'neg', 'imag', 'asym'}
            T.([q{1} '_' tags{t}]) = zeros(numel(pvec), 1);
        end
    end
    for i = 1:numel(pvec)
        p = pvec(i);
        model = usemBuildModel(nodes, topo, p);
        for t = 1:3
            if t == 1
                [K, M] = assembleClassical(model, C);
            else
                model.p.metricMethod = t - 1;
                [K, M] = assembleUnsymmetric(model, C);
            end
            [~, ~, T.(['patch_' tags{t}])(i)] = usemDirichletSolve(K, model, ...
                @(x, y) harmonicField(x - 1.3, y - 0.7, p));
            [~, ~, T.(['smooth_' tags{t}])(i)] = usemDirichletSolve(K, model, ...
                @(x, y) deal(exp(x).*cos(y), -exp(x).*sin(y)));
            [lam, lamAll] = lowestEigenvalues(K, M, model, opt.nEig);
            nl = numel(lam);
            T.(['eig_'  tags{t}])(i) = max(abs(real(lam) - lamRef(1:nl))./lamRef(1:nl));
            T.(['neg_'  tags{t}])(i) = nnz(real(lamAll) < 0);
            T.(['imag_' tags{t}])(i) = max(abs(imag(lamAll))./abs(lamAll));
            T.(['asym_' tags{t}])(i) = norm(full(K - K'), 'fro')/norm(full(K), 'fro');
        end
        fprintf(['%-10s p = %2d | patch %8.1e %8.1e %8.1e | smooth %8.1e %8.1e %8.1e' ...
                 ' | eig %8.1e %8.1e %8.1e | neg %3d %3d %3d | imag %8.1e %8.1e %8.1e\n'], opt.mesh{im}, p, ...
                T.patch_cla(i), T.patch_skew(i), T.patch_hyb(i), ...
                T.smooth_cla(i), T.smooth_skew(i), T.smooth_hyb(i), ...
                T.eig_cla(i), T.eig_skew(i), T.eig_hyb(i), ...
                T.neg_cla(i), T.neg_skew(i), T.neg_hyb(i), ...
                T.imag_cla(i), T.imag_skew(i), T.imag_hyb(i));
    end
    R.(opt.mesh{im}) = T;
    if opt.writeData
        writeTable(fullfile(root, 'docs', 'figures', 'data', ...
            ['metric-hybrid-' opt.mesh{im} '.dat']), T);
    end
end

end

function nodes = distort(nodes, topo, name)
switch name
    case 'curved15',   nodes = usemDistortMesh(nodes, topo, 'curved', 0.15);
    case 'curved20',   nodes = usemDistortMesh(nodes, topo, 'curved', 0.20);
    case 'combined30', nodes = usemDistortMesh(nodes, topo, 'combined', 0.30);
    otherwise, error('study_hybrid_space:mesh', 'Unknown mesh "%s".', name);
end
end

function [K, M] = assembleUnsymmetric(model, C)
% AssemblyFEM with consistent mass (evalc with one output for Octave)
material_elset = {(1:model.num_Elem)'};
evalc(['[K, M] = AssemblyFEM(model.Topology_Node, model.p, model.xCoor, ' ...
    'model.yCoor, model.num_Elem, model.num_DoF, {C}, {1}, 1, 1, ' ...
    'material_elset, 0, 0, {0}, 1);']);
end

function [K, M] = assembleClassical(model, C)
% Bubnov-Galerkin spectral element (parametric test and trial functions)
n = model.num_DoF;
K = zeros(n);
M = zeros(n);
T = Topology_DoF_2D(model.Topology_Node);
P = model.p;
[g, w] = lgwt(2*P.pxi, -1, 1);
[Xi, Eta] = meshgrid(g, g);
[Wx, We] = meshgrid(w, w);
W = Wx(:).*We(:);
N = ShapeFunctions2D_SEM(Xi(:), Eta(:), P);
for e = 1:model.num_Elem
    nd = model.Topology_Node(e,:);
    x = model.xCoor(nd); y = model.yCoor(nd);
    Ke = usemKeSymmetric(P, C, x, y);
    Me = zeros(2*numel(nd));
    for k = 1:numel(W)
        J = Jacobi_Matrix(P, [Xi(k) Eta(k)], x, y);
        Nk = kron(N(:,k)', eye(2));
        Me = Me + (Nk'*Nk)*det(J)*W(k);
    end
    K(T(e,:), T(e,:)) = K(T(e,:), T(e,:)) + Ke;
    M(T(e,:), T(e,:)) = M(T(e,:), T(e,:)) + Me;
end
end

function [lam, lamAll] = lowestEigenvalues(K, M, model, nEig)
% u = 0 on the boundary of the rectangle; the nEig eigenvalues of smallest
% modulus (spurious eigenvalues with negative real part have large moduli)
% and the whole spectrum
x = model.xCoor; y = model.yCoor;
tol = 1e-10*max(max(x) - min(x), max(y) - min(y));
bnd = abs(x - min(x)) < tol | abs(x - max(x)) < tol | ...
      abs(y - min(y)) < tol | abs(y - max(y)) < tol;
act = setdiff((1:model.num_DoF)', sort([2*find(bnd)-1; 2*find(bnd)]));
lamAll = eig(full(K(act, act)), full(M(act, act)));
[~, k] = sort(abs(lamAll));
lam = lamAll(k(1:min(nEig, numel(lamAll))));
end

function [ux, uy] = harmonicField(x, y, k)
z = complex(x, y).^k;
ux = real(z);
uy = -imag(z);
end

function writeTable(file, T)
names = fieldnames(T);
fid = fopen(file, 'w');
fprintf(fid, '%s ', names{:});
fprintf(fid, '\n');
for i = 1:numel(T.p)
    for k = 1:numel(names)
        fprintf(fid, '%.6e ', T.(names{k})(i));
    end
    fprintf(fid, '\n');
end
fclose(fid);
end
