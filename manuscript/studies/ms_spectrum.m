function R = ms_spectrum(pvec, meshes)
%MS_SPECTRUM Patch test, smooth field and generalized eigenvalues on the
%distorted 3 x 2 meshes (manuscript, Secs. 5.4-5.5).
%
%   R = ms_spectrum(pvec, meshes)
%
%   Domain [0,3] x [0,2], 3 x 2 elements from an 8-node mesh:
%     'curved15'  : curved-edge distortion 0.15
%     'curved20'  : curved-edge distortion 0.20
%     'combined30': angular 0.3 + curved 0.15
%   For each p and formulation (cla: PP, skew: PM, hyb: PH):
%     patch   max relative nodal error, harmonic field of degree p
%     smooth  max relative nodal error, u = (e^x cos y, -e^x sin y)
%     eig     max relative error of the real parts of the 8 eigenvalues of
%             smallest modulus (u = 0 on the boundary, consistent mass,
%             rho = 1); reference: PP on the undistorted mesh, p = 12
%     neg     number of eigenvalues with Re(lambda) < 0
%     imag    max |Im(lambda)|/|lambda| over the whole spectrum
%     asym    ||K - K^T||_F/||K||_F of the assembled stiffness
%   Stiffness and mass matrices are integrated with the calibrated Gauss
%   rules (ms_rule). Same quantities as studies/metric_basis/
%   study_hybrid_space.m, which uses the former 2p x 2p rule.
%   Writes manuscript/data/spectrum-<mesh>.dat.
if nargin < 1 || isempty(pvec), pvec = 1:10; end
if nargin < 2 || isempty(meshes), meshes = {'curved15', 'curved20', 'combined30'}; end
[~, dataDir] = ms_setup();
C = HookeMatrix2d(1, 0.3, 1);
nEig = 8;
[nodes0, topo] = usemTestMesh(3, 2, 3, 2, 2);
model = usemBuildModel(nodes0, topo, 12);
[K, M] = ms_assemble_KM(model, C, 'PP', 'rect3x2');
lamRef = real(lowestEigenvalues(K, M, model, nEig));

forms = {'PP', 'PM', 'PH'};
tags = {'cla', 'skew', 'hyb'};
for im = 1:numel(meshes)
    nodes = nodes0;
    switch meshes{im}
        case 'curved15',   nodes = usemDistortMesh(nodes, topo, 'curved', 0.15);
        case 'curved20',   nodes = usemDistortMesh(nodes, topo, 'curved', 0.20);
        case 'combined30', nodes = usemDistortMesh(nodes, topo, 'combined', 0.30);
    end
    T = struct('p', pvec(:));
    for t = 1:3
        for q = {'patch', 'smooth', 'eig', 'neg', 'imag', 'asym'}
            T.([q{1} '_' tags{t}]) = zeros(numel(pvec), 1);
        end
    end
    for i = 1:numel(pvec)
        p = pvec(i);
        model = usemBuildModel(nodes, topo, p);
        for t = 1:3
            [K, M] = ms_assemble_KM(model, C, forms{t}, meshes{im});
            [~, ~, T.(['patch_' tags{t}])(i)] = usemDirichletSolve(K, model, ...
                @(x, y) harmonicField(x - 1.3, y - 0.7, p));
            [~, ~, T.(['smooth_' tags{t}])(i)] = usemDirichletSolve(K, model, ...
                @(x, y) deal(exp(x).*cos(y), -exp(x).*sin(y)));
            [lam, lamAll] = lowestEigenvalues(K, M, model, nEig);
            nl = numel(lam);
            T.(['eig_'  tags{t}])(i) = max(abs(real(lam) - lamRef(1:nl))./lamRef(1:nl));
            T.(['neg_'  tags{t}])(i) = nnz(real(lamAll) < 0);
            T.(['imag_' tags{t}])(i) = max(abs(imag(lamAll))./abs(lamAll));
            T.(['asym_' tags{t}])(i) = norm(K - K', 'fro')/norm(K, 'fro');
        end
        fprintf(['%-10s p = %2d | patch %8.1e %8.1e %8.1e | smooth %8.1e %8.1e %8.1e' ...
                 ' | eig %8.1e %8.1e %8.1e | neg %3d %3d %3d | imag %8.1e %8.1e %8.1e\n'], ...
                meshes{im}, p, T.patch_cla(i), T.patch_skew(i), T.patch_hyb(i), ...
                T.smooth_cla(i), T.smooth_skew(i), T.smooth_hyb(i), ...
                T.eig_cla(i), T.eig_skew(i), T.eig_hyb(i), ...
                T.neg_cla(i), T.neg_skew(i), T.neg_hyb(i), ...
                T.imag_cla(i), T.imag_skew(i), T.imag_hyb(i));
    end
    R.(meshes{im}) = T;
    ms_write_table(fullfile(dataDir, ['spectrum-' meshes{im} '.dat']), T);
end
end

function [lam, lamAll] = lowestEigenvalues(K, M, model, nEig)
% u = 0 on the boundary; the nEig eigenvalues of smallest modulus and the
% whole spectrum (as in study_hybrid_space)
x = model.xCoor; y = model.yCoor;
tol = 1e-10*max(max(x) - min(x), max(y) - min(y));
bnd = abs(x - min(x)) < tol | abs(x - max(x)) < tol | ...
      abs(y - min(y)) < tol | abs(y - max(y)) < tol;
act = setdiff((1:model.num_DoF)', sort([2*find(bnd)-1; 2*find(bnd)]));
lamAll = eig(K(act, act), M(act, act));
[~, k] = sort(abs(lamAll));
lam = lamAll(k(1:min(nEig, numel(lamAll))));
end

function [ux, uy] = harmonicField(x, y, k)
z = complex(x, y).^k;
ux = real(z);
uy = -imag(z);
end
