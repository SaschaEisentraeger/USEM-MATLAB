function R = study_metric_basis(varargin)
%STUDY_METRIC_BASIS Basis of the generalized Vandermonde matrix of the
%metric (trial) shape functions: monomials vs. Legendre, Chebyshev and
%Lagrange polynomials in (bounding-box scaled) skew coordinates.
%
%   R = study_metric_basis()                 full study (p = 1..16)
%   R = study_metric_basis('pmax', 12, 'pmaxMesh', 10, 'writeData', false)
%
%   All bases span the same space Q_p(s,t), so the metric shape functions
%   are identical in exact arithmetic. The study measures what differs in
%   floating-point arithmetic:
%
%   Element level (test elements of usemTestElement):
%     kappa   2-norm condition number of the Vandermonde matrix V
%     eKron   max |M_j(x_k) - delta_jk|                (Kronecker property)
%     eRep    max |sum_k f(x_k) M_k - f| over the monomials f of total
%             degree <= p in normalized Cartesian coordinates, at the
%             (2p)^2 Gauss points used by Ke_2D        (completeness)
%     eGrad   the same for h*grad(f), divided by p
%     LebM    Lebesgue constant of the metric interpolation (basis
%             independent); LebN that of the parametric interpolation
%     asym    ||K_e - K_e^T||_F/||K_e||_F on the parallelogram (affine
%             element: K_e is symmetric in exact arithmetic)
%     time    run time of one call with (2p)^2 query points, relative to
%             the monomial basis ('general' element)
%   Mesh level (3x2 mesh, angular + curved-edge distortion 0.3):
%     patch   max relative nodal error of the harmonic patch test of
%             degree p (exact in exact arithmetic)
%     smooth  max relative nodal error for u = (e^x cos y, -e^x sin y)
%
%   With 'writeData' the tables are written to docs/figures/data/ for the
%   plots in docs/chapters/metric-basis.tex.
%
%   Runs in MATLAB and in GNU Octave (studies/metric_basis/octave provides
%   a stand-in for MATLAB's decomposition).

ip = inputParser;
ip.addParameter('pmax', 16);
ip.addParameter('pmaxMesh', 12);
ip.addParameter('writeData', true);
ip.parse(varargin{:});
opt = ip.Results;

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(fullfile(root, 'code', 'subroutines')));
addpath(fullfile(root, 'tests', 'helpers'));
if exist('OCTAVE_VERSION', 'builtin')
    addpath(fullfile(here, 'octave'));
    warning('off', 'all');
end

bases  = {'monomial', 'legendre', 'chebyshev', 'lagrange'};
tags   = {'mon', 'leg', 'che', 'lag'};
shapes = {'parallelogram', 'trapezoid', 'curved', 'general'};
pvec   = 1:opt.pmax;
nB = numel(bases);

%% Element level
for s = 1:numel(shapes)
    T = struct('p', pvec(:));
    for b = 1:nB
        T.(['kappa_' tags{b}]) = zeros(numel(pvec), 1);
        T.(['eKron_' tags{b}]) = zeros(numel(pvec), 1);
        T.(['eRep_'  tags{b}]) = zeros(numel(pvec), 1);
        T.(['eGrad_' tags{b}]) = zeros(numel(pvec), 1);
        T.(['time_'  tags{b}]) = zeros(numel(pvec), 1);
    end
    T.LebM = zeros(numel(pvec), 1);
    T.LebN = zeros(numel(pvec), 1);
    T.asym_mon = nan(numel(pvec), 1);
    T.asym_lag = nan(numel(pvec), 1);

    for ii = 1:numel(pvec)
        p = pvec(ii);
        [x, y, P] = usemTestElement(shapes{s}, p);
        nN = numel(x);

        % Gauss points of Ke_2D, mapped to the element
        [g, ~] = lgwt(2*p, -1, 1);
        [Xi, Eta] = meshgrid(g, g);
        Xg = Local2GlobalCoordinatesQuad(P, [Xi(:), Eta(:)], [x, y]);
        xq = Xg(:,1); yq = Xg(:,2);

        % Normalized Cartesian coordinates for the test polynomials
        xc = mean(x); yc = mean(y);
        h = max(sqrt((x - xc).^2 + (y - yc).^2));
        X = (x - xc)/h; Y = (y - yc)/h;
        Xq = (xq - xc)/h; Yq = (yq - yc)/h;

        for b = 1:nB
            P.metricBasis = bases{b};
            [Nk, ~, ~, V] = ShapeFunctions2D_Metric(x, y, x, y, P);
            [N, dNx, dNy] = ShapeFunctions2D_Metric(xq, yq, x, y, P);
            eRep = 0; eGrad = 0;
            for a = 0:p
                for c = 0:p-a
                    f = X.^a.*Y.^c;
                    fq = Xq.^a.*Yq.^c;
                    fx = a*Xq.^max(a-1, 0).*Yq.^c;
                    fy = c*Xq.^a.*Yq.^max(c-1, 0);
                    eRep  = max(eRep, max(abs(N'*f - fq)));
                    eGrad = max(eGrad, max(abs([h*dNx'*f - fx; h*dNy'*f - fy]))/p);
                end
            end
            T.(['kappa_' tags{b}])(ii) = cond(V);
            T.(['eKron_' tags{b}])(ii) = max(max(abs(Nk - eye(nN))));
            T.(['eRep_'  tags{b}])(ii) = eRep;
            T.(['eGrad_' tags{b}])(ii) = eGrad;
            T.(['time_'  tags{b}])(ii) = runTime(@() ...
                ShapeFunctions2D_Metric(xq, yq, x, y, P));
        end

        % Lebesgue constants (dense grid in the element)
        [Xi, Eta] = meshgrid(linspace(-1, 1, 4*p + 9));
        Xd = Local2GlobalCoordinatesQuad(P, [Xi(:), Eta(:)], [x, y]);
        P.metricBasis = 'lagrange';
        M = ShapeFunctions2D_Metric(Xd(:,1), Xd(:,2), x, y, P);
        T.LebM(ii) = max(sum(abs(M), 1));
        T.LebN(ii) = max(sum(abs(ShapeFunctions2D_SEM(Xi(:), Eta(:), P)), 1));

        % Symmetry of K_e on the affine element
        if strcmp(shapes{s}, 'parallelogram')
            C = HookeMatrix2d(1, 0.3, 1);
            for b = [1 4]
                P.metricBasis = bases{b};
                K = Ke_2D(P, C, x, y);
                T.(['asym_' tags{b}])(ii) = norm(K - K', 'fro')/norm(K, 'fro');
            end
        end
    end
    for b = nB:-1:1
        T.(['time_' tags{b}]) = T.(['time_' tags{b}])./T.time_mon;
    end
    R.(shapes{s}) = T;
    printElementTable(shapes{s}, T, tags);
end

%% Mesh level
pm = 1:opt.pmaxMesh;
M = struct('p', pm(:));
for b = 1:nB
    M.(['patch_'  tags{b}]) = zeros(numel(pm), 1);
    M.(['smooth_' tags{b}]) = zeros(numel(pm), 1);
end
[nodes, topo] = usemTestMesh(3, 2, 3, 2, 2);
nodes = usemDistortMesh(nodes, topo, 'combined', 0.3);
C = HookeMatrix2d(1, 0.3, 1);
for ii = 1:numel(pm)
    model = usemBuildModel(nodes, topo, pm(ii));
    for b = 1:nB
        model.p.metricBasis = bases{b};
        K = assembleK(model, C);
        k = pm(ii);
        [~, ~, M.(['patch_' tags{b}])(ii)] = usemDirichletSolve(K, model, ...
            @(x, y) harmonicField(x - 1.3, y - 0.7, k));
        [~, ~, M.(['smooth_' tags{b}])(ii)] = usemDirichletSolve(K, model, ...
            @(x, y) deal(exp(x).*cos(y), -exp(x).*sin(y)));
    end
end
R.mesh = M;
fprintf('\nMesh level (3x2 mesh, combined distortion 0.3)\n');
fprintf('%3s', 'p');
for b = 1:nB, fprintf('  patch_%s', tags{b}); end
for b = 1:nB, fprintf(' smooth_%s', tags{b}); end
fprintf('\n');
for ii = 1:numel(pm)
    fprintf('%3d', pm(ii));
    for b = 1:nB, fprintf('  %9.2e', M.(['patch_' tags{b}])(ii)); end
    for b = 1:nB, fprintf('  %9.2e', M.(['smooth_' tags{b}])(ii)); end
    fprintf('\n');
end

%% Data files for the documentation
if opt.writeData
    out = fullfile(root, 'docs', 'figures', 'data');
    if ~exist(out, 'dir'), mkdir(out); end
    for s = 1:numel(shapes)
        writeTable(fullfile(out, ['metric-basis-' shapes{s} '.dat']), R.(shapes{s}));
    end
    writeTable(fullfile(out, 'metric-basis-mesh.dat'), R.mesh);
end

end

function K = assembleK(model, C)
% usemAssemble without the multiple-output evalc (not supported by Octave)
material_elset = {(1:model.num_Elem)'};
evalc(['K = AssemblyFEM(model.Topology_Node, model.p, model.xCoor, ' ...
    'model.yCoor, model.num_Elem, model.num_DoF, {C}, {1}, 1, 1, ' ...
    'material_elset, 0, 0, {0}, 0);']);
end

function t = runTime(fun)
% Median of several runs
n = 7;
tt = zeros(n, 1);
for k = 1:n
    t0 = tic;
    fun();
    tt(k) = toc(t0);
end
t = median(tt);
end

function [ux, uy] = harmonicField(x, y, k)
% u = (Re z^k, -Im z^k): divergence-free and harmonic (equilibrium)
z = complex(x, y).^k;
ux = real(z);
uy = -imag(z);
end

function printElementTable(name, T, tags)
fprintf('\nElement "%s"\n', name);
fprintf('%3s', 'p');
cols = {'kappa', 'eKron', 'eRep', 'eGrad'};
for c = 1:numel(cols)
    for b = 1:numel(tags)
        fprintf(' %9s', [cols{c} '_' tags{b}]);
    end
end
fprintf(' %9s %9s %9s\n', 'LebM', 'LebN', 'time_lag');
for i = 1:numel(T.p)
    fprintf('%3d', T.p(i));
    for c = 1:numel(cols)
        for b = 1:numel(tags)
            fprintf(' %9.2e', T.([cols{c} '_' tags{b}])(i));
        end
    end
    fprintf(' %9.2e %9.2f %9.2f\n', T.LebM(i), T.LebN(i), T.time_lag(i));
end
if any(~isnan(T.asym_mon))
    fprintf('asym(K_e) monomial: %s\n', mat2str(T.asym_mon', 2));
    fprintf('asym(K_e) lagrange: %s\n', mat2str(T.asym_lag', 2));
end
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
