function R = ms_convergence(delta, pvec, Nvec, tag, rule)
%MS_CONVERGENCE h-convergence on meshes with persistent curved-edge
%distortion (manuscript, Sec. 5.3).
%
%   R = ms_convergence(delta, pvec, Nvec, tag, rule)
%
%   Unit square, N x N 8-node elements lifted to degree p (GLL). Every
%   interior mid-side node is moved normal to its edge by delta*h (h = 1/N),
%   so the element shapes are the same on every mesh (self-similar, not
%   asymptotically affine). Exact field u = (Re f, -Im f), f(z) = exp(2z)
%   (equilibrium without body force for any isotropic material), prescribed
%   at all boundary nodes. Formulations:
%     PP  classical (Bubnov-Galerkin) spectral element
%     PM  unsymmetric element, trial space Q_p(s,t) (method 1, Lagrange basis)
%     PH  unsymmetric element, hybrid trial space T_p (method 2)
%   Errors (relative): energy norm and L2 norm with the strains/displacements
%   of the trial functions of each formulation, maximum nodal error.
%   The stiffness matrices are integrated with the calibrated Gauss rule of
%   the mesh family conv<100*delta> (ms_rule), the errors with nG + 4
%   points per direction (as in the companion paper). With rule = '2p', the
%   stiffness matrices use the former 2p rule instead.
%   Writes manuscript/data/conv-<tag>.dat (one row per (p, N)).

if nargin < 1 || isempty(delta), delta = 0.2; end
if nargin < 2 || isempty(pvec), pvec = 2:6; end
if nargin < 3 || isempty(Nvec), Nvec = [1 2 4 8 16]; end
if nargin < 4 || isempty(tag), tag = sprintf('d%03d', round(100*delta)); end
if nargin < 5, rule = 'calibrated'; end
[~, dataDir] = ms_setup();

C = HookeMatrix2d(1, 0.3);
forms = {'PP', 'PM', 'PH'};
rows = [];
for p = pvec
    for N = Nvec
        [nodes, topo] = usemTestMesh(N, N, 1, 1, 2);
        if N > 1
            nodes = usemDistortMesh(nodes, topo, 'curved', delta/N);
        end
        model = usemBuildModel(nodes, topo, p);
        row = [p, N, model.num_DoF, minDetJ(model)];
        for f = 1:numel(forms)
            t0 = tic;
            if strcmp(rule, '2p')
                nG = 2*p;
            else
                nG = ms_rule(sprintf('conv%03d', round(100*delta)), forms{f}, p);
            end
            K = assemble(model, C, forms{f}, nG);
            u = dirichletSolve(K, model);
            [eE, eL, en] = errors(model, u, C, forms{f}, nG + 4);
            row = [row, eE, eL, en]; %#ok<AGROW>
            fprintf('%s p=%d N=%2d %s: energy %9.3e  L2 %9.3e  nodal %9.3e  (%.1f s)\n', ...
                    tag, p, N, forms{f}, eE, eL, en, toc(t0));
        end
        rows = [rows; row]; %#ok<AGROW>
    end
end
names = {'p', 'N', 'ndof', 'minDetJ'};
for f = 1:numel(forms)
    names = [names, {['energy_' forms{f}], ['L2_' forms{f}], ['nodal_' forms{f}]}]; %#ok<AGROW>
end
R = struct();
for k = 1:numel(names)
    R.(names{k}) = rows(:, k);
end
ms_write_table(fullfile(dataDir, ['conv-' tag '.dat']), R);
end

function [ux, uy, exx, eyy, gxy] = exactField(x, y)
e = exp(2*x);
ux = e.*cos(2*y);
uy = -e.*sin(2*y);
exx = 2*e.*cos(2*y);
eyy = -2*e.*cos(2*y);
gxy = -4*e.*sin(2*y);
end

function K = assemble(model, C, form, nq)
P = model.p;
nE = model.num_Elem;
T = Topology_DoF_2D(model.Topology_Node);
nd = size(T, 2);
I = zeros(nd^2, nE); J = I; V = I;
for e = 1:nE
    nodes = model.Topology_Node(e, :);
    x = model.xCoor(nodes); y = model.yCoor(nodes);
    Ke = ms_Ke(P, C, x, y, nq, form);
    [c, r] = meshgrid(T(e, :), T(e, :));
    I(:, e) = r(:); J(:, e) = c(:); V(:, e) = Ke(:);
end
K = sparse(I(:), J(:), V(:), model.num_DoF, model.num_DoF);
end

function u = dirichletSolve(K, model)
x = model.xCoor; y = model.yCoor;
[ux, uy] = exactField(x, y);
uEx = zeros(model.num_DoF, 1);
uEx(1:2:end) = ux; uEx(2:2:end) = uy;
tol = 1e-10;
bnd = abs(x) < tol | abs(x - 1) < tol | abs(y) < tol | abs(y - 1) < tol;
bD = sort([2*find(bnd) - 1; 2*find(bnd)]);
aD = setdiff((1:model.num_DoF)', bD);
u = uEx;
u(aD) = -K(aD, aD) \ (K(aD, bD)*uEx(bD));
end

function [eE, eL, en] = errors(model, u, C, form, nq)
P = model.p;
[g, w] = lgwt(nq, -1, 1);
[Xi, Eta] = meshgrid(g); [W1, W2] = meshgrid(w);
IP = [Xi(:), Eta(:)]; W = W1(:).*W2(:);
[Np, dNxi, dNeta] = ShapeFunctions2D_SEM(IP(:,1), IP(:,2), P);
num = [0 0]; den = [0 0];
x = model.xCoor; y = model.yCoor;
[uxN, uyN] = exactField(x, y);
en = max(abs([u(1:2:end) - uxN; u(2:2:end) - uyN]))/max(abs([uxN; uyN]));
for e = 1:model.num_Elem
    nodes = model.Topology_Node(e, :);
    xe = x(nodes); ye = y(nodes);
    ue = u(2*nodes - 1); ve = u(2*nodes);
    xg = Np'*xe; yg = Np'*ye;
    % Jacobian J = [x_xi y_xi; x_eta y_eta] at the quadrature points
    a = dNxi'*xe; b = dNxi'*ye; c = dNeta'*xe; d = dNeta'*ye;
    dJ = a.*d - b.*c;
    if strcmp(form, 'PP')
        N = Np;
        dNx = ( d'.*dNxi - b'.*dNeta)./dJ';
        dNy = (-c'.*dNxi + a'.*dNeta)./dJ';
    else
        P.metricMethod = 1 + strcmp(form, 'PH');
        [N, dNx, dNy] = ShapeFunctions2D_Metric(xg, yg, xe, ye, P, [], IP);
    end
    [ux, uy, exx, eyy, gxy] = exactField(xg, yg);
    dux = ux - N'*ue; duy = uy - N'*ve;
    de = [exx - dNx'*ue, eyy - dNy'*ve, gxy - dNy'*ue - dNx'*ve];
    ee = [exx, eyy, gxy];
    wd = W.*dJ;
    num(1) = num(1) + sum(sum((de*C).*de, 2).*wd);
    den(1) = den(1) + sum(sum((ee*C).*ee, 2).*wd);
    num(2) = num(2) + sum((dux.^2 + duy.^2).*wd);
    den(2) = den(2) + sum((ux.^2 + uy.^2).*wd);
end
eE = sqrt(abs(num(1))/den(1));
eL = sqrt(abs(num(2))/den(2));
end

function m = minDetJ(model)
P = model.p;
[g, ~] = lgwt(2*P.pxi, -1, 1);
[Xi, Eta] = meshgrid(g);
[~, dNxi, dNeta] = ShapeFunctions2D_SEM(Xi(:), Eta(:), P);
m = inf;
for e = 1:model.num_Elem
    nodes = model.Topology_Node(e, :);
    xe = model.xCoor(nodes); ye = model.yCoor(nodes);
    dJ = (dNxi'*xe).*(dNeta'*ye) - (dNxi'*ye).*(dNeta'*xe);
    m = min(m, min(dJ)*model.num_Elem);   % scaled by 1/h^2
end
end
