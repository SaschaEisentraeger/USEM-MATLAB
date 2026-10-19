function [Xnew, stats] = distort_mesh(X, T, type, varargin)
%DISTORT_MESH Randomly distorts a structured 2D FE mesh (Q4/Q8) without inversion.
%
%   [Xnew, stats] = DISTORT_MESH(X, T, type, Name,Value,...)
%
%   Inputs
%   ------
%   X     : (nNode x 2) array of node coordinates [x y].
%   T     : (nElem x nen) integer connectivity (1-based), Q4: nen=4, Q8: nen=8 (serendipity).
%   type  : 'Q4' or 'Q8'.
%
%   Name-Value options
%   ------------------
%   'Level'       : Distortion level in [0,1]. Default 0.35. (0 = none, 1 = aggressive)
%   'Seed'        : RNG seed (default [] = don't touch RNG).
%   'MinJac'      : Strictly positive lower bound for det(J) at Gauss points after distortion.
%                   Default 1e-6 * median positive Jacobian before.
%   'MaxBacktrack': Max halving steps if invalid elements are detected (default 12).
%   'BoundaryTol' : Tolerance to detect outer boundary via min/max of x/y (default 1e-12).
%   'Verbose'     : true/false (default true).
%
%   Outputs
%   -------
%   Xnew  : (nNode x 2) distorted coordinates.
%   stats : struct with fields:
%           .scale        : final global scale applied to initial random field
%           .nBacktracks  : number of halving steps taken
%           .minJacBefore : min detJ before distortion
%           .minJacAfter  : min detJ after distortion
%           .movedMask    : logical mask of moved (interior) nodes
%
%   Notes
%   -----
%   * Only interior nodes (strictly inside the rectangular bounding box) are moved.
%   * For robustness, displacement magnitudes are locally bounded by the characteristic
%     edge length around each node and globally backtracked until all detJ > MinJac.
%   * Q8 is treated as 8-node serendipity (no center node). Node order expected:
%       [ (4)----(7)----(3)
%          |             |
%         (8)           (6)
%          |             |
%         (1)----(5)----(2) ]
%     but the code only relies on standard shape functions; reorder if needed.
%
%   Example
%   -------
%   % X: n-by-2, T4: m-by-4, from a structured grid generator
%   [Xnew, stats] = distort_mesh(X, T4, 'Q4', 'Level', 0.4, 'Seed', 1);
%
%   See also: delaunayTriangulation
%
%   Sascha-safe, Aug 2025.

% Parse
p = inputParser; p.FunctionName = 'distort_mesh';
p.addRequired('X', @(v)isnumeric(v)&&size(v,2)==2);
p.addRequired('T', @(v)isnumeric(v)&&ismatrix(v));
p.addRequired('type', @(s)ischar(s)||isstring(s));
p.addParameter('Level', 0.35, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('Seed', [], @(v)isempty(v)||(isscalar(v)&&isfinite(v)));
p.addParameter('MinJac', [], @(v)isempty(v)||(isscalar(v)&&v>0));
%-----Changes-----%
% p.addParameter('MaxBacktrack', 12, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('MaxBacktrack', 100, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
%-----Changes-----%
p.addParameter('BoundaryTol', 1e-12, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('Verbose', true, @(v)islogical(v)&&isscalar(v));
p.parse(X,T,type,varargin{:});
opt = p.Results; type = upper(string(type));

nNode = size(X,1);
Xnew = X; %#ok<NASGU> % set later

if ~any(type==["Q4","Q8"]) %#ok<*ISMAT>
    error('type must be ''Q4'' or ''Q8''');
end
nen = size(T,2);
if (type=="Q4" && nen~=4) || (type=="Q8" && nen~=8)
    error('Connectivity size (columns) does not match element type.');
end

% Build nodal adjacency to estimate characteristic length per node
% Precompute all edge pairs across elements
nElem = size(T,1);
% Each element contributes nen*(nen-1)/2 undirected edges
nEdgesMax = nElem * nen*(nen-1)/2;
I = zeros(nEdgesMax*2,1); % doubled for symmetric (i,j) and (j,i)
J = zeros(nEdgesMax*2,1);
k = 0;
for e=1:nElem
    nodes = T(e,:);
    for a=1:nen-1
        for b=a+1:nen
            k = k + 1;
            i = nodes(a); j = nodes(b);
            I(2*k-1:2*k) = [i; j];
            J(2*k-1:2*k) = [j; i];
        end
    end
end
I = I(1:2*k);
J = J(1:2*k);
V = ones(size(I));
adj = sparse(I,J,V,nNode,nNode);

[ii,jj] = find(adj);
if isempty(ii)
    error('Adjacency empty. Check connectivity.');
end
% characteristic length as median distance to adjacent nodes
Lij = sqrt(sum((X(ii,:)-X(jj,:)).^2,2));
Ls = accumarray(ii, Lij, [nNode 1], @(v) median(v));
Ls(Ls==0) = median(Ls(Ls>0));

% Detect rectangular outer boundary (axis-aligned) via min/max of x,y
xmin = min(X(:,1)); xmax = max(X(:,1));
ymin = min(X(:,2)); ymax = max(X(:,2));
BT = opt.BoundaryTol;
boundaryMask = abs(X(:,1)-xmin)<=BT | abs(X(:,1)-xmax)<=BT | ...
               abs(X(:,2)-ymin)<=BT | abs(X(:,2)-ymax)<=BT;
interiorMask = ~boundaryMask;  % move only these

% Precompute baseline Jacobians
[minJac0, jacAll0] = min_detJ_all(X, T, type);
if isempty(opt.MinJac)
    medPos = median(jacAll0(jacAll0>0));
    opt.MinJac = max(1e-6*max(medPos,eps), 1e-12);
end

if any(jacAll0<=0)
    warning('Some elements are already inverted or singular before distortion (min detJ = %g). Proceeding anyway.', minJac0);
end

% Random displacement field (smoothed by averaging over neighbors)
if ~isempty(opt.Seed); rng(opt.Seed); end
% unit random directions per node
R = randn(nNode,2);
R = R ./ max(hypot(R(:,1),R(:,2)), eps);
% smooth with n iterations of neighbor averaging to avoid wild oscillations
nSmooth = 0;
S = R;
for k=1:nSmooth
    S = (S + spmv_rownorm(adj, S))/2; % lightweight smoothing
end
% scale by local length and global Level
%-----Changes-----%
% beta = 0.6;
% alpha = (1-beta) + beta*opt.Level; % map Level ∈ [0,1] → alpha ∈ [1-beta,1.0]
% mag = alpha * beta*Ls;             % conservative fraction of local size
beta = 0.24;
alpha = beta*opt.Level;     % map Level ∈ [0,1] → alpha ∈ [0,beta]
mag = alpha*Ls;             % conservative fraction of local size
%-----Changes-----%
D = S .* mag;               % proposed displacements
D(~interiorMask,:) = 0;     % keep outer boundary fixed

% Backtracking to enforce positive Jacobians everywhere
Xtry = X + D;
% [minJac, ~, bad] = min_detJ_all(Xtry, T, type, opt.MinJac);
[~, ~, bad] = min_detJ_all(Xtry, T, type, opt.MinJac);
scale = 1; nbt = 0;
while any(bad) && nbt < opt.MaxBacktrack
    %-----Changes-----%
    % scale = scale * 0.50;
    scale = scale * (1-nbt/100);
    %-----Changes-----%
    Xtry = X + scale*D;
    % [minJac, ~, bad] = min_detJ_all(Xtry, T, type, opt.MinJac);
    [~, ~, bad] = min_detJ_all(Xtry, T, type, opt.MinJac);
    nbt = nbt + 1;
end

if any(bad)
    warning('Reached MaxBacktrack=%d with some elements near singular. Consider reducing Level.', opt.MaxBacktrack);
end

Xnew = X + scale*D;
% [minJac1, jacAll1] = min_detJ_all(Xnew, T, type);
[minJac1, ~] = min_detJ_all(Xnew, T, type);

if opt.Verbose
    fprintf('distort\t type=%s nodes=%d elems=%d\n', type, nNode, size(T,1));
    fprintf('  Moved %d interior nodes (of %d).\n', nnz(interiorMask), nNode);
    fprintf('  Backtracks = %d, final scale = %.4g\n', nbt, scale);
    fprintf('  min detJ: before = %.4g, after = %.4g (target MinJac = %.4g)\n', minJac0, minJac1, opt.MinJac);
end

stats = struct('scale',scale,'nBacktracks',nbt,'minJacBefore',minJac0, ...
               'minJacAfter',minJac1,'movedMask',interiorMask);
end

% ----------------- helpers -----------------
function Y = spmv_rownorm(A, X)
% Row-normalized sparse mat-vec: (D^{-1} A) * X with D(ii)=deg(i), adding self-loop
n = size(A,1);
A = A + speye(n); % include self for stability
Dinv = 1./max(sum(A,2),1);
Y = (spdiags(Dinv,0,n,n) * (A*X));
end

function [minDetJ, allDetJ, badMask] = min_detJ_all(X, T, type, minJac)
if nargin<4, minJac = -inf; end
nE = size(T,1);
allDetJ = zeros(nE,1);
badMask = false(nE,1);
if type=="Q4"
    % 2x2 Gauss (±1/sqrt(3))
    g = 1/sqrt(3);
    gps = [-g -g; +g -g; +g +g; -g +g];
    for e=1:nE
        Xe = X(T(e,:),:);
        minj = inf;
        for q=1:4
            xi = gps(q,1); eta = gps(q,2);
            [dNdxi, dNdeta] = dN_Q4(xi,eta);
            J = [dNdxi'*Xe(:,1), dNdxi'*Xe(:,2); dNdeta'*Xe(:,1), dNdeta'*Xe(:,2)];
            dj = det(J);
            if dj < minj, minj = dj; end
        end
        allDetJ(e) = minj;
        badMask(e) = minj <= minJac;
    end
elseif type=="Q8"
    % 3x3 Gauss for safety
    g = [-sqrt(3/5), 0, +sqrt(3/5)];
    gps = combvec(g,g)'; % requires Neural Network TB in old MATLAB; roll our own:
    if size(gps,1)~=9
        gps = zeros(9,2); k=0; for a=1:3, for b=1:3, k=k+1; gps(k,:)=[g(a),g(b)]; end, end
    end
    for e=1:nE
        Xe = X(T(e,:),:);
        minj = inf;
        for q=1:size(gps,1)
            xi = gps(q,1); eta = gps(q,2);
            [dNdxi, dNdeta] = dN_Q8(xi,eta);
            J = [dNdxi'*Xe(:,1), dNdxi'*Xe(:,2); dNdeta'*Xe(:,1), dNdeta'*Xe(:,2)];
            dj = det(J);
            if dj < minj, minj = dj; end
        end
        allDetJ(e) = minj;
        badMask(e) = minj <= minJac;
    end
else
    error('Unsupported type');
end
minDetJ = min(allDetJ);
end

function [dNdxi, dNdeta] = dN_Q4(xi,eta)
% Q4 bilinear derivatives wrt (xi,eta) at the given point
% Node order assumed: [(-1,-1),(+1,-1),(+1,+1),(-1,+1)]
dNdxi = 0.25*[-(1-eta); +(1-eta); +(1+eta); -(1+eta)];
dNdeta= 0.25*[-(1-xi);  -(1+xi); +(1+xi); +(1-xi)];
end

function [dNdxi, dNdeta] = dN_Q8(xi,eta)
% Q8 serendipity (8-node) derivatives wrt (xi,eta)
% Node order assumed:
% 1: (-1,-1), 2: (1,-1), 3: (1,1), 4: (-1,1),
% 5: (0,-1), 6: (1,0), 7: (0,1), 8: (-1,0)
% (mid-side nodes 5-8)

% Shape function derivatives (see standard serendipity Q8)
% Corner nodes
N1_xi  =  0.25*(1-eta)*(eta + 2*xi);
N1_eta =  0.25*(1-xi)*(2*eta + xi);
N2_xi  =  0.25*(1-eta)*(2*xi - eta);
N2_eta =  0.25*(1+xi)*(2*eta - xi);
N3_xi  =  0.25*(1+eta)*(eta + 2*xi);
N3_eta =  0.25*(1+xi)*(xi + 2*eta);
N4_xi  =  0.25*(1+eta)*(2*xi - eta);
N4_eta =  0.25*(1-xi)*(2*eta - xi);

% Mid-side nodes
N5_xi  =  xi*(eta - 1);
N5_eta =  0.5*(xi^2 - 1);
N6_xi  =  0.5*(1 - eta^2);
N6_eta = -eta*(1 + xi);
N7_xi  = -xi*(1 + eta);
N7_eta =  0.5*(1 - xi^2);
N8_xi  =  0.5*(eta^2 - 1);
N8_eta =  eta*(xi - 1);

dNdxi  = [N1_xi; N2_xi; N3_xi; N4_xi; N5_xi; N6_xi; N7_xi; N8_xi];
dNdeta = [N1_eta;N2_eta;N3_eta;N4_eta;N5_eta;N6_eta;N7_eta;N8_eta];
end