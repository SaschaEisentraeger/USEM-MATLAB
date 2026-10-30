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
[minJac0, jacAll0] = MinJacobianDeterminant(X, T, type);
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
    S = (S + RowNormalizedProduct(adj, S))/2; % lightweight smoothing
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
% [minJac, ~, bad] = MinJacobianDeterminant(Xtry, T, type, opt.MinJac);
[~, ~, bad] = MinJacobianDeterminant(Xtry, T, type, opt.MinJac);
scale = 1; nbt = 0;
while any(bad) && nbt < opt.MaxBacktrack
    %-----Changes-----%
    % scale = scale * 0.50;
    scale = scale * (1-nbt/100);
    %-----Changes-----%
    Xtry = X + scale*D;
    % [minJac, ~, bad] = MinJacobianDeterminant(Xtry, T, type, opt.MinJac);
    [~, ~, bad] = MinJacobianDeterminant(Xtry, T, type, opt.MinJac);
    nbt = nbt + 1;
end

if any(bad)
    warning('Reached MaxBacktrack=%d with some elements near singular. Consider reducing Level.', opt.MaxBacktrack);
end

Xnew = X + scale*D;
% [minJac1, jacAll1] = MinJacobianDeterminant(Xnew, T, type);
[minJac1, ~] = MinJacobianDeterminant(Xnew, T, type);

if opt.Verbose
    fprintf('distort\t type=%s nodes=%d elems=%d\n', type, nNode, size(T,1));
    fprintf('  Moved %d interior nodes (of %d).\n', nnz(interiorMask), nNode);
    fprintf('  Backtracks = %d, final scale = %.4g\n', nbt, scale);
    fprintf('  min detJ: before = %.4g, after = %.4g (target MinJac = %.4g)\n', minJac0, minJac1, opt.MinJac);
end

stats = struct('scale',scale,'nBacktracks',nbt,'minJacBefore',minJac0, ...
               'minJacAfter',minJac1,'movedMask',interiorMask);
end