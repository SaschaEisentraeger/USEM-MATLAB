function [Xnew, stats] = distort_single_element(X, T, type, varargin)
%DISTORT_SINGLE_ELEMENT Distorts one interior element in a 2D FE mesh with backtracking.
%
%   Boundary nodes are not moved.

% ---------------- Parse input ----------------
p = inputParser;
p.addRequired('X', @(v)isnumeric(v)&&size(v,2)==2);
p.addRequired('T', @(v)isnumeric(v)&&ismatrix(v));
p.addRequired('type', @(s)ischar(s)||isstring(s));
p.addParameter('Level',0.35,@(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('Seed',[],@(v)isempty(v)||(isscalar(v)&&isfinite(v)));
p.addParameter('MinJac',[],@(v)isempty(v)||(isscalar(v)&&v>0));
%-----Changes-----%
% p.addParameter('MaxBacktrack', 12, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('MaxBacktrack', 100, @(v)isnumeric(v)&&isscalar(v)&&v>=0);
%-----Changes-----%
p.addParameter('BoundaryTol',1e-12,@(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('Verbose',true,@(v)islogical(v)&&isscalar(v));
p.parse(X,T,type,varargin{:});
opt = p.Results; type = upper(string(type));

Xnew = X;
nNode = size(X,1);
nElem = size(T,1);
nen = size(T,2);

if (type=="Q4" && nen~=4) || (type=="Q8" && nen~=8)
    error('Connectivity does not match element type');
end

% ---------------- Select element by centroid ----------------
elemCentroids = zeros(nElem,2);
for e=1:nElem
    elemCentroids(e,:) = mean(X(T(e,:),:),1);
end
meshCenter = mean(X,1);
[~, midElem] = min(vecnorm(elemCentroids - meshCenter,2,2));
elementNodes = T(midElem,:);
elementMask = false(nNode,1);
elementMask(elementNodes) = true;

% ---------------- Detect boundary nodes ----------------
xmin = min(X(:,1)); xmax = max(X(:,1));
ymin = min(X(:,2)); ymax = max(X(:,2));
BT = opt.BoundaryTol;
boundaryMask = abs(X(:,1)-xmin)<=BT | abs(X(:,1)-xmax)<=BT | ...
               abs(X(:,2)-ymin)<=BT | abs(X(:,2)-ymax)<=BT;

% Remove boundary nodes from elementMask
elementMask(boundaryMask) = false;

% ---------------- Characteristic length ----------------
coords = X(elementMask,:);
if isempty(coords)
    warning('All element nodes are on the boundary. No displacement applied.');
    stats = struct('scale',0,'nBacktracks',0,'minJacBefore',NaN,'minJacAfter',NaN,'movedMask',elementMask,'elementNodes',elementNodes);
    return;
end
%-----Changes-----%
% Lchar = mean(sqrt(sum((coords - mean(coords)).^2,2)));
Lchar = 2* max([max(coords(:,1))-min(coords(:,1)),...
                max(coords(:,2))-min(coords(:,2))]);
%-----Changes-----%

% ---------------- Random displacement ----------------
if ~isempty(opt.Seed); rng(opt.Seed); end
R = randn(sum(elementMask),2);
R = R ./ max(hypot(R(:,1),R(:,2)), eps);

%-----Changes-----%
% beta = 0.6;
% alpha = (1-beta) + beta*opt.Level; % map Level ∈ [0,1] → alpha ∈ [1-beta,1.0]
% mag = alpha * beta*Lchar;          % conservative fraction of local size
beta = 0.24;
alpha = beta*opt.Level;     % map Level ∈ [0,1] → alpha ∈ [0,beta]
mag = alpha*Lchar;             % conservative fraction of local size
%-----Changes-----%
D = zeros(nNode,2);
D(elementMask,:) = R .* mag;

% ---------------- Minimum Jacobian ----------------
[minJac0, jacAll0] = MinJacobianDeterminant(X, T, type);
if isempty(opt.MinJac)
    medPos = median(jacAll0(jacAll0>0));
    opt.MinJac = max(1e-6*max(medPos,eps),1e-12);
end

% ---------------- Backtracking ----------------
Xtry = X + D;
[~, ~, bad] = MinJacobianDeterminant(Xtry, T, type, opt.MinJac);
scale = 1; nbt = 0;
while any(bad) && nbt < opt.MaxBacktrack
    %-----Changes-----%
    % scale = scale * 0.50;
    scale = scale * (1-nbt/100);
    %-----Changes-----%
    Xtry = X + scale*D;
    [~, ~, bad] = MinJacobianDeterminant(Xtry, T, type, opt.MinJac);
    nbt = nbt + 1;
end
if any(bad) && opt.Verbose
    warning('Reached MaxBacktrack with some elements near singular.');
end
Xnew = X + scale*D;
[minJac1, ~] = MinJacobianDeterminant(Xnew, T, type);

% ---------------- Output stats ----------------
stats = struct('scale',scale, 'nBacktracks',nbt, ...
               'minJacBefore',minJac0,'minJacAfter',minJac1, ...
               'movedMask',elementMask, 'elementNodes',elementNodes);

if opt.Verbose
    fprintf('Distorted element %d with %d interior nodes\n', midElem, nnz(elementMask));
    fprintf('Backtracks = %d, final scale = %.4g\n', nbt, scale);
    fprintf('min detJ before = %.4g, after = %.4g (target MinJac=%.4g)\n', ...
        minJac0, minJac1, opt.MinJac);
end
end