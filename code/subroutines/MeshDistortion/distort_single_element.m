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
[minJac0, jacAll0] = min_detJ_all(X, T, type);
if isempty(opt.MinJac)
    medPos = median(jacAll0(jacAll0>0));
    opt.MinJac = max(1e-6*max(medPos,eps),1e-12);
end

% ---------------- Backtracking ----------------
Xtry = X + D;
[~, ~, bad] = min_detJ_all(Xtry, T, type, opt.MinJac);
scale = 1; nbt = 0;
while any(bad) && nbt < opt.MaxBacktrack
    %-----Changes-----%
    % scale = scale * 0.50;
    scale = scale * (1-nbt/100);
    %-----Changes-----%
    Xtry = X + scale*D;
    [~, ~, bad] = min_detJ_all(Xtry, T, type, opt.MinJac);
    nbt = nbt + 1;
end
if any(bad) && opt.Verbose
    warning('Reached MaxBacktrack with some elements near singular.');
end
Xnew = X + scale*D;
[minJac1, ~] = min_detJ_all(Xnew, T, type);

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

% ----------------- helpers -----------------
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