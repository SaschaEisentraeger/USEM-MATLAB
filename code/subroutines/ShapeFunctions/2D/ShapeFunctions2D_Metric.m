%% Metric shape functions and derivatives (based on skew coordinates)
%-------------------------------------------------------------------------%
%
% Sascha Eisenträger, Otto von Guericke University, 15.03.2024
%
% Metric shape functions (i.e., shape functions that are expressed in
% global coordinates) are constructed including their derivatives. These
% shape functions generally fulfill completeness requirements, but violate
% compatibility.
% To ensure frame invariance and isotropy of the element properties, the
% the metric shape functions are formulated in non-dimensional corotational
% skew coordinates.
%
% Reference:
% [1] Rajendran, S. / Liew, K. M.: 
% A Novel Unsymmetric 8-Node Plane Element Immune to Mesh Distortion Under
% a Quadratic Displacement Field, 2003, International Journal for Numerical
% Methods in Engineering, Vol. 58, No. 11, Wiley, p. 1713-1748
%
% [2] Ooi, E. T.: 
% Unsymmetric Finite Elements for Mesh-Distortion Tolerance Performance,
% 2006, Nanyang Technological University, Nanyang Technological University
%
% [3] Xie, Q. / Sze, K. Y. / Zhou, Y. X.:
% Modified and Trefftz Unsymmetric Finite Element Models, 2014,
% International Journal of Mechanics and Materials in Design, Vol. 12,
% No. 1, Springer Science and Business Media LLC, p. 53-70
%
%-------------------------------------------------------------------------%
%% Input
%
% (xQ,yQ):      Global coordinates of the query point(s)
% (xN,yN):      Global coordinates of the nodes
% p:            Polynomial degree template
% p.pxi:     	Polynomial degree in x
% p.peta:     	Polynomial degree in y
% p.ndist:      Nodal distribution (reference nodes of the 'lagrange' basis)
% p.metricBasis:Optional, basis of the generalized Vandermonde matrix in
%               skew coordinates (method = 1): 'lagrange' (default),
%               'legendre', 'chebyshev' or 'monomial' (see PolynomialBasis1D)
% p.metricMethod: Optional, default for method (1 if absent)
% method:       Optional, 0: centroid-shifted global coordinates (always
%               monomials), 1: corotational skew coordinates (default),
%               2: hybrid space P_p(x,y) + {L_a(xi)L_b(eta): a+b > p}
%               (see metricHybrid below); empty: p.metricMethod
% xiQ:          Optional (method = 2), local coordinates [xi, eta] of the
%               query points; computed by Newton's method if absent
%-------------------------------------------------------------------------%
%% Output
%
% N:    Function value of the shape functions at (xQ,yQ)
% dNx:  Function value of the derivatives of the shape functions w.r.t. x
%       at (xQ,yQ)
% dNy:  Function value of the derivatives of the shape functions w.r.t. y
%       at (xQ,yQ)
% V:    Generalized Vandermonde matrix (transpose), V(c,k) = phi_c(x_k)
%-------------------------------------------------------------------------%
function [N,dNx,dNy,V] = ShapeFunctions2D_Metric(xQ,yQ,xN,yN,p,method,xiQ)

if nargin < 6 || isempty(method)
    % Global (Cartesian) coordinates: 0
    % Corotational skew coordinates:  1
    % Hybrid space:                   2
    if isfield(p,'metricMethod')
        method = p.metricMethod;
    else
        method = 1;
    end
end

if method == 2
    if nargin < 7 || isempty(xiQ)
        xiQ = Global2LocalNewton(xQ,yQ,xN,yN,p);
    end
    [N,dNx,dNy,V] = metricHybrid(xQ(:),yQ(:),xN(:),yN(:),p,xiQ);
    return
end

% Polynomial degree
px = p.pxi;
py = p.peta;

% Number of nodes (element)
nN = (px+1)*(py+1);

% Number of query points
nQ = numel(xQ);

if method == 1
    % Non-dimensional skew coordinates (query points)
    [xQ,yQ,J] = CorotationalSkewCoordinates(xQ,yQ,xN,yN,p);

    % Non-dimensional skew coordinates (nodes)
    [xN,yN] = CorotationalSkewCoordinates(xN,yN,xN,yN,p);
else
	% Standard definition according to Refs. [1,2]
	% Centroid of the element
	xC = mean(xN);
	yC = mean(yN);

	% Nodal coordinates relative to centroid
	xN = xN - repmat(xC,nN,1);
	yN = yN - repmat(yC,nN,1);

	% Query point coordinates relative to centroid
	xQ = xQ - repmat(xC,nQ,1);
	yQ = yQ - repmat(yC,nQ,1);
end

% Polynomial basis of the metric space Q_p. The shape functions do not
% depend on this choice (in exact arithmetic), but the condition number of
% the generalized Vandermonde matrix does: it grows exponentially with p for
% monomials and is 1 for the Lagrange basis of the reference element on
% affine elements (see PolynomialBasis1D)
if method == 1
    if isfield(p,'metricBasis')
        basis = p.metricBasis;
    else
        basis = 'lagrange';
    end
else
    % Dimensional global coordinates: the reference interval [-1,1] has
    % no meaning, keep the monomials of Refs. [1,2]
    basis = 'monomial';
end
if isfield(p,'ndist')
    ndist = p.ndist;
else
    ndist = 2;
end

% The Legendre, Chebyshev and Lagrange bases are defined on [-1,1]. On
% distorted elements the skew coordinates of the nodes leave [-1,1]^2,
% where these polynomials grow quickly; map the bounding box of the nodes
% onto [-1,1]^2 (affine, so neither the space Q_p nor frame invariance is
% affected; identity for parallelograms)
sx = 1;
sy = 1;
if ~strcmpi(basis,'monomial')
    xc = (max(xN) + min(xN))/2;
    yc = (max(yN) + min(yN))/2;
    sx = (max(xN) - min(xN))/2;
    sy = (max(yN) - min(yN))/2;
    xN = (xN - xc)/sx;
    yN = (yN - yc)/sy;
    xQ = (xQ - xc)/sx;
    yQ = (yQ - yc)/sy;
end

% One-dimensional bases at the nodes and at the query points
[Vx,~]    = PolynomialBasis1D(xN,px,basis,ndist);
[Vy,~]    = PolynomialBasis1D(yN,py,basis,ndist);
[Qx,dQx]  = PolynomialBasis1D(xQ,px,basis,ndist);
[Qy,dQy]  = PolynomialBasis1D(yQ,py,basis,ndist);

% Generalized Vandermonde matrix (transpose): V(c,k) = phi_c(x_k), with
% phi_c(x,y) = phi_j(x)*phi_i(y), c = i*(px+1) + j + 1
V = kronRows(Vy,Vx);

% Polynomial vector/matrix for the query point and its derivatives
P  = kronRows(Qy,Qx);
Px = kronRows(Qy,dQx);
Py = kronRows(dQy,Qx);

% Stop if the nodes are not unisolvent for the polynomial space
checkUnisolvent(V,px,py);

% Shape functions
dV = decomposition(V);
N  = dV\P;

if method == 1
    % Derivatives w.r.t. to the skew corotational coordinates
    dNxi = (dV\Px)/sx;
    dNeta = (dV\Py)/sy;
    % Derivatives w.r.t. to the global coordinates
    gradN = J\[reshape(dNxi,1,[]);reshape(dNeta,1,[])];
    dNx = reshape(gradN(1,:),nN,nQ);
    dNy = reshape(gradN(2,:),nN,nQ);
else
    % Derivatives w.r.t. to the global coordinates
    dNx = dV\Px;
    dNy = dV\Py;
end

end

%% Hybrid trial space (method = 2)
% T_p = P_m(x,y) + span{L_a(xi) L_b(eta) : a <= px, b <= py, a + b > m},
% m = min(px,py), L_a: Legendre polynomials. The complete polynomials are
% written in the bounding-box scaled skew coordinates (affine in x, so the
% space is P_m(x,y)), the complement in the parametric coordinates, in
% which the nodes form a tensor grid. T_p contains P_m(x,y) (completeness,
% patch tests), is frame invariant and equals Q_p on affine elements. The
% complement depends on the choice of L_a, so the basis is fixed here
% (p.metricBasis is ignored).
function [N,dNx,dNy,V] = metricHybrid(xQ,yQ,xN,yN,p,xiQ)

px = p.pxi;
py = p.peta;
m  = min(px,py);

% Skew coordinates, bounding box of the nodes mapped onto [-1,1]^2
[sQ,tQ,J0] = CorotationalSkewCoordinates(xQ,yQ,xN,yN,p);
[sN,tN]    = CorotationalSkewCoordinates(xN,yN,xN,yN,p);
cs = (max(sN) + min(sN))/2;  rs = (max(sN) - min(sN))/2;
ct = (max(tN) + min(tN))/2;  rt = (max(tN) - min(tN))/2;
sN = (sN - cs)/rs;  sQ = (sQ - cs)/rs;
tN = (tN - ct)/rt;  tQ = (tQ - ct)/rt;

% Local coordinates of the nodes and of the query points
[xiN,etaN] = ReferenceNodes2D(p);
xiQq  = xiQ(:,1);
etaQq = xiQ(:,2);

% One-dimensional Legendre polynomials
[Ls,~]      = PolynomialBasis1D(sN,px,'legendre');
[Lt,~]      = PolynomialBasis1D(tN,py,'legendre');
[Lxi,~]     = PolynomialBasis1D(xiN,px,'legendre');
[Leta,~]    = PolynomialBasis1D(etaN,py,'legendre');
[Qs,dQs]    = PolynomialBasis1D(sQ,px,'legendre');
[Qt,dQt]    = PolynomialBasis1D(tQ,py,'legendre');
[Qxi,dQxi]  = PolynomialBasis1D(xiQq,px,'legendre');
[Qeta,dQeta]= PolynomialBasis1D(etaQq,py,'legendre');

% Jacobian of the geometry mapping at the query points
[~,dNxi,dNeta] = ShapeFunctions2D_SEM(xiQq,etaQq,p);
xxi  = (dNxi'*xN)';   yxi  = (dNxi'*yN)';
xeta = (dNeta'*xN)';  yeta = (dNeta'*yN)';
detJ = xxi.*yeta - xeta.*yxi;

% The functions are evaluated through (xi,eta), so the sign of det J does
% not matter (completeness holds pointwise in (xi,eta), also on folded
% elements); the gradients of the complement need J^-1 at the query points
if any(~(abs(detJ) > 0))
    error('ShapeFunctions2D_Metric:singularJacobian', ...
         ['det J = 0 at a query point: the gradients of the hybrid ' ...
          'metric shape functions (method = 2) are undefined there.']);
end

% Inverse of J0 (skew coordinates are affine: grad_x = J0^-1 grad_(s,t))
J0inv = inv(J0);

nN = (px+1)*(py+1);
nQ = numel(xQ);
V  = zeros(nN);
P  = zeros(nN,nQ);
Px = zeros(nN,nQ);
Py = zeros(nN,nQ);
c = 0;
for b = 0:py
    for a = 0:px
        c = c + 1;
        if a + b <= m
            % Complete polynomials in (x,y)
            V(c,:) = Ls(a+1,:).*Lt(b+1,:);
            P(c,:) = Qs(a+1,:).*Qt(b+1,:);
            Ps = dQs(a+1,:).*Qt(b+1,:)/rs;
            Pt = Qs(a+1,:).*dQt(b+1,:)/rt;
            Px(c,:) = J0inv(1,1)*Ps + J0inv(1,2)*Pt;
            Py(c,:) = J0inv(2,1)*Ps + J0inv(2,2)*Pt;
        else
            % Complement in the parametric coordinates
            V(c,:) = Lxi(a+1,:).*Leta(b+1,:);
            P(c,:) = Qxi(a+1,:).*Qeta(b+1,:);
            Pxi  = dQxi(a+1,:).*Qeta(b+1,:);
            Peta = Qxi(a+1,:).*dQeta(b+1,:);
            Px(c,:) = ( yeta.*Pxi - yxi.*Peta)./detJ;
            Py(c,:) = (-xeta.*Pxi + xxi.*Peta)./detJ;
        end
    end
end

checkUnisolvent(V,px,py);

dV  = decomposition(V);
N   = dV\P;
dNx = dV\Px;
dNy = dV\Py;

end

%% Unisolvence check
% The metric shape functions are undefined if the nodes are not unisolvent
% for the polynomial space (e.g. coincident nodes or strongly distorted
% elements): stop instead of returning meaningless shape functions
function checkUnisolvent(V,px,py)
rc = rcond(V);
if ~(rc >= eps)
    error('ShapeFunctions2D_Metric:singular', ...
         ['The interpolation matrix of the metric shape functions is ' ...
          'singular to working precision (rcond = %.2e): the element ' ...
          'nodes are (nearly) not unisolvent. The element is too ' ...
          'strongly distorted for the metric shape functions of degree ' ...
          '(%d,%d).'], rc, px, py);
end
end

%% Column-wise Kronecker product: C(i*nB+j+1,:) = A(i+1,:).*B(j+1,:)
function C = kronRows(A,B)
nA = size(A,1);
nB = size(B,1);
C = reshape(reshape(B,nB,1,[]).*reshape(A,1,nA,[]),nA*nB,[]);
end
%------------------------------Visualization------------------------------%
% % Grid points
% nGP = 101;
% xi = linspace(-1,1,nGP);
% eta= linspace(-1,1,nGP);
% 
% % Save integration points in a vector
% [Xi,Eta] = meshgrid(xi,eta);
% IP = [reshape(Xi,[],1),reshape(Eta,[],1)];
% 
% % Nodal coordinates
% nE = 1;
% xE = xCoor(Topology_Node(nE,:));
% yE = yCoor(Topology_Node(nE,:));
% 
% % Global coordinates of the integration points
% IPg = Local2GlobalCoordinatesQuad(p,IP,[xE,yE]);
% 
% % Shape functions and derivatives
% [M,dMx,dMy] = ShapeFunctions2D_Metric(IPg(:,1),IPg(:,2),xE,yE,p);
% 
% % Visualization
% nSF = 1;
% f1 = figure('Name',['Metric Shape Function #' num2str(nSF)]);
% f1.WindowState = 'maximized';
% set(0, 'CurrentFigure', f1)
% set(gca,'FontName','Times New Roman','FontSize',24);
% grid on
% box on
% hold on
% surf(reshape(IP(:,1),nGP,nGP), reshape(IP(:,2),nGP,nGP), reshape(M(nSF,:),nGP,nGP))
% plot3([-1,1,1,-1,-1],[-1,-1,1,1,-1],[0,0,0,0,0],'r')
% % surf(reshape(IPg(:,1),nGP,nGP), reshape(IPg(:,2),nGP,nGP), reshape(M(nSF,:),nGP,nGP))
% % plot3(xCoor([1,2,3,4,1]),yCoor([1,2,3,4,1]),[0,0,0,0,0],'r')
% view(3)
%-------------------------------------------------------------------------%