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
%               2: hybrid space P_m(x,y) + {L_a(xi)L_b(eta): a+b > m},
%               m = min(pxi,peta)
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
    [N,dNx,dNy,V] = ShapeFunctions2D_MetricHybrid(xQ(:),yQ(:),xN(:),yN(:),p,xiQ);
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
V = KronRows(Vy,Vx);

% Polynomial vector/matrix for the query point and its derivatives
P  = KronRows(Qy,Qx);
Px = KronRows(Qy,dQx);
Py = KronRows(dQy,Qx);

% Stop if the nodes are not unisolvent for the polynomial space
CheckUnisolvent(V,px,py);

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