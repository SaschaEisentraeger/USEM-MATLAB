%% Element mass matrix for the unsymmetric SEM
%-------------------------------------------------------------------------%
%
% Sascha Eisenträger, Otto von Guericke University, 13.03.2024
%
% Numerical computation of the element mass matrix of a two-dimensional
% unsymmetric spectral element.
%
% Unsymmetric FEM - Reference:
% Rajendran, S.:
% A Technique to Develop Mesh-Distortion Immune Finite Elements, 2010, 
% Computer Methods in Applied Mechanics and Engineering, Vol. 199,
% No. 17–20, Elsevier BV, p. 1044-1063
%
%-------------------------------------------------------------------------%
%% Input
%
% p:            Polynomial degree structure
% p.pxi:      	Polynomial degree in xi
% p.peta:      	Polynomial degree in eta
% p.ndist:   	Nodal distribution
% rho:         	Mass density
% xCoor:        x coordinates of the nodes
% yCoor:        y coordinates of the nodes
% lumping:     	Mass matrix flags 
%             	1: Consistent mass matrix
%              	2: Lumped mass matrix (nodal quadrature)
%              	3: Lumped mass matrix (row-sum method)
%             	4: Lumped mass matrix (HRZ-method/diagonal scaling)
% nIPxi:		Number of integration points for the xi-direction
% nIPeta:		Number of integration points for the xi-direction
%-------------------------------------------------------------------------%
%% Output
%
% M_element:   	Element mass matrix
%-------------------------------------------------------------------------%
function [M_element] = Me_2D(p,rho,xCoor,yCoor,lumping,nIPxi,nIPeta)

% Polynomial degrees in xi and eta
pxi = p.pxi;
peta = p.peta;

if nargin == 5
    % Number of integration points
    d = 2;
    nIPxi = d*pxi;
    nIPeta = d*peta;
    % nIPxi = pxi + 1;
    % nIPeta = peta + 1;
    % nIPxi = 10;
    % nIPeta = 10;
end

% Distribution of the nodes
ndist = p.ndist;

% Number of degrees of freedom per element
ndof = 2*numel(xCoor);

% Element mass matrix
M_element = zeros(ndof,ndof);

% Nodal quadrature is only used for spectral elements on GLL points (only
% in this case the nodal quadrature achieves optimal results)
if lumping == 2 && ndist == 1
    % Use the HRZ-method for elements with equidistant points
    lumping = 4;
elseif lumping == 2 && ndist == 3
    % Use the row-sum method for elements with GLC points
    lumping = 3;
end

if lumping == 2
    % Gauss-Lobatto points (numerical integration)
    [xi,wxi,~]   = lglnodes(pxi);
    [eta,weta,~] = lglnodes(peta);
else
    % Gauss-Legendre points (numerical integration)
    [xi,wxi]   = lgwt(nIPxi,-1,1);
    [eta,weta] = lgwt(nIPeta,-1,1);
end

% Save integration points in a vector
[Xi,Eta] = meshgrid(xi,eta);
IP = [reshape(Xi,[],1),reshape(Eta,[],1)];

% Number of integration points
NumIP = numel(IP(:,1));

% Save integration weights in a vector
[Wxi,Weta] = meshgrid(wxi,weta);
IW = reshape(Wxi,[],1).*reshape(Weta,[],1);

% Identity matrix
I = eye(2);

% Parametric shape functions (xi,eta)
Naux = ShapeFunctions2D_SEM(IP(:,1),IP(:,2),p);

% Metric shape functions (x,y)
IPg = Local2GlobalCoordinatesQuad(p,IP,[xCoor,yCoor]);
Maux = ShapeFunctions2D_Metric(IPg(:,1),IPg(:,2),xCoor,yCoor,p);

% Loop over all integration points
for i = 1:NumIP
    % Shape function matrix (parametric)
    N = kron(Naux(:,i)',I);
	
	% Shape function matrix (metric)
    M = kron(Maux(:,i)',I);

    % Jacobi matrix (geometry mapping)
    [J] = Jacobi_Matrix(p,IP(i,:),xCoor,yCoor);
    detJ = det(J);
	
	% Contribution to the element mass matrix at the ith integration point
	M_element = M_element + rho*(N'*M)*detJ*IW(i);
end

% Nodal quadrature: the parametric and the metric shape functions are
% interpolatory at the GLL points, so the mass matrix is diagonal; remove
% the round-off in the off-diagonal entries
if lumping == 2
    M_element = diag(diag(M_element));
end

% Mass lumping procedures
if lumping == 3
    % Row-sum method
	[M_element] = RowSum_Method(M_element);
elseif lumping == 4
    % HRZ-method/diagonal scaling
	[M_element] = HRZ_Method(M_element);
end

end