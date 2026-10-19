%% Element stiffness matrix for the unsymmetric SEM
%-------------------------------------------------------------------------%
%
% Sascha Eisenträger, Otto von Guericke University, 13.03.2024
%
% Numerical computation of the element stiffness matrix of a two-
% dimensional unsymmetric spectral element.
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
% C:            Elasticity matrix
% xCoor:        x coordinates of the nodes
% yCoor:        y coordinates of the nodes
% nIPxi:		Number of integration points for the xi-direction
% nIPeta:		Number of integration points for the xi-direction
%-------------------------------------------------------------------------%
%% Output
%
% K_element:   	Element stiffness matrix
%-------------------------------------------------------------------------%
function [K_element] = Ke_2D(p,C,xCoor,yCoor,nIPxi,nIPeta)

if nargin == 4
    % Polynomial degrees in xi and eta
    pxi = p.pxi;
    peta = p.peta;
    % Number of integration points
    d = 2;
    nIPxi = d*pxi;
    nIPeta = d*peta;
    % nIPxi = pxi + 1;
    % nIPeta = peta + 1;
    % nIPxi = 6;
    % nIPeta = 6;
end

% Number of degrees of freedom per element
ndof = 2*numel(xCoor);

% Element stiffness matrix
K_element = zeros(ndof,ndof);

% Gauss-Legendre points (numerical integration)
[xi,wxi]   = lgwt(nIPxi,-1,1);
[eta,weta] = lgwt(nIPeta,-1,1);

% Save integration points in a vector
[Xi,Eta] = meshgrid(xi,eta);
IP = [reshape(Xi,[],1),reshape(Eta,[],1)];
IP = flipud(IP);

% Number of integration points
NumIP = numel(IP(:,1));

% Save integration weights in a vector
[Wxi,Weta] = meshgrid(wxi,weta);
IW = reshape(Wxi,[],1).*reshape(Weta,[],1);
IW = flipud(IW);

% Derivatives of the parametric shape functions (xi,eta)
[~,dNxi,dNeta] = ShapeFunctions2D_SEM(IP(:,1),IP(:,2),p);

% Derivatives of the metric shape functions (x,y)
IPg = Local2GlobalCoordinatesQuad(p,IP,[xCoor,yCoor]);
[~,dMx,dMy] = ShapeFunctions2D_Metric(IPg(:,1),IPg(:,2),xCoor,yCoor,p);
    
% Loop over all integration points
for i = 1:NumIP
    % Jacobi matrix (geometry mapping)
    [J] = Jacobi_Matrix(p,IP(i,:),xCoor,yCoor);

    % Adjoint of the Jacobi matrix
    Jadj = [J(2,2), -J(1,2); -J(2,1), J(1,1)];

    % B-Matrix: Strain-displacement matrix
    % Derivative of the parametric shape function w.r.t. x and y
    dN = Jadj*[dNxi(:,i)';dNeta(:,i)'];
    dNx = dN(1,:);
    dNy = dN(2,:);
	% BL: parametric shape functions
    BL = zeros(3,ndof);
    BL(1,1:2:ndof) = dNx;
	BL(2,2:2:ndof) = dNy;
	BL(3,1:2:ndof) = dNy;
	BL(3,2:2:ndof) = dNx;
	% BR: metric shape functions
	BR = zeros(3,ndof);
	BR(1,1:2:ndof) = dMx(:,i);
	BR(2,2:2:ndof) = dMy(:,i);
	BR(3,1:2:ndof) = dMy(:,i);
	BR(3,2:2:ndof) = dMx(:,i);

    % Contribution to the element stiffness matrix at the ith integration
	% point
	K_element = K_element + BL'*C*BR*IW(i);
end

end