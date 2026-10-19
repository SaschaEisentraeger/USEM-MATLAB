%% Element stiffness and mass matrices for the unsymmetric SEM
%-------------------------------------------------------------------------%
%
% Sascha Eisenträger, Otto von Guericke University, 13.03.2024
%
% Numerical computation of the element stiffness and mass matrices of a
% two-dimensional unsymmetric spectral element.
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
%-------------------------------------------------------------------------%
%% Output
%
% K_element:   	Element stiffness matrix
% M_element:   	Element mass matrix
%-------------------------------------------------------------------------%
function [K_element,M_element] = KMe_2D(p,C,rho,xCoor,yCoor,lumping)

% Polynomial degrees in xi and eta
pxi = p.pxi;
peta = p.peta;

% Distribution of the nodes
ndist = p.ndist;

% Number of degrees of freedom per element
ndof = 2*numel(xCoor);

% Element stiffness matrix
K_element = zeros(ndof,ndof);

% Gauss-Legendre points (numerical integration)
[xiK,wxiK]     = lgwt(2*pxi,-1,1);
[etaK,wetaK]   = lgwt(2*peta,-1,1);

% Save integration points in a vector - stiffness matrix
[XiK,EtaK] = meshgrid(xiK,etaK);
IPK = [reshape(XiK,[],1),reshape(EtaK,[],1)];

% Save integration weights in a vector - stiffness matrix
[WxiK,WetaK] = meshgrid(wxiK,wetaK);
IWK = reshape(WxiK,[],1).*reshape(WetaK,[],1);

% Number of integration points (stiffness matrix)
NumIP = numel(IPK(:,1));

% Derivatives of the parametric shape functions (xi,eta)
[~,dNxi,dNeta] = ShapeFunctions2D_SEM(IPK(:,1),IPK(:,2),p);

% Derivatives of the metric shape functions (x,y)
IPKg = Local2GlobalCoordinatesQuad(p,IPK,[xCoor,yCoor]);
[~,dMx,dMy] = ShapeFunctions2D_Metric(IPKg(:,1),IPKg(:,2),xCoor,yCoor,p);

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
    [xiM,wxiM,~]   = lglnodes(pxi);
    [etaM,wetaM,~] = lglnodes(peta);
else
    % Gauss-Legendre points (numerical integration): 2p points per
    % direction as in Me_2D. On distorted elements the metric shape
    % functions are not polynomials of degree p in (xi,eta), so p+1 points
    % would under-integrate the unsymmetric mass matrix.
    [xiM,wxiM]   = lgwt(2*pxi,-1,1);
    [etaM,wetaM] = lgwt(2*peta,-1,1);
end

% Identity matrix
I = eye(2);

% Save integration points in a vector - mass matrix
[XiM,EtaM] = meshgrid(xiM,etaM);
IPM = [reshape(XiM,[],1),reshape(EtaM,[],1)];

% Save integration weights in a vector - mass matrix
[WxiM,WetaM] = meshgrid(wxiM,wetaM);
IWM = reshape(WxiM,[],1).*reshape(WetaM,[],1);

% Number of integration points (mass matrix)
NumIPM = numel(IPM(:,1));

% Parametric shape functions (xi,eta)
Naux = ShapeFunctions2D_SEM(IPM(:,1),IPM(:,2),p);

% Metric shape functions (x,y)
IPMg = Local2GlobalCoordinatesQuad(p,IPM,[xCoor,yCoor]);
Maux = ShapeFunctions2D_Metric(IPMg(:,1),IPMg(:,2),xCoor,yCoor,p);
 
% Loop over all integration points
for i = 1:NumIP
	%---------------------------Stiffness Matrix--------------------------%
    % Jacobi matrix (geometry mapping)
    [JK] = Jacobi_Matrix(p,IPK(i,:),xCoor,yCoor);

    % Adjoint of the Jacobi matrix
	Jadj = [JK(2,2), -JK(1,2); -JK(2,1), JK(1,1)];

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
	K_element = K_element + BL'*C*BR*IWK(i);
	
	%-----------------------------Mass Matrix-----------------------------%
	if i <= NumIPM
		% Shape function matrix (parametric)
		N = kron(Naux(:,i)',I);
		
		% Shape function matrix (metric)
		M = kron(Maux(:,i)',I);

		% Jacobi matrix (geometry mapping)
		[JM] = Jacobi_Matrix(p,IPM(i,:),xCoor,yCoor);
		detJM = det(JM);
		
		% Contribution to the element mass matrix at the ith integration
		% point
		M_element = M_element + rho*(N'*M)*detJM*IWM(i);
	end
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