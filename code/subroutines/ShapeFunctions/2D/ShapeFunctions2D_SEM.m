%% Shape functions and their derivative for spectral elements (2D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek, Otto von Guericke University, 09.03.2018
%
% The shape functions for a spectral element and the corresponding 
% derivatives with respect to the local coordinates xi and eta are 
% computed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% xiIP, etaIP:  Local coordinates
% p:            Polynomial degree template
% p.pxi:     	Polynomial degree in xi
% p.peta:     	Polynomial degree in eta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% N:        Function value of the shape functions at (xi,eta)
% dNxi:     Function value of the derivatives of the shape functions w.r.t
%           xi at (xi,eta)
% dNeta:    Function value of the derivatives of the shape functions w.r.t
%           eta at (xi,eta)
%------------------------------Visualization------------------------------%
% xi = (-1:0.05:1);
% eta = (-1:0.05:1);
% ndiv.xi  = numel(xi);
% ndiv.eta = numel(eta);
% % Grid data: Local coordinates
% [Xi,Eta] = meshgrid(xi,eta);
% % Save all points in a vector
% IP = [reshape(Xi,[],1),reshape(Eta,[],1)];
% 
% % Shape functions 2d
% [N,dNxi,dNeta] = ShapeFunctions2D_SEM(IP(:,1),IP(:,2),p,3);
% 
% % Number of shape functions
% n_sf = size(N,1);
% 
% for i = 1:n_sf
%     % Visualization of the shape functions
%     var = N(i,:);
%     Res = reshape(var,ndiv.xi,ndiv.eta);
%     figure('Name','Visualization of the shape functions of a xNy-element')
%     set(gcf,'Units','normalized','Position',[0 0 1 1]);
%     surf(Xi,Eta,Res);
%     axis equal
%     rotate3d on
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,dNxi,dNeta] = ShapeFunctions2D_SEM(xiIP,etaIP,p)

% Convert the input to row vectors
xiIP  = reshape(xiIP,1,[]);
etaIP = reshape(etaIP,1,[]);

% Polynomial degree
pxi = p.pxi;
peta = p.peta;

% Nodal distribution --> supporting points of the Lagrange polynomials
xi  = nodal_distribution(pxi);
eta = nodal_distribution(peta);

% The shape functions are stored as Matlab's anonymous functions in 
% Phi_mat_xi and Phi_mat_eta. See LagrangePolynomials.m for more 
% information!
[Phi_mat_xi , dPhi_mat_xi]  = LagrangePolynomials(xi);
[Phi_mat_eta, dPhi_mat_eta] = LagrangePolynomials(eta);

% Number of integration points
nIP = numel(xiIP);

% Ordering of the 1d shape functions; the two nodes at the interval 
% boundaries are assigned the first two degrees of freedom while the 
% interior degrees of freedom are number consecutively
vec_pxi  = [1 pxi+1  (2:pxi)] ;
vec_peta = [1 peta+1 (2:peta)];

% Shape function vector 1d
Phi_j  = zeros(pxi+1,nIP);
dPhi_j = zeros(pxi+1,nIP);
Phi_i  = zeros(peta+1,nIP);
dPhi_i = zeros(peta+1,nIP);

% Evaluate the 1d shape functions at the current integration point
for i = 1:pxi+1
	 Phi_j(i,:) =  Phi_mat_xi{vec_pxi(i)}(xiIP);
	dPhi_j(i,:) = dPhi_mat_xi{vec_pxi(i)}(xiIP);
end
for i = 1:peta+1
	 Phi_i(i,:) =  Phi_mat_eta{vec_peta(i)}(etaIP);
	dPhi_i(i,:) = dPhi_mat_eta{vec_peta(i)}(etaIP);
end

% Nodal modes
N_node1 = Phi_j(1,:).*Phi_i(1,:);
N_node2 = Phi_j(2,:).*Phi_i(1,:);
N_node3 = Phi_j(2,:).*Phi_i(2,:);
N_node4 = Phi_j(1,:).*Phi_i(2,:);
N_node = [N_node1; N_node2; N_node3; N_node4];

% Derivative of the nodal modes w.r.t. xi
dNxi_node1 = dPhi_j(1,:).*Phi_i(1,:);
dNxi_node2 = dPhi_j(2,:).*Phi_i(1,:);
dNxi_node3 = dPhi_j(2,:).*Phi_i(2,:);
dNxi_node4 = dPhi_j(1,:).*Phi_i(2,:);
dNxi_node = [dNxi_node1; dNxi_node2; dNxi_node3; dNxi_node4];

% Derivative of the nodal modes w.r.t. eta
dNeta_node1 = Phi_j(1,:).*dPhi_i(1,:);
dNeta_node2 = Phi_j(2,:).*dPhi_i(1,:);
dNeta_node3 = Phi_j(2,:).*dPhi_i(2,:);
dNeta_node4 = Phi_j(1,:).*dPhi_i(2,:);
dNeta_node = [dNeta_node1; dNeta_node2; dNeta_node3; dNeta_node4];

% Edge modes
% xi-Edges
if pxi > 1
    N_edge1     = zeros((pxi-1),nIP);
    N_edge3     = zeros((pxi-1),nIP);
    dNxi_edge1  = zeros((pxi-1),nIP);
    dNxi_edge3  = zeros((pxi-1),nIP);
    dNeta_edge1 = zeros((pxi-1),nIP);
    dNeta_edge3 = zeros((pxi-1),nIP);
    for k = 2:pxi
        % Edge 1
        N_edge1(k-1,:)     =  Phi_i(1,:).* Phi_j(k+1,:);
        dNxi_edge1(k-1,:)  =  Phi_i(1,:).*dPhi_j(k+1,:);
        dNeta_edge1(k-1,:) = dPhi_i(1,:).* Phi_j(k+1,:);
        % Edge 3
        N_edge3(k-1,:)     =  Phi_i(2,:).* Phi_j(k+1,:);
        dNxi_edge3(k-1,:)  =  Phi_i(2,:).*dPhi_j(k+1,:);
        dNeta_edge3(k-1,:) = dPhi_i(2,:).* Phi_j(k+1,:);
    end
else
    N_edge1     = [];
    N_edge3     = [];
    dNxi_edge1  = [];
    dNxi_edge3  = [];
    dNeta_edge1 = [];
    dNeta_edge3 = [];
end

% eta-Edges
if peta > 1
    N_edge2     = zeros((peta-1),nIP);
    N_edge4     = zeros((peta-1),nIP);
    dNxi_edge2  = zeros((peta-1),nIP);
    dNxi_edge4  = zeros((peta-1),nIP);
    dNeta_edge2 = zeros((peta-1),nIP);
    dNeta_edge4 = zeros((peta-1),nIP);
    for k = 2:peta
        % Edge 2
        N_edge2(k-1,:)     =  Phi_j(2,:).* Phi_i(k+1,:);
        dNxi_edge2(k-1,:)  = dPhi_j(2,:).* Phi_i(k+1,:);
        dNeta_edge2(k-1,:) =  Phi_j(2,:).*dPhi_i(k+1,:);
        % Edge 4
        N_edge4(k-1,:)     =  Phi_j(1,:).* Phi_i(k+1,:);
        dNxi_edge4(k-1,:)  = dPhi_j(1,:).* Phi_i(k+1,:);
        dNeta_edge4(k-1,:) =  Phi_j(1,:).*dPhi_i(k+1,:);
    end
else
    N_edge2 = [];
    N_edge4 = [];
    dNxi_edge2 = [];
    dNxi_edge4 = [];
    dNeta_edge2 = [];
    dNeta_edge4 = [];
end 
N_edge     = [N_edge1; N_edge2; N_edge3; N_edge4];
dNxi_edge  = [dNxi_edge1; dNxi_edge2; dNxi_edge3; dNxi_edge4];
dNeta_edge = [dNeta_edge1; dNeta_edge2; dNeta_edge3; dNeta_edge4];

if pxi > 1 && peta > 1
    % Internal/bubble modes
    N_int     = zeros((pxi-1)*(peta-1),nIP);
    dNxi_int  = zeros((pxi-1)*(peta-1),nIP);
    dNeta_int = zeros((pxi-1)*(peta-1),nIP);
    m = 1;
    for k = 2:peta
        for l = 2:pxi
            N_int(m,:)     =  Phi_j(l+1,:).* Phi_i(k+1,:);
            dNxi_int(m,:)  = dPhi_j(l+1,:).* Phi_i(k+1,:);
            dNeta_int(m,:) =  Phi_j(l+1,:).*dPhi_i(k+1,:);
           	m = m + 1;
        end
    end
else
    N_int     = [];
    dNxi_int  = [];
    dNeta_int = [];
end

% Shape function matrix
N     = [N_node; N_edge; N_int];
dNxi  = [dNxi_node; dNxi_edge; dNxi_int];
dNeta = [dNeta_node; dNeta_edge; dNeta_int];

% Delete small values
N(abs(N) < eps)         = 0;
dNxi(abs(dNxi) < eps)   = 0;
dNeta(abs(dNeta) < eps) = 0;
end