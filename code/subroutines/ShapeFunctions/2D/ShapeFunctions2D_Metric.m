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
%-------------------------------------------------------------------------%
%% Output
%
% N:    Function value of the shape functions at (xQ,yQ)
% dNx:  Function value of the derivatives of the shape functions w.r.t. x
%       at (xQ,yQ)
% dNy:  Function value of the derivatives of the shape functions w.r.t. y
%       at (xQ,yQ)
%-------------------------------------------------------------------------%
function [N,dNx,dNy] = ShapeFunctions2D_Metric(xQ,yQ,xN,yN,p,method)

if nargin == 5
    % Global (Cartesian) coordinates: 0
    % Corotational skew coordinates:  1
    method = 1;
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

% Generalized Vandermonde matrix (transpose)
V = zeros(nN);
c = 0;
for i = 0:py
    for j = 0:px
        c = c + 1;
        V(c,:) = (xN.^j).*(yN.^i);
    end
end

% Polynomial vector/matrix for the query point
P = zeros(nN,nQ);
c = 0;
for i = 0:py
    for j = 0:px
        c = c + 1;
        P(c,:) = (xQ.^j).*(yQ.^i);
    end
end

% Derivative w.r.t. x
Px = zeros(nN,nQ);
c = 0;
for i = 0:py
    for j = 0:px
        c = c + 1;
        if j == 0
            Px(c,:) = zeros(nQ,1);
        else
            Px(c,:) = (j*xQ.^(j-1)).*(yQ.^i);
        end
    end
end

% Derivative w.r.t. y
Py = zeros(nN,nQ);
c = 0;
for i = 0:py
    for j = 0:px
        c = c + 1;
        if i == 0
            Py(c,:) = zeros(nQ,1);
        else
            Py(c,:) = (xQ.^j).*(i*yQ.^(i-1));
        end
    end
end

% The metric shape functions are undefined if the nodes are not unisolvent
% for the polynomial space (e.g. coincident nodes or strongly distorted
% elements): stop instead of returning meaningless shape functions
rc = rcond(V);
if ~(rc >= eps)
    error('ShapeFunctions2D_Metric:singular', ...
         ['The interpolation matrix of the metric shape functions is ' ...
          'singular to working precision (rcond = %.2e): the element ' ...
          'nodes are (nearly) not unisolvent. The element is too ' ...
          'strongly distorted for the metric shape functions of degree ' ...
          '(%d,%d).'], rc, px, py);
end

% Shape functions
dV = decomposition(V);
N  = dV\P;

if method == 1
    % Derivatives w.r.t. to the skew corotational coordinates
    dNxi = dV\Px;
    dNeta = dV\Py;
    % Derivatives w.r.t. to the global coordinates
    dNx = zeros(size(dNxi));
    dNy = zeros(size(dNeta));
    for n = 1:nN
        gradN = J\[dNxi(n,:);dNeta(n,:)];
        dNx(n,:) = gradN(1,:);
        dNy(n,:) = gradN(2,:);
    end
else
    % Derivatives w.r.t. to the global coordinates
    dNx = dV\Px;
    dNy = dV\Py;
end

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