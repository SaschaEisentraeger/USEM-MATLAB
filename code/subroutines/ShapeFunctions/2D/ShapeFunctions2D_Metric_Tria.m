%% Metric shape functions and derivatives (triangular element)
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
%-------------------------------------------------------------------------%
%% Input
%
% (xQ,yQ):      Global coordinates of the query point(s)
% (xN,yN):      Global coordinates of the nodes
% p:            Polynomial degree
%-------------------------------------------------------------------------%
%% Output
%
% N:    Function value of the shape functions at (xQ,yQ)
% dNx:  Function value of the derivatives of the shape functions w.r.t. x
%       at (xQ,yQ)
% dNy:  Function value of the derivatives of the shape functions w.r.t. y
%       at (xQ,yQ)
%-------------------------------------------------------------------------%
function [N,dNx,dNy] = ShapeFunctions2D_Metric_Tria(xQ,yQ,xN,yN,p)

% Number of nodes (element)
nN = (p+1)*(p+2)/2;

% Number of query points
nQ = numel(xQ);

% Centroid of the element
xC = mean(xN);
yC = mean(yN);

% Nodal coordinates relative to centroid
xN = xN - repmat(xC,nN,1);
yN = yN - repmat(yC,nN,1);

% Query point coordinates relative to centroid
xQ = xQ - repmat(xC,nQ,1);
yQ = yQ - repmat(yC,nQ,1);

% Generalized Vandermonde matrix (transpose)
V = zeros(nN);
c = 0;
for i = 0:p
    for j = 0:p
		if (i <= 2) && (j <= 2) && (i + j <= 2)
			c = c + 1;
			V(c,:) = (xN.^j).*(yN.^i);
		end
    end
end

% Polynomial vector/matrix for the query point
P = zeros(nN,nQ);
c = 0;
for i = 0:p
    for j = 0:p
		if (i <= 2) && (j <= 2) && (i + j <= 2)
			c = c + 1;
			P(c,:) = (xQ.^j).*(yQ.^i);
		end
    end
end

% Derivative w.r.t. x
Px = zeros(nN,nQ);
c = 0;
for i = 0:p
    for j = 0:p
		if (i <= 2) && (j <= 2) && (i + j <= 2)
			c = c + 1;
			if j == 0
				Px(c,:) = zeros(nQ,1);
			else
				Px(c,:) = (j*xQ.^(j-1)).*(yQ.^i);
			end
		end
    end
end

% Derivative w.r.t. y
Py = zeros(nN,nQ);
c = 0;
for i = 0:p
    for j = 0:p
		if (i <= 2) && (j <= 2) && (i + j <= 2)
			c = c + 1;
			if i == 0
				Py(c,:) = zeros(nQ,1);
			else
				Py(c,:) = (xQ.^j).*(i*yQ.^(i-1));
			end
		end
    end
end

% Shape functions and derivatives
dV = decomposition(V);
N   = dV\P;
dNx = dV\Px;
dNy = dV\Py;

end