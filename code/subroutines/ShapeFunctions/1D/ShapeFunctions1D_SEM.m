%% Shape functions and their derivatives for the SEM (1D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The shape functions and shape function derivatives for the spectral
% element method (SEM) are computed. Spectral shape functions are generated
% by Lagrangian interpolation polynomials through non-equidistant nodal
% distributions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% n:        Polynomial degree (integer value)
% xi:       Local coordinate(s) of the evaluation point
% ndist:    Nodal distribution; 1: GLL, 2: CGL, 3: Equidistant
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% N:        Function value of the shape functions at xi
% dN:       Function value of the derivatives of the shape functions at xi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,dN] = ShapeFunctions1D_SEM(n,xi,ndist)

% Determine whether input is column vector or a matrix and convert to a row
% vector
if iscolumn(xi)
    xi  = xi';
elseif ismatrix(xi)
    xi = reshape(xi,1,[]);
end

% Number of evaluation points
NP = numel(xi);

% Nodal distribution - equidistant
xii = nodal_distribution(n,ndist);
[Naux,dNaux] = LagrangePolynomials(xii);

% Ordering of the nodes --> the first to shape functions correspond to the
% "corner" nodes of the 1d element for xi = -1 and xi = 1; the bubble nodes
% are ordered in positive xi-direction.
vec_order = [1, n+1, 2:(n)];

% Matrix of shape functions
N  = zeros(n+1,NP);
% Matrix of derivatives of the shape functions
dN = zeros(n+1,NP);

% Function value of the shape functions and their derivatives
for i = 1:n+1
    N(i,:)  = Naux{vec_order(i)}(xi);
    dN(i,:) = dNaux{vec_order(i)}(xi);
end

% Delete small values
N(abs(N) < eps)   = 0;
dN(abs(dN) < eps) = 0;

%-------------------------------------------------------------------------%
% [N,dN] = ShapeFunctions1D_SEM(n,(-1:0.001:1),ndist);
% % Visualization of the shape functions
% figure('Name','Shape functions')
% plot((-1:0.001:1),N)
% % Axis labels
% xlabel('x','FontName','Times','FontSize',24,'FontWeight','bold','Color','k')
% ylabel('N','FontName','Times','FontSize',24,'FontWeight','bold','Color','k')
% % Get handle to current axis
% hAxes = gca;
% set(hAxes,'FontName','Times','FontSize',24);
% % Get handle to current figure
% hFig =  gcf;
% set(hFig,'Units','normalized','Position',[0 0 1 1]);
% 
% % Visualization of the derivatives of the shape functions
% figure('Name','Derivatives')
% plot((-1:0.001:1),dN)
% % Axis labels
% xlabel('x','FontName','Times','FontSize',24,'FontWeight','bold','Color','k')
% ylabel('dN/dxi','FontName','Times','FontSize',24,'FontWeight','bold','Color','k')
% % Get handle to current axis
% hAxes = gca;
% set(hAxes,'FontName','Times','FontSize',24);
% % Get handle to current figure
% hFig =  gcf;
% set(hFig,'Units','normalized','Position',[0 0 1 1]);
%-------------------------------------------------------------------------%

end