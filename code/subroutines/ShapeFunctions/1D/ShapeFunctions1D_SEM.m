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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% N:        Function value of the shape functions at xi
% dN:       Function value of the derivatives of the shape functions at xi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,dN] = ShapeFunctions1D_SEM(n,xi)

% Convert the input to a row vector
xi = reshape(xi,1,[]);

% Number of evaluation points
NP = numel(xi);

% Nodal distribution - equidistant
xii = nodal_distribution(n);
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

end