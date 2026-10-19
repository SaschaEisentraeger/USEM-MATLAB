%% Shape functions and their derivative for the h-version of the FEM (2D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The shape functions of a 4-node finite element and the corresponding
% derivatives with respect to the local coordinates xi and eta are
% computed.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% xi, eta:  local coordinates 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% N:        Function value of the shape functions at (xi,eta)
% dNxi:     Function value of the derivatives of the shape functions w.r.t
%           xi at (xi,eta)
% dNeta:    Function value of the derivatives of the shape functions w.r.t
%           eta at (xi,eta)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N,dNxi,dNeta] = ShapeFunctions2D_4nodes(xi,eta)

% Determine whether input is column vector or a matrix and convert to a row
% vector
if iscolumn(xi)
    xi  = xi';
    eta = eta';
elseif ismatrix(xi)
    xi  = reshape(xi,1,[]);
    eta = reshape(eta,1,[]);
end

% Shape function values
N1 = 1/4*(1-xi).*(1-eta);
N2 = 1/4*(1+xi).*(1-eta);
N3 = 1/4*(1+xi).*(1+eta);
N4 = 1/4*(1-xi).*(1+eta);
N = [N1; N2; N3; N4];

% Derivatives w.r.t xi
dN1xi = -1/4*(1-eta);
dN2xi =  1/4*(1-eta);
dN3xi =  1/4*(1+eta);
dN4xi = -1/4*(1+eta);
dNxi = [dN1xi; dN2xi; dN3xi; dN4xi];

% Derivatives w.r.t eta
dN1eta = -1/4*(1-xi);
dN2eta = -1/4*(1+xi);
dN3eta =  1/4*(1+xi);
dN4eta =  1/4*(1-xi);
dNeta = [dN1eta; dN2eta; dN3eta; dN4eta];

end