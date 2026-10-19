%% Shape functions and their derivative for the h-version of the FEM (2D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The shape functions of a 8-node finite element and the corresponding
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
function [N,dNxi,dNeta] = ShapeFunctions2D_8nodes(xi,eta)

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
N1 = -(1-xi).*(1-eta).*(1+xi+eta)/4;
N2 = -(1+xi).*(1-eta).*(1-xi+eta)/4;
N3 = -(1+xi).*(1+eta).*(1-xi-eta)/4;
N4 = -(1-xi).*(1+eta).*(1+xi-eta)/4;
N5 =  (1-xi.*xi).*(1-eta)/2;
N6 =  (1+xi).*(1-eta.*eta)/2;
N7 =  (1-xi.*xi).*(1+eta)/2;
N8 =  (1-xi).*(1-eta.*eta)/2;
N = [N1; N2; N3; N4; N5; N6; N7; N8];

% Derivatives w.r.t xi
dN1xi = -(-(1-eta).*(1+xi+eta)+(1-xi).*(1-eta))/4;
dN2xi = -( (1-eta).*(1-xi+eta)-(1+xi).*(1-eta))/4;
dN3xi = -( (1+eta).*(1-xi-eta)-(1+xi).*(1+eta))/4;
dN4xi = -(-(1+eta).*(1+xi-eta)+(1-xi).*(1+eta))/4;
dN5xi = -xi.*(1-eta);
dN6xi = (1-eta.*eta)/2;
dN7xi = -xi.*(1+eta);
dN8xi = -(1-eta.*eta)/2;
dNxi = [dN1xi; dN2xi; dN3xi; dN4xi; dN5xi; dN6xi; dN7xi; dN8xi];

% Derivatives w.r.t eta
dN1eta = -(-(1-xi).*(1+xi+eta)+(1-xi).*(1-eta))/4;
dN2eta = -(-(1+xi).*(1-xi+eta)+(1+xi).*(1-eta))/4;
dN3eta = -( (1+xi).*(1-xi-eta)-(1+xi).*(1+eta))/4;
dN4eta = -( (1-xi).*(1+xi-eta)-(1-xi).*(1+eta))/4;
dN5eta = -(1-xi.*xi)/2;
dN6eta = -eta.*(1+xi);
dN7eta =  (1-xi.*xi)/2;
dN8eta = -eta.*(1-xi);
dNeta = [dN1eta; dN2eta; dN3eta; dN4eta; dN5eta; dN6eta; dN7eta; dN8eta];

end