function [dNdxi, dNdeta] = ShapeDerivativesQ4(xi,eta)
%SHAPEDERIVATIVESQ4 Derivatives of the bilinear shape functions (4-node element).

% Q4 bilinear derivatives wrt (xi,eta) at the given point
% Node order assumed: [(-1,-1),(+1,-1),(+1,+1),(-1,+1)]
dNdxi = 0.25*[-(1-eta); +(1-eta); +(1+eta); -(1+eta)];
dNdeta= 0.25*[-(1-xi);  -(1+xi); +(1+xi); +(1-xi)];
end
