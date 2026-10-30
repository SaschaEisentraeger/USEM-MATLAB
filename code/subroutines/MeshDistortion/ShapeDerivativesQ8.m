function [dNdxi, dNdeta] = ShapeDerivativesQ8(xi,eta)
%SHAPEDERIVATIVESQ8 Derivatives of the Serendipity shape functions (8-node element).

% Q8 serendipity (8-node) derivatives wrt (xi,eta)
% Node order assumed:
% 1: (-1,-1), 2: (1,-1), 3: (1,1), 4: (-1,1),
% 5: (0,-1), 6: (1,0), 7: (0,1), 8: (-1,0)
% (mid-side nodes 5-8)

% Shape function derivatives (see standard serendipity Q8)
% Corner nodes
N1_xi  =  0.25*(1-eta)*(eta + 2*xi);
N1_eta =  0.25*(1-xi)*(2*eta + xi);
N2_xi  =  0.25*(1-eta)*(2*xi - eta);
N2_eta =  0.25*(1+xi)*(2*eta - xi);
N3_xi  =  0.25*(1+eta)*(eta + 2*xi);
N3_eta =  0.25*(1+xi)*(xi + 2*eta);
N4_xi  =  0.25*(1+eta)*(2*xi - eta);
N4_eta =  0.25*(1-xi)*(2*eta - xi);

% Mid-side nodes
N5_xi  =  xi*(eta - 1);
N5_eta =  0.5*(xi^2 - 1);
N6_xi  =  0.5*(1 - eta^2);
N6_eta = -eta*(1 + xi);
N7_xi  = -xi*(1 + eta);
N7_eta =  0.5*(1 - xi^2);
N8_xi  =  0.5*(eta^2 - 1);
N8_eta =  eta*(xi - 1);

dNdxi  = [N1_xi; N2_xi; N3_xi; N4_xi; N5_xi; N6_xi; N7_xi; N8_xi];
dNdeta = [N1_eta;N2_eta;N3_eta;N4_eta;N5_eta;N6_eta;N7_eta;N8_eta];
end
