function [N,dNx,dNy,V] = ShapeFunctions2D_MetricHybrid(xQ,yQ,xN,yN,p,xiQ)
%SHAPEFUNCTIONS2D_METRICHYBRID Metric shape functions of the hybrid trial space (method 2).
%
% Hybrid trial space (method = 2)
% T_p = P_m(x,y) + span{L_a(xi) L_b(eta) : a <= px, b <= py, a + b > m},
% m = min(px,py), L_a: Legendre polynomials. The complete polynomials are
% written in the bounding-box scaled skew coordinates (affine in x, so the
% space is P_m(x,y)), the complement in the parametric coordinates, in
% which the nodes form a tensor grid. T_p contains P_m(x,y) (completeness,
% patch tests), is frame invariant and equals Q_p on affine elements. The
% complement depends on the choice of L_a, so the basis is fixed here
% (p.metricBasis is ignored).


px = p.pxi;
py = p.peta;
m  = min(px,py);

% Skew coordinates, bounding box of the nodes mapped onto [-1,1]^2
[sQ,tQ,J0] = CorotationalSkewCoordinates(xQ,yQ,xN,yN,p);
[sN,tN]    = CorotationalSkewCoordinates(xN,yN,xN,yN,p);
cs = (max(sN) + min(sN))/2;  rs = (max(sN) - min(sN))/2;
ct = (max(tN) + min(tN))/2;  rt = (max(tN) - min(tN))/2;
sN = (sN - cs)/rs;  sQ = (sQ - cs)/rs;
tN = (tN - ct)/rt;  tQ = (tQ - ct)/rt;

% Local coordinates of the nodes and of the query points
[xiN,etaN] = ReferenceNodes2D(p);
xiQq  = xiQ(:,1);
etaQq = xiQ(:,2);

% One-dimensional Legendre polynomials
[Ls,~]      = PolynomialBasis1D(sN,px,'legendre');
[Lt,~]      = PolynomialBasis1D(tN,py,'legendre');
[Lxi,~]     = PolynomialBasis1D(xiN,px,'legendre');
[Leta,~]    = PolynomialBasis1D(etaN,py,'legendre');
[Qs,dQs]    = PolynomialBasis1D(sQ,px,'legendre');
[Qt,dQt]    = PolynomialBasis1D(tQ,py,'legendre');
[Qxi,dQxi]  = PolynomialBasis1D(xiQq,px,'legendre');
[Qeta,dQeta]= PolynomialBasis1D(etaQq,py,'legendre');

% Jacobian of the geometry mapping at the query points
[~,dNxi,dNeta] = ShapeFunctions2D_SEM(xiQq,etaQq,p);
xxi  = (dNxi'*xN)';   yxi  = (dNxi'*yN)';
xeta = (dNeta'*xN)';  yeta = (dNeta'*yN)';
detJ = xxi.*yeta - xeta.*yxi;

% The functions are evaluated through (xi,eta), so the sign of det J does
% not matter (completeness holds pointwise in (xi,eta), also on folded
% elements); the gradients of the complement need J^-1 at the query points
if any(~(abs(detJ) > 0))
    error('ShapeFunctions2D_Metric:singularJacobian', ...
         ['det J = 0 at a query point: the gradients of the hybrid ' ...
          'metric shape functions (method = 2) are undefined there.']);
end

% Inverse of J0 (skew coordinates are affine: grad_x = J0^-1 grad_(s,t))
J0inv = inv(J0);

nN = (px+1)*(py+1);
nQ = numel(xQ);
V  = zeros(nN);
P  = zeros(nN,nQ);
Px = zeros(nN,nQ);
Py = zeros(nN,nQ);
c = 0;
for b = 0:py
    for a = 0:px
        c = c + 1;
        if a + b <= m
            % Complete polynomials in (x,y)
            V(c,:) = Ls(a+1,:).*Lt(b+1,:);
            P(c,:) = Qs(a+1,:).*Qt(b+1,:);
            Ps = dQs(a+1,:).*Qt(b+1,:)/rs;
            Pt = Qs(a+1,:).*dQt(b+1,:)/rt;
            Px(c,:) = J0inv(1,1)*Ps + J0inv(1,2)*Pt;
            Py(c,:) = J0inv(2,1)*Ps + J0inv(2,2)*Pt;
        else
            % Complement in the parametric coordinates
            V(c,:) = Lxi(a+1,:).*Leta(b+1,:);
            P(c,:) = Qxi(a+1,:).*Qeta(b+1,:);
            Pxi  = dQxi(a+1,:).*Qeta(b+1,:);
            Peta = Qxi(a+1,:).*dQeta(b+1,:);
            Px(c,:) = ( yeta.*Pxi - yxi.*Peta)./detJ;
            Py(c,:) = (-xeta.*Pxi + xxi.*Peta)./detJ;
        end
    end
end

CheckUnisolvent(V,px,py);

dV  = decomposition(V);
N   = dV\P;
dNx = dV\Px;
dNy = dV\Py;

end
