function prob = rp_exact(kind, E, nu)
%RP_EXACT Exact solutions of the panel [0,2] x [-1,1] (Rajendran 2010,
%   Sec. 5.5), supports u_x(0,+-1) = 0, u(0,0) = 0.
%   kind = 'tension' : sigma_x = 1000
%          'moment'  : sigma_x = -1500 y      (net moment 1000)
%          'shear'   : sigma_x = -1500 (2-x) y, tau_xy = 750 (1-y^2)
%                      (net shear force 1000 at x = 2)
%   prob.ux, uy, exx, eyy, gxy, sxx, syy, sxy : function handles @(x,y)
G = E/(2*(1 + nu));
z = @(x, y) zeros(size(x));
switch kind
    case 'tension'
        sx = @(x, y) 1000 + z(x, y);  txy = z;
        prob.ux = @(x, y) 1000*x/E;
        prob.uy = @(x, y) -1000*nu*y/E;
    case 'moment'
        sx = @(x, y) -1500*y + z(x, y);  txy = z;
        prob.ux = @(x, y) -1500*x.*y/E;
        prob.uy = @(x, y) 750*(nu*y.^2 + x.^2)/E;
    case 'shear'
        sx = @(x, y) -1500*(2 - x).*y;
        txy = @(x, y) 750*(1 - y.^2) + z(x, y);
        prob.ux = @(x, y) 250*y.*(-nu*y.^2 + nu + 3*x.^2 - 12*x - 2*y.^2 + 2)/E;
        prob.uy = @(x, y) -250*(3*nu*x.*y.^2 - 5*nu*x - 6*nu*y.^2 + x.^3 - 6*x.^2 - 4*x)/E;
end
prob.kind = kind;
prob.sxx = sx;  prob.syy = z;  prob.sxy = txy;
prob.exx = @(x, y) sx(x, y)/E;
prob.eyy = @(x, y) -nu*sx(x, y)/E;
prob.gxy = @(x, y) txy(x, y)/G;
end
