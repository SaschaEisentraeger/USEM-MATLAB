function ref = cyl_reference_norms(prob, C, n)
%CYL_REFERENCE_NORMS Exact energy and L2 norm on the exact quarter annulus
%   (tensor Gauss quadrature in polar coordinates, n x n points).
if nargin < 3, n = 80; end
[r, wr] = lgwt(n, 1, 3);
[t, wt] = lgwt(n, 0, pi/2);
[Rg, Tg] = meshgrid(r, t);  W = (wt*wr.').*Rg;
x = Rg(:).*cos(Tg(:));  y = Rg(:).*sin(Tg(:));
e = [prob.exx(x, y), prob.eyy(x, y), prob.gxy(x, y)];
ref.energy = sum(sum((e*C).*e, 2).*W(:));
ref.L2 = sum((prob.ux(x, y).^2 + prob.uy(x, y).^2).*W(:));
end
