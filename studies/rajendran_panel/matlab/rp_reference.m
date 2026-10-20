function ref = rp_reference(prob, C)
%RP_REFERENCE Exact energy and L2 norm on the panel [0,2] x [-1,1].
[x, wx] = lgwt(20, 0, 2);
[y, wy] = lgwt(20, -1, 1);
[X, Y] = meshgrid(x, y);  W = wy*wx.';
e = [prob.exx(X(:), Y(:)), prob.eyy(X(:), Y(:)), prob.gxy(X(:), Y(:))];
ref.energy = sum(sum((e*C).*e, 2).*W(:));
ref.L2 = sum((prob.ux(X(:), Y(:)).^2 + prob.uy(X(:), Y(:)).^2).*W(:));
end
