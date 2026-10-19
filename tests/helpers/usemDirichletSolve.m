function [u, uExact, err] = usemDirichletSolve(K, model, field)
%USEMDIRICHLETSOLVE Solve a displacement-driven patch test.
%
%   [u, uExact, err] = usemDirichletSolve(K, model, field) prescribes the
%   exact displacement field on every node on the boundary of the
%   (rectangular) domain, sets all other loads to zero, and solves
%   K_II u_I = -K_IB u_B for the interior degrees of freedom.
%
%   field : function handle [ux, uy] = field(x, y) of an equilibrium
%           field (zero body force) for the material used in K.
%   err   : max|u - uExact| / max|uExact|

x = model.xCoor;
y = model.yCoor;
[ux, uy] = field(x, y);
uExact = zeros(model.num_DoF, 1);
uExact(1:2:end) = ux;
uExact(2:2:end) = uy;

tol = 1e-10*max(max(x) - min(x), max(y) - min(y));
bnd = abs(x - min(x)) < tol | abs(x - max(x)) < tol | ...
      abs(y - min(y)) < tol | abs(y - max(y)) < tol;
bndDoF = sort([2*find(bnd)-1; 2*find(bnd)]);
actDoF = setdiff((1:model.num_DoF)', bndDoF);

u = zeros(model.num_DoF, 1);
u(bndDoF) = uExact(bndDoF);
u(actDoF) = -K(actDoF, actDoF) \ (K(actDoF, bndDoF)*u(bndDoF));

err = max(abs(u - uExact))/max(abs(uExact));

end
