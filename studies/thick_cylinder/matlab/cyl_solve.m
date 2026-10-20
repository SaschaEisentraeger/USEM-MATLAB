function [u, uEx, act] = cyl_solve(model, K, F, prob)
%CYL_SOLVE Static solution with Dirichlet conditions on the straight edges.
%   'lame'        : symmetry conditions u_x = 0 (x = 0), u_y = 0 (y = 0)
%   'poly'/'smooth': both components prescribed (exact nodal values) on
%                    x = 0 and y = 0
x = model.xCoor;  y = model.yCoor;
uEx = zeros(model.num_DoF, 1);
uEx(1:2:end) = prob.ux(x, y);
uEx(2:2:end) = prob.uy(x, y);
if strcmp(prob.kind, 'lame')     % symmetry conditions
    fix = [2*model.set.xsym - 1; 2*model.set.ysym];
else
    b = union(model.set.xsym, model.set.ysym);
    fix = [2*b - 1; 2*b];
end
fix = unique(fix);
act = setdiff((1:model.num_DoF).', fix);
u = zeros(model.num_DoF, 1);
u(fix) = uEx(fix);
u(act) = K(act, act) \ (F(act) - K(act, fix)*u(fix));
end
