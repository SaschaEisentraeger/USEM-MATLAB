function [u, act] = rp_solve(model, K, F)
%RP_SOLVE Static solution with the three point supports of Fig. 8.
act = setdiff((1:model.num_DoF).', model.fixed);
u = zeros(model.num_DoF, 1);
u(act) = K(act, act)\F(act);
end
