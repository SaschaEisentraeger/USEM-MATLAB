function K = rp_assemble(model, C, nq, form, method)
%RP_ASSEMBLE Global stiffness matrix (see rp_Ke).
T = model.Topology_Node;
[nEl, nen] = size(T);
nd = 2*nen;
I = zeros(nd*nd, nEl);  J = I;  V = I;
for e = 1:nEl
    ids = T(e, :);
    Ke = rp_Ke(model.p, C, model.xCoor(ids), model.yCoor(ids), nq, form, method);
    dof = reshape([2*ids - 1; 2*ids], [], 1);
    [Ie, Je] = ndgrid(dof, dof);
    I(:, e) = Ie(:);  J(:, e) = Je(:);  V(:, e) = Ke(:);
end
K = sparse(I(:), J(:), V(:), model.num_DoF, model.num_DoF);
end
