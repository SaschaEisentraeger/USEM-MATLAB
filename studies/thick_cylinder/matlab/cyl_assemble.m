function K = cyl_assemble(model, C, form, nq)
%CYL_ASSEMBLE Global stiffness matrix.
%   form = 'PM' : unsymmetric spectral element of the code base (Ke_2D,
%                 vectorized in cyl_Ke_PM:
%                 parametric test / metric trial functions, |J| cancelled)
%          'PP' : conventional symmetric spectral element (cyl_Ke_PP)
%   nq   : Gauss points per direction (default 2p, as in Ke_2D)
p = model.p;
if nargin < 4 || isempty(nq), nq = 2*p.pxi; end
T = model.Topology_Node;
[nEl, nen] = size(T);
nd = 2*nen;
I = zeros(nd*nd, nEl);  J = I;  V = I;
for e = 1:nEl
    ids = T(e, :);
    x = model.xCoor(ids);  y = model.yCoor(ids);
    if strcmp(form, 'PM')
        Ke = cyl_Ke_PM(p, C, x, y, nq);
    else
        Ke = cyl_Ke_PP(p, C, x, y, nq);
    end
    dof = reshape([2*ids - 1; 2*ids], [], 1);
    [Ie, Je] = ndgrid(dof, dof);
    I(:, e) = Ie(:);  J(:, e) = Je(:);  V(:, e) = Ke(:);
end
K = sparse(I(:), J(:), V(:), model.num_DoF, model.num_DoF);
end
