function K = usemAssemble(model, C)
%USEMASSEMBLE Global stiffness matrix of the unsymmetric element (Ke_2D).
%
%   K = usemAssemble(model, C)
%   model : struct from usemBuildModel
%   C     : 3x3 plane elasticity matrix

n = model.num_DoF;
K = sparse(n, n);
T = Topology_DoF_2D(model.Topology_Node);
for e = 1:model.num_Elem
    nd = model.Topology_Node(e,:);
    Ke = Ke_2D(model.p, C, model.xCoor(nd), model.yCoor(nd));
    K(T(e,:), T(e,:)) = K(T(e,:), T(e,:)) + Ke; %#ok<SPRIX>
end

end
