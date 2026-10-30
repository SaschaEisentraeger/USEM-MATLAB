function model = usemBuildModel(nodes, topo, pdeg)
%USEMBUILDMODEL Lift a linear/quadratic mesh to spectral elements.
%
%   model = usemBuildModel(nodes, topo, pdeg) runs
%   pSpectralElement_NodesTopology and returns a struct with the fields
%       p              polynomial degree structure (pxi, peta)
%       p_org          degree of the input mesh (1: 4-node, 2: 8-node)
%       Topology_Node  connectivity of the spectral elements
%       xCoor, yCoor   nodal coordinates
%       num_Elem, num_Nodes, num_DoF

p.pxi   = pdeg;
p.peta  = pdeg;

if size(topo, 2) == 4
    p_org = 1;
else
    p_org = 2;
end

[Topology_Node, xCoor, yCoor] = pSpectralElement_NodesTopology( ...
    topo, nodes(:,1), nodes(:,2), p, p_org);

model.p             = p;
model.p_org         = p_org;
model.Topology_Node = Topology_Node;
model.xCoor         = xCoor;
model.yCoor         = yCoor;
model.num_Elem      = size(Topology_Node, 1);
model.num_Nodes     = max(Topology_Node(:));
model.num_DoF       = 2*model.num_Nodes;

end
