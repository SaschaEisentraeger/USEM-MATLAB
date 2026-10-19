function S = usemReadAbaqus(inpFile, p)
%USEMREADABAQUS Read an ABAQUS input file like Start_Sim.m does.
%
%   S = usemReadAbaqus(inpFile, p) runs the input part of Start_Sim.m
%   (set renaming, mesh, material, step, loads, boundary conditions,
%   lifting to spectral elements) without the parallel pool, plots and
%   clearvars of the script. p is the polynomial degree structure
%   (pxi, peta, ndist). The adjusted copy of the input file is written to
%   and deleted from the current folder, so call this function from a
%   temporary folder.

[~, name, ext] = fileparts(inpFile);
local = [name ext];
if ~isfile(local)
    copyfile(inpFile, local);
end
Adjust_Set_names_input_file_ABAQUS(local);
adj = [name '_adjusted' ext];
cleanup = onCleanup(@() delete(adj));

[S.nodes, S.topology] = Read_mesh_ABAQUS(adj);
[S.material_names, S.material_elset, S.material_properties] = Read_material_ABAQUS(adj);
[S.step_name, S.step_type, S.t_end] = Read_step_ABAQUS(adj);
[S.cload, S.dsload, S.body_force, S.gravity] = Read_load_ABAQUS(adj);
S.boundary_conditions = Read_boundary_condition_ABAQUS(adj);

S.p = p;
S.p_org = 1 + (size(S.topology, 2) == 8);
[S.Topology_Node, S.xCoor, S.yCoor] = pSpectralElement_NodesTopology( ...
    S.topology, S.nodes(:,1), S.nodes(:,2), p, S.p_org);
[S.dsload, S.boundary_conditions] = pSpectralElement_BCs(S.dsload, ...
    S.boundary_conditions, S.Topology_Node, p, S.p_org);

S.num_Elem  = size(S.Topology_Node, 1);
S.num_Nodes = max(S.Topology_Node(:));
S.num_DoF   = 2*S.num_Nodes;

end
