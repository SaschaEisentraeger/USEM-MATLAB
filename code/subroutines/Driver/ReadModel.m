function model = ReadModel(filename,p,verbose)
%READMODEL Read an ABAQUS input file and lift the mesh to spectral elements.
%
%   model = ReadModel(filename, p) reads mesh, materials, step, loads and
%   boundary conditions from the ABAQUS input file and creates the nodes of
%   the spectral elements of degree p.pxi x p.peta (p: polynomial degree
%   structure of DefaultSettings). filename is a full path or the name of
%   a file on the MATLAB path.
%
%   model = ReadModel(filename, p, verbose) prints progress information if
%   verbose is true (default).
%
%   Fields of model:
%     filename, name             input file (full path) and its base name
%     nodes, topology            mesh of the input file (4- or 8-node)
%     element_type               ABAQUS element type
%     material_names, material_elset, material_properties
%     step_name, step_type, t_end, numSteps
%     cload, dsload, body_force, gravity, boundary_conditions
%                                loads and boundary conditions; dsload and
%                                boundary_conditions include the
%                                high-order edge nodes
%     p                          polynomial degree structure (pxi = peta
%                                if the element axes are not aligned)
%     p_org                      degree of the input mesh (1: 4-node,
%                                2: 8-node)
%     Topology_Node              connectivity of the spectral elements
%     xCoor, yCoor               nodal coordinates
%     num_Elem, num_Nodes, num_DoF

if nargin < 3
    verbose = true;
end

model.filename = ResolveInputFile(filename);
[~,model.name] = fileparts(model.filename);

% Since ABAQUS 6.12, element and node set names can be used several times
% in different parts of the input file, while the readers rely on unique
% names. The sets are therefore renamed in a temporary copy of the file.
adjusted = [tempname '.inp'];
Adjust_Set_names_input_file_ABAQUS(model.filename,adjusted);
cleanup = onCleanup(@() delete(adjusted));

PrintStep(verbose,'Read mesh')
[model.nodes,model.topology,model.element_type] = Read_mesh_ABAQUS(adjusted);
PrintStep(verbose,'Read material data')
[model.material_names,model.material_elset,model.material_properties] = ...
    Read_material_ABAQUS(adjusted);
PrintStep(verbose,'Read step definition')
[model.step_name,model.step_type,model.t_end,model.numSteps] = ...
    Read_step_ABAQUS(adjusted);
PrintStep(verbose,'Read loads')
[model.cload,model.dsload,model.body_force,model.gravity] = ...
    Read_load_ABAQUS(adjusted);
PrintStep(verbose,'Read boundary conditions')
model.boundary_conditions = Read_boundary_condition_ABAQUS(adjusted);

% Degree of the input mesh
if size(model.topology,2) == 4
    model.p_org = 1;
else
    model.p_org = 2;
end

% Elements whose local axes are not aligned need pxi = peta
p = AlignPolynomialDegrees(model.topology,p,verbose);
model.p = p;

% Edge and interior nodes of the spectral elements (tensor-product space);
% the node sets of loads and boundary conditions get the new edge nodes
PrintStep(verbose,'Create the nodes of the spectral elements')
[model.Topology_Node,model.xCoor,model.yCoor] = ...
    pSpectralElement_NodesTopology(model.topology,model.nodes(:,1),...
    model.nodes(:,2),p,model.p_org);
[model.dsload,model.boundary_conditions] = pSpectralElement_BCs(...
    model.dsload,model.boundary_conditions,model.Topology_Node,p,model.p_org);

model.num_Elem  = size(model.Topology_Node,1);
model.num_Nodes = max(model.Topology_Node(:));
model.num_DoF   = 2*model.num_Nodes;

end
