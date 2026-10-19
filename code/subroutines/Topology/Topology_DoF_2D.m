%% Degree of freedom topology (connectivity/coincidence) matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 28.04.2010 (University of Magdeburg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% Topology_Node:	Nodal topology (connectivity/coincidence) matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% Topo:             DOF topology (connectivity/coincidence) matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Topo] = Topology_DoF_2D(Topology_Node)

% Number of elements in the model
NumberOfElements = size(Topology_Node,1);

% Number of degrees of freedom per element
ndof = 2*size(Topology_Node,2);

% Initialize the topology matrix for the mechanical degrees of freedom
Topo = zeros(NumberOfElements,ndof);

% DOF connectivity matrix
Topo(:,1:2:ndof) = 2*Topology_Node-1;
Topo(:,2:2:ndof) = 2*Topology_Node;

end