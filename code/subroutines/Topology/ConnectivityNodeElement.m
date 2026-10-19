%% For each node all elements are listed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 09.06.2011 (University of Magdeburg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% Topology_Node: 	Nodal topology (connectivity/coincidence) matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% Res:              Cell array containing all finite elements that share
% 					a node; the nodal number corresponds to the cell index
% NumNodes:         Number of nodes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Res,NumNodes] = ConnectivityNodeElement(Topology_Node)

% Determine the number of nodes in the finite element model
NumNodes = numel(unique(Topology_Node));

% Cell array containing the element numbers sharing a specific node
Res = cell(NumNodes,1);

% Determine which elements share a specific node
parfor i = 1:NumNodes
	[row,~] = find(Topology_Node == i);
	Res{i} = unique(row);
end

end