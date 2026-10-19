%% Edge connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 30.06.2017 (University of Magdeburg)
%
% For each edge in the model the (two) elements that share that edge are
% listed. If only one element is connected to a certain edge this edge has
% to be on the boundary of the structure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% Topology_Edge: 	Edge topology (connectivity/coincidence) matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% Res:              Matrix of element numbers (each row contains the two
% 					elements that share the edge with the global edge
% 					number corresponding to the row number)
% Res2:             Local edge numbers
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Res,Res2] = ConnectivityEdgeElement(Topology_Edge)

% Determine the number of edges in the finite element model
NumEdges = numel(unique(Topology_Edge));

% Cell array containing the element numbers sharing a specific edge
Res = cell(NumEdges,1);

% Local edge number
Res2 = cell(NumEdges,1);

% Determine which elements share a specific edge
parfor i = 1:NumEdges
	[row,col] = find(Topology_Edge == i);
	[Res{i},idx] = unique(row);
    Res2{i} = col(idx);
end

end