%% Boundary edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Eisentraeger 25/02/2021 (UNSW Sydney)
%
% Determine all element edges that are located at the boundary of the
% model.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% Topology_Edge:	Edge topology (connectivity/coincidence) matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% Boundary:     Edge data for all element edges that are at the boundary of
%               the structure (element number, local edge number, global
%               edge number)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Boundary] = Boundary_Edges(Topology_Edge)

% Number of edges
NumEdges = max(max(Topology_Edge));

% Size of the topology matrix
sz = size(Topology_Edge);

% Multiplicity of edges
mult = histcounts(Topology_Edge,NumEdges);

% Find edges with multiplicity 1; if an edge belongs to only one element it
% has to be located on the boundary of the structure
Edge = find(mult == 1);

% Find edges in the topology data
[~,idx] = ismember(Edge,Topology_Edge);

% Convert linear index to matrix indices:
% [number of element, local number of edge]
[row,col] = ind2sub(sz,idx);

% Sort by element
[row,idx_sort] = sort(row);

% Data of the boundary edges: 
% [number of element, local number of edge, global number of edge]
Boundary = [row', col(idx_sort)', Edge(idx_sort)'];

end