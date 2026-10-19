%% Connectivity/topology data of the edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 18.12.2016 (University of Magdeburg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% Topology_Node:	Nodal topology (connectivity/coincidence) matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% Topology_Edge:  	Edge topology (connectivity/coincidence) matrix
% NumEdges:        	Total number of edges in the FE model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Topology_Edge,NumEdges] = Topology_Edge_2D(Topology_Node)

if size(Topology_Node,1) == 1
    EdgeNumbers = [1;2;3;4];
else
% Nodal numbers that make up the local edges of an element. The first edge
% consists of nodes 1 and 2. The second edge consists of nodes 2 and 3. The
% third edge consists of nodes 3 and 4. The fourth and last edge consist of 
% nodes 4 and 1.
idx = [1 2; 2 3; 3 4; 4 1];

% Matrix containing the edges for each element.
% ordering: edges = [Edge 1 of element 1; 
%                    Edge 1 of element 2;
%                    ...;
%                    Edge 1 of element n;
%                    Edge 2 of element 1; 
%                    Edge 2 of element 2;
%                    ...; 
%                    Edge 2 of element n;
%                    Edge 3 of element 1; 
%                    Edge 3 of element 2;
%                    ...; 
%                    Edge 4 of element n;
edges = sort(cat(1,Topology_Node(:,idx(1,:)), Topology_Node(:,idx(2,:)),...
                   Topology_Node(:,idx(3,:)), Topology_Node(:,idx(4,:))),2);

% Determine the unique edges and their global numbers
[~,~,EdgeNumbers] = unique(edges,'rows');
end

% Number of edges
NumEdges = max(EdgeNumbers);

% Edge topology - global edge numbers for each element
Topology_Edge = reshape(EdgeNumbers,[],4);

end