%% Sorts a node set according to the edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 12.11.2014 (University of Magdeburg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% p:                Polynomial degree structure containing the field pxi
%                   and peta
% p.pxi:            Polynomial degree in xi
% p.peta:           Polynomial degree in eta
% Topology_Nodes:	Nodal topology (connectivity/coincidence) matrix
% Nodes:            Node set of the nodes with boundary conditions
% Elements:         Element set of elements with boundary conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% Node_matrix: 		Matrix of nodes; each row contains the nodes for an
% 					edge with applied boundary conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Node_matrix] = BC_Edge_sub(p,Topology_Nodes,Nodes,Elements)

% Polynomial degrees in xi and eta
pxi = p.pxi;
peta = p.peta;

% Maximal polynomial degree (determines the number of interior nodes)
pmax = max([pxi,peta]);

% Local nodal numbers for the nodes on each edge
Index{1} = [1 2 (5:4+(pxi-1))];
Index{2} = [2 3 (5+(pxi-1):4+(pxi-1)+(peta-1))];
Index{3} = [4 3 (5+(pxi-1)+(peta-1):4+2*(pxi-1)+(peta-1))];
Index{4} = [1 4 (5+2*(pxi-1)+(peta-1):4+2*(pxi-1)+2*(peta-1))];
NumberOfNodesPerEdge = [pxi+1,peta+1,pxi+1,peta+1];

num_elem = length(Elements);
NodeMat = cell(num_elem,1);

for i = 1:num_elem
    % Initialize temporary variables
    Temp{1,1} = zeros(NumberOfNodesPerEdge(1),1);
    Temp{2,1} = zeros(NumberOfNodesPerEdge(2),1);
    Temp{3,1} = zeros(NumberOfNodesPerEdge(3),1);
    Temp{4,1} = zeros(NumberOfNodesPerEdge(4),1);
    
    [Nodes_vec, idx, ~] = intersect(Topology_Nodes(Elements(i),:),Nodes);
    % Only one edge of the ith element carries a load
	if length(Nodes_vec) == NumberOfNodesPerEdge(1) || ...
       length(Nodes_vec) == NumberOfNodesPerEdge(2)
        for j = 1:4
            if isempty(setdiff(idx,Index{j}))
                % Save the nodal topology data for the jth edge of the ith
                % element
                Temp{j} = Topology_Nodes(Elements(i),Index{j})';
                break
            end
        end
        NodeMat{i} = Temp;
    elseif length(Nodes_vec) > NumberOfNodesPerEdge(1) && ...
           length(Nodes_vec) > NumberOfNodesPerEdge(2)
        % Maximum number of edges of the ith element that carry a load 
        NumberOfEdges = length(Nodes_vec);
        % Loop over all edges
        Temp{1,1} = zeros(NumberOfNodesPerEdge(1),1);
        Temp{2,1} = zeros(NumberOfNodesPerEdge(2),1);
        Temp{3,1} = zeros(NumberOfNodesPerEdge(3),1);
        Temp{4,1} = zeros(NumberOfNodesPerEdge(4),1);
        for j = 1:4
        	if numel(intersect(idx,Index{j})) == NumberOfNodesPerEdge(j)
                % Save the nodal topology data for the jth edge of the ith
                % element
                Temp{j} = Topology_Nodes(Elements(i),Index{j})';
				% Count the number of edges that carry a load down to zero
                NumberOfEdges = NumberOfEdges - 1;
                if NumberOfEdges == 0
                    break
                end
        	end
        end
        NodeMat{i} = Temp;
	end
end

count = 0;
Node_matrix = zeros(4*num_elem,pmax+1);
for i = 1:num_elem
    for j = 1:4
        count = count + 1;
        Node_matrix(count,1:NumberOfNodesPerEdge(j)) = NodeMat{i}{j};
    end
end

% Delete zero rows
Idx = Node_matrix(:,1) ~= 0;
Node_matrix = Node_matrix(Idx,:);

% To delete duplicate rows we sort the variable column-wise
[Node_matrixS,~] = sort(Node_matrix,2);

% Retain only unique rows
[~,~,Idx2] = unique(Node_matrixS,'rows','stable');
Idx2 = unique(Idx2,'stable');

% Undo the sorting
Node_matrix_Unsort = Node_matrix(Idx2,:);

% Save the output variable
Node_matrix = Node_matrix_Unsort;

end