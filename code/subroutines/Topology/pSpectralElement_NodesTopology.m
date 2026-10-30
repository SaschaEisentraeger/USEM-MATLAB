%% Adjust nodes and topology for higher order spectral elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 30.06.2017 (University of Magdeburg)
%
% New nodal coordinates and a new topology for higher order spectral
% elements is generated from the original linear/quadratic discretization. 
%
% All new nodes are placed at the images of the reference nodes under the
% geometry map of the original element (bilinear for 4-node, Serendipity
% for 8-node elements). The degree-p geometry is therefore that of the
% input mesh for every p: the elements are subparametric for p > 2 and the
% geometry description does not improve with p. Curved boundaries are only
% as accurate as in the input mesh; to obtain isoparametric elements, the
% new nodes would have to be placed on the exact geometry, which is not
% done here (see studies/thick_cylinder/matlab/cyl_model.m for a
% Gordon-Hall correction of circular boundary edges). For strongly
% distorted 8-node elements this map can send an interior node onto an edge
% node of the same element; such a mesh cannot be represented by distinct
% nodes and an error is raised (identifier
% pSpectralElement_NodesTopology:coincidentNodes).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% Topology_Node:        Nodal topology (connectivity/coincidence) matrix
% xCoor:                x-coordinates of all nodes (original mesh)
% yCoor:                y-coordinates of all nodes (original mesh)
% p:                    Polynomial degree structure containing the field 
%                       pxi and peta
% p.pxi:                Polynomial degree in xi
% p.peta:               Polynomial degree in eta
% p.ndist:              Nodal distribution
% p.nodeShift:          Optional; if > 0, single nodes of elements with a
%                       singular metric interpolation matrix (e.g.
%                       coincident nodes) are moved by p.nodeShift times
%                       the element size (ShiftSingularNodes)
% p_org:              	Polynomial degree (original mesh/discretization)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% New data for the high order spectral element:
% Topology_Node_p:      Nodal topology (connectivity/coincidence) matrix
% xCoor_p:              x-coordinates of all nodes
% yCoor_p:              y-coordinates of all nodes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Topology_Node_p,xCoor_p,yCoor_p] = ...
          pSpectralElement_NodesTopology(Topology_Node,xCoor,yCoor,p,...
          p_org)

% Polynomial degrees in xi and eta
pxi  = p.pxi;
peta = p.peta;

% Nodal distribution
ndist = p.ndist;

% Number of finite elements
NumElem = size(Topology_Node,1);

% Number of corner nodes
NumNodes = numel(unique(Topology_Node(:,1:4)));

% The corner nodes are kept with their numbers, all other nodes are
% renumbered. This requires the corner nodes to be numbered 1..NumNodes.
CornerNodes = unique(Topology_Node(:,1:4));
if ~isequal(double(CornerNodes(:)), (1:NumNodes)')
    error('pSpectralElement_NodesTopology:cornerNodes', ...
         ['The corner nodes must be numbered 1..%d (before the ' ...
          'mid-side nodes of 8-node elements).'], NumNodes);
end

% Do not change the input data if: p_xi = p_eta = p_org = 1. 
if pxi == p_org && peta == p_org && p_org == 1
    Topology_Node_p = Topology_Node;
    xCoor_p = xCoor;
    yCoor_p = yCoor;
    return
end

%-------------------------Adjust the topology data------------------------%
% Number of unique edges in the discretization
[Topology_Edge,NumEdges] = Topology_Edge_2D(Topology_Node);

% Number of new (edge) nodes (for each element 2 of 4 edges are xi-edges
% and the other 2 of 4 edges are eta-edges)
NumNewNodesEdge = NumEdges*((pxi-1) + (peta-1))/2;

% Number of interior node
nNodesInterior = (pxi-1)*(peta-1);
NumNewNodesInternal = nNodesInterior*NumElem;

% For each edge the corresponding elements are listed
[EdgeList,EdgeListLocal] = ConnectivityEdgeElement(Topology_Edge);

% Number of nodes for spectral element
NumNodesElem = 4 + 2*(pxi-1)+ 2*(peta-1) + nNodesInterior;

% Initialize the result variables
Topology_Node_p = zeros(NumElem,NumNodesElem);
xCoor_p = zeros(NumNodes+NumNewNodesEdge+NumNewNodesInternal,1);
yCoor_p = zeros(NumNodes+NumNewNodesEdge+NumNewNodesInternal,1);

% Assign the topology data corresponding to the corner nodes
Topology_Node_p(:,1:4) = Topology_Node(:,1:4);

% Assign the corner nodes
xCoor_p(1:NumNodes) = xCoor(1:NumNodes);
yCoor_p(1:NumNodes) = yCoor(1:NumNodes);

% If only linear shape functions should be used we have already assigned
% all necessary values to the result variables
if pxi == 1 && peta == 1
    return
end

% Local coordinates of the (p-1) new edge nodes
xi  = nodal_distribution(pxi,ndist);
eta = nodal_distribution(peta,ndist);
% Only interior nodes are needed
xi  = xi(2:end-1);
eta = eta(2:end-1);

% Node counter
counter = NumNodes;
% Compute the new nodes and adjust the topology data
for i = 1:NumEdges
	% Number of the element that the ith edge belongs to
	Elem = EdgeList{i}(1);
	% Indices to the nodes belonging to an edge (original input data)
	idx = [1 2 4+1+0*(p_org-1):4+1*(p_org-1); ...
		   2 3 4+1+1*(p_org-1):4+2*(p_org-1); ...
		   4 3 4+1+2*(p_org-1):4+3*(p_org-1); ...
		   1 4 4+1+3*(p_org-1):4+4*(p_org-1)];
	% Local edge number
	Edge = find(Topology_Edge(Elem,:) == i,1,'first');
	% Nodal coordinates of the ith edge
	xCoorElem = xCoor(Topology_Node(Elem,idx(Edge,:)));
	yCoorElem = yCoor(Topology_Node(Elem,idx(Edge,:)));
	% xi-edges
    if Edge == 1 || Edge == 3
        % Polynomial degree of the edge
        pe = pxi;
        % Nodal coordinates on the edge
        xie = xi;
    % eta-edges
    elseif Edge == 2 || Edge == 4
        % Polynomial degree of the edge
        pe = peta;
        % Nodal coordinates on the edge
        xie = eta;
    end
    % Compute global coordinates of the nodes
    % Linear or quadratic geometry description of the original edge of
	% the finite element
	[N,~] = ShapeFunctions1D_SEM(p_org,xie,ndist);
	xnew = N'*xCoorElem;
 	ynew = N'*yCoorElem;
	% Assign the new edge nodes
 	xCoor_p(counter+1:counter+pe-1) = xnew;
 	yCoor_p(counter+1:counter+pe-1) = ynew;

	% Adjust the topology data
	idx_EdgeNodes{1} = (5:4+(pxi-1));
    idx_EdgeNodes{2} = (5+(pxi-1):4+(pxi-1)+(peta-1));
    idx_EdgeNodes{3} = (5+(pxi-1)+(peta-1):4+2*(pxi-1)+(peta-1));
    idx_EdgeNodes{4} = (5+2*(pxi-1)+(peta-1):4+2*(pxi-1)+2*(peta-1));
	
	for j = 1:numel(EdgeList{i})
        if j == 1
            Topology_Node_p(Elem,idx_EdgeNodes{Edge}) = counter+(1:pe-1);
            Nodes = Topology_Node(Elem,idx(Edge,1:2));
        else
            % Number of the element that the ith edge belongs to
            Elem2 = EdgeList{i}(j);
            % Local edge number
            Edge2 = EdgeListLocal{i}(j);
            % Depending on the orientation of the edge the nodes need
        	% to be assigned differently
          	if sum(Nodes == Topology_Node(Elem2,idx(Edge2,1:2))) == 2
            	Topology_Node_p(Elem2,idx_EdgeNodes{Edge2}) = ...
               	counter+(1:pe-1);
            else
                Topology_Node_p(Elem2,idx_EdgeNodes{Edge2}) = ...
             	fliplr(counter+(1:pe-1));
            end
      	end
	end
	% Adjust the node counter
	counter = counter + pe - 1;
end

% Add the topology data for the interior nodes for pxi & peta >= 2
DoF_Mat = zeros(NumElem,nNodesInterior);
for i = 1:NumElem
    DoF_Mat(i,:) = (1:nNodesInterior) + counter + nNodesInterior*(i-1);
end
Topology_Node_p(:,NumNodesElem-nNodesInterior+1:NumNodesElem) = DoF_Mat;

% The positions of the interior nodes have to be adjusted
if pxi >= 2 && peta >= 2
    % Add the nodal coordinates of the interior nodes for p >= 2
    for i = 1:NumElem
        % Nodal coordinates of the ith element
        xCoorElem = xCoor(Topology_Node(i,:));
        yCoorElem = yCoor(Topology_Node(i,:));

        % Number of internal nodes per element
        % Interior nodes
        xi  = zeros(1,nNodesInterior);
        eta = zeros(1,nNodesInterior);
      	% Nodal distribution (1d)
    	Nodes_xi  = nodal_distribution(pxi,ndist);
        Nodes_eta = nodal_distribution(peta,ndist);
		% Interior nodes
     	counter_int = 0;
        for j = 2:peta
            for k = 2:pxi
                counter_int = counter_int + 1;
            	xi(counter_int)  = Nodes_xi(k);
             	eta(counter_int) = Nodes_eta(j);
            end
        end
        
        % Compute the new nodal coordinates
        % Linear or quadratic geometry description of the original finite
        % element
        if p_org == 1
          	[N,~,~] = ShapeFunctions2D_4nodes(xi,eta);
       	elseif p_org == 2
        	[N,~,~] = ShapeFunctions2D_8nodes(xi,eta);
        end
       	xnew = N'*xCoorElem;
      	ynew = N'*yCoorElem;

    	% Assign the new edge nodes
     	xCoor_p(counter+1:counter+nNodesInterior) = xnew;
      	yCoor_p(counter+1:counter+nNodesInterior) = ynew;
            
        % Adjust the node counter
        counter = counter + nNodesInterior;
    end
end

% Optional (p.nodeShift > 0): move single nodes of elements whose metric
% interpolation matrix is singular, e.g. because of coincident nodes
if isfield(p,'nodeShift') && p.nodeShift > 0
    [xCoor_p,yCoor_p] = ShiftSingularNodes(Topology_Node_p,xCoor_p,...
                        yCoor_p,p);
end

% Coincident nodes make the element (and its metric shape functions)
% degenerate; stop instead of continuing with a singular element
CheckCoincidentNodes(Topology_Node_p,xCoor_p,yCoor_p,pxi,peta);

end