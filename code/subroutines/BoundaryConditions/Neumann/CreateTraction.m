%% Two-dimensional distributed load - Surface traction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 20.01.2011 (University of Magdeburg)
%
% Compute the (energetically equivalent) nodal forces for an element edge
% loaded by a surface traction.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% p:            	Polynomial degree structure
% p.pxi:        	Polynomial degree in xi
% p.peta:        	Polynomial degree in eta
% p.ndist:      	Nodal distribution
% x:                x-coordinates of the nodes
% y:                y-coordinates of the nodes
% num_DoF:          Number of degrees of freedom
% Topology_Node: 	Nodal topology (connectivity/coincidence) matrix
% Nodes:            Nodal numbers of the corner nodes of the edges with the
%                   applied traction load
% amp:            	Amplitude of the load
% DirVec:           Direction vector of the load
% elset:            Element set of the load
% flag_uni:         Uniform distribution of the surface traction (1) or a
%                   spatially varying on (0)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% F_mech:          	External load vector due to a surface traction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F_mech] = CreateTraction(p,x,y,num_DoF,Topology_Nodes,Nodes,...
                    amp,DirVec,elset,flag_uni)

% Polynomial degrees in xi and eta
pxi = p.pxi;
peta = p.peta;

% Local nodal numbers for the nodes on each edge
idxEdge{1} = [1 2 (5:4+(pxi-1))];
idxEdge{2} = [2 3 (5+(pxi-1):4+(pxi-1)+(peta-1))];
idxEdge{3} = [4 3 (5+(pxi-1)+(peta-1):4+2*(pxi-1)+(peta-1))];
idxEdge{4} = [1 4 (5+2*(pxi-1)+(peta-1):4+2*(pxi-1)+2*(peta-1))];
NumberOfNodesPerEdge = [pxi+1,peta+1,pxi+1,peta+1];

% Number of faces with applied surface traction load
loop = size(Nodes,1);

% Initialize the external load vector
F_mech = zeros(num_DoF,1);

% Elemental force vector
Fe = cell(loop,1);

% Elemental degrees of freedom
DoFe = cell(loop,1);

% Nodal coordinates - element
xCell = cell(loop,1);
yCell = cell(loop,1);

% Nodal coordinates - edge
xCellEdge = cell(loop,1);
yCellEdge = cell(loop,1);

% Local edge number
Edge_vec = zeros(loop,1);

% Polynomial degree distribution for each edge
p_edge = cell(loop,1);

% Integration points for each edge
xi_edge   = cell(loop,1);
wxi_edge  = cell(loop,1);

% Initialize traction vector for each edge
tract_cell = cell(loop,1);

% Data for each edge with a traction load
for i = 1:loop
    % (Local) nodal numbers of the corner nodes of the element edges
    Index = [1 2; 2 3; 4 3; 1 4];
	
	for k = 1:4
		mat = sort([Topology_Nodes(elset,Index(k,1)),...
					Topology_Nodes(elset,Index(k,2))],2);
		vec = sort([Nodes(i,1),Nodes(i,2)]);
		idx_logical = ismember(mat,vec,'rows');
        idx = find(idx_logical == 1,1,'first');
        if ~isempty(idx)
            break
        end
	end
	
	% Index to the element in the element set elset
    h = idx;
    
    % Local number of the edge
    Edge = k;
    % Save local edge number
    Edge_vec(i) = Edge;
    
    % Nodal coordinates - element
    xCell{i} = x(Topology_Nodes(elset(h),:));
    yCell{i} = y(Topology_Nodes(elset(h),:));
    
    % Nodal coordinates - edge
    xCellEdge{i} = x(Topology_Nodes(elset(h),idxEdge{Edge}));
    yCellEdge{i} = y(Topology_Nodes(elset(h),idxEdge{Edge}));
    
    % Integration points
    % Matrix of shape functions
	if Edge == 1 || Edge == 3
        pe.p = pxi;
	elseif Edge == 2 || Edge == 4
     	pe.p = peta;
	end
    [xi,wxi]   = lgwt(pe.p+1,-1,1);
    % nodal distribution
    pe.ndist = p.ndist;
    
    % Save data in cell array
    p_edge{i}    = pe;
    xi_edge{i}   = xi;
    wxi_edge{i}  = wxi;
    
    % Press vector (initialize) - Element load vector
	tract_cell{i} = zeros(2*NumberOfNodesPerEdge(Edge),1);
    
    % Assign element load vector to the (global) external load vector
    vec_dof = zeros(1,2*NumberOfNodesPerEdge(Edge));
    % Find possible zero components
    idx = find(Nodes(i,:)==0,1,'first');
    if isempty(idx)
        idx = numel(Nodes(i,:)) + 1;
    end
    vec_dof(1:2:end) = 2*Nodes(i,1:idx-1)-1;
    vec_dof(2:2:end) = 2*Nodes(i,1:idx-1);
    DoFe{i} = vec_dof;
end

% Spatial loading function
% Amplitude
Param(1) = amp;
% Identifier for traction load
Param(2) = 2;

% Determine the element force vectors due to the applied Neumann boundary
% conditions
parfor counter = 1:loop
	% Auxiliary variable
	PL = [];
	
    % Nodal coordinates
    xCoor = xCell{counter};
    yCoor = yCell{counter};
    
    % Integration points
    xi  = xi_edge{counter};
    wxi = wxi_edge{counter};
    
    % Number of integration points
    NumIP = numel(xi);

    % Element load vector
	tract_elem = tract_cell{counter};
    
    % Polynomial degree
    pe = p_edge{counter};
    
    % Identity matrix
    I = eye(2);

    % Shape functions
    Naux = ShapeFunctions1D_SEM(pe.p,xi,pe.ndist);
    
    % Local edge number
    Edge = Edge_vec(counter);
    
    % Compute element load vector
    for i = 1:NumIP
        % Matrix of shape functions
        N = kron(Naux(:,i)',I);

        % Determine the unit normal vectors of the local coordinate system 
        % at point P
        xE = xCellEdge{counter};
        yE = yCellEdge{counter};
        [~,V2,~] = LocalCoordinateSystem(xi(i),xE,yE,pe.ndist);
        
        % Determine the traction value at the integration point
        if flag_uni == 1
            amp_val = amp;
        else
            % Local coordinates of the integration point
            if Edge == 1
                PL = [xi(i), -1];
            elseif Edge == 2
                PL = [1, xi(i)];
            elseif Edge == 3
                PL = [xi(i), 1];
            elseif Edge == 4
                PL = [-1, xi(i)];
            end
            % Compute global coordinates of the integration point
            [CoorPG] = Local2GlobalCoordinatesQuad(p,PL,[xCoor,yCoor]);
            % Amplitude
            amp_val = SpatialLoadingFunction(Param,CoorPG);
        end
        
        % Mapping
        map     = norm(V2);
        amp_vec = DirVec*amp_val;
                   
        % Element load vector
        tract_elem(:,1) = tract_elem(:,1) + N'*map*amp_vec*wxi(i);
    end

    % Save the load
    Fe{counter} = tract_elem;

end

% Assemble the system force vector
for i = 1:loop
    F_mech(DoFe{i}) = F_mech(DoFe{i}) + Fe{i};
end

end