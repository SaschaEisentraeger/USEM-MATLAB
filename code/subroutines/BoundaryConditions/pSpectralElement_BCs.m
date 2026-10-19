%% Adjust node sets (BCs) for higher order spectral elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 30.06.2017 (University of Magdeburg)
%
% The node sets of the boundary conditions are adjusted for higher order 
% spectral elements.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% dsload:              	Data for the distributed load (struct)
% boundary_conditions:	Data for the Dirichlet boundary conditions (struct)
% Topology_Node:        Nodal topology data (matrix)
% p:                    Polynomial degree structure containing the field 
%                       pxi and peta
% p.pxi:                Polynomial degree in xi
% p.peta:               Polynomial degree in eta
% p_org:              	Polynomial degree (original mesh/discretization)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% dsload:              	Data for the distributed load (struct) adjusted for
%                       higher order spectral elements
% boundary_conditions:	Data for the Dirichlet boundary conditions (struct)
%                       adjusted for higher order spectral elements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dsload,boundary_conditions] = pSpectralElement_BCs(...
          dsload,boundary_conditions,Topology_Node,p,p_org)

% Polynomial degrees in xi and eta
pxi = p.pxi;
peta = p.peta;
% Maximum polynomial degree
pmax = max([pxi,peta]);

% Do not change the input data if: p_xi = p_eta = p_org = 1. 
if pxi == p_org && peta == p_org && p_org == 1
    return
end

%% ------------------Adjust Dirichlet boundary conditions-----------------%
% Number of different boundary conditions that are applied to the system
number_of_ds = size(boundary_conditions,1);

% Adjust BCs for the spectral elements
for i = 1:number_of_ds
    % If a node of interest has been selected for the analysis of the
    % results this boundary condition is not adjusted
    if strcmp(boundary_conditions(i).type,'noi')
        continue
    end
        
    % Node set for the ith boundary condition
	nset = boundary_conditions(i).nset;
	if size(nset,1) == 3
        nset1 = (nset(1):nset(3):nset(2))';
        if nset1(end) == nset(2)
            nset = nset1;
        end
	end
    % Element set for the ith boundary condition
    elset = boundary_conditions(i).elset;
	if size(elset,1) == 3
        elset1 = (elset(1):elset(3):elset(2))';
        if elset1(end) == elset(2)
           	elset = elset1;
        end
	end
    
    if ~isempty(elset)
        % Matrix containing the nodes for an element edge
        pl.pxi  = 1;
        pl.peta = 1;
        [NSET] = BC_Edge_sub(pl,Topology_Node(:,1:4),nset,elset);

        % Number of edges that contain Dirichlet boundary conditions
        NumEdgesBC = size(NSET,1);

        % New constraint node set
        nset_new = zeros(NumEdgesBC*(pmax+1),1);

        for j = 1:NumEdgesBC
            % Local numbering of the nodes for the four edges of the linear 
            % element
            Index = [1 2; 2 3; 4 3; 1 4];
            % Local numbering of the nodes for the four edges of the higher
            % order spectral element
            Index_p{1} = [1 2 (5:4+(pxi-1))];
            Index_p{2} = [2 3 (5+(pxi-1):4+(pxi-1)+(peta-1))];
            Index_p{3} = [4 3 (5+(pxi-1)+(peta-1):4+2*(pxi-1)+(peta-1))];
            Index_p{4} = [1 4 (5+2*(pxi-1)+(peta-1):4+2*(pxi-1)+2*(peta-1))];

            for k = 1:4
                mat = sort([Topology_Node(elset,Index(k,1)),...
                            Topology_Node(elset,Index(k,2))],2);
                vec = sort([NSET(j,1),NSET(j,2)]);
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

            if Edge == 1 || Edge == 3
                % Polynomial degree of the edge
                pe = pxi;
            % eta-edges
            elseif Edge == 2 || Edge == 4
                % Polynomial degree of the edge
                pe = peta;
            end

            % Nodes belonging to a constraint edge
            nset_new(1+(pe+1)*(j-1):(pe+1)*(j-1)+(pe+1)) = ...
            Topology_Node(elset(h),Index_p{Edge});
        end
        % Adjust the boundary condition data
        nset_new = nset_new(nset_new ~= 0);
        boundary_conditions(i).nset = unique(nset_new);
    else
        % If only single nodes are constraint there is no need to adjust
        % the data
        boundary_conditions(i).nset = nset;
    end 
end

%% -----------------------Adjust line loads (dsload)----------------------%
% Number of different distributed loads that are applied to the system
number_of_ds = size(dsload,1);

% Adjust dsloads for the higher order finite elements
for i = 1:number_of_ds
    % Node set for the ith distributed load
	nset = dsload(i).nset;
	if size(nset,1) == 3
        nset1 = (nset(1):nset(3):nset(2))';
        if nset1(end) == nset(2)
            nset = nset1;
        end
	end
    % Element set for the ith distributed load
    elset = dsload(i).elset;
	if size(elset,1) == 3
        elset1 = (elset(1):elset(3):elset(2))';
        if elset1(end) == elset(2)
           	elset = elset1;
        end
	end
    % Matrix containing the nodes for an element edge
    pl.pxi  = 1;
    pl.peta = 1;
    [NSET] = BC_Edge_sub(pl,Topology_Node(:,1:4),nset,elset);
    
    % Number of edges that contain distributed loads
    NumEdgesBC = size(NSET,1);
    
    % New loaded node set
    nset_new = zeros(NumEdgesBC*(pmax+1),1);
    
    for j = 1:NumEdgesBC
        % Local numbering of the nodes for the four edges of the linear 
        % element
        Index = [1 2; 2 3; 4 3; 1 4];
        % Local numbering of the nodes for the four edges of the higher
        % order spectral element
        Index_p{1} = [1 2 (5:4+(pxi-1))];
        Index_p{2} = [2 3 (5+(pxi-1):4+(pxi-1)+(peta-1))];
        Index_p{3} = [4 3 (5+(pxi-1)+(peta-1):4+2*(pxi-1)+(peta-1))];
        Index_p{4} = [1 4 (5+2*(pxi-1)+(peta-1):4+2*(pxi-1)+2*(peta-1))];
       
	   for k = 1:4
	    	mat = sort([Topology_Node(elset,Index(k,1)),...
	    				Topology_Node(elset,Index(k,2))],2);
	    	vec = sort([NSET(j,1),NSET(j,2)]);
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
        
        if Edge == 1 || Edge == 3
            % Polynomial degree of the edge
            pe = pxi;
        % eta-edges
        elseif Edge == 2 || Edge == 4
            % Polynomial degree of the edge
            pe = peta;
        end
        
		% Nodes belonging to a loaded edge
        nset_new(1+(pe+1)*(j-1):(pe+1)*(j-1)+(pe+1)) = ...
        Topology_Node(elset(h),Index_p{Edge});
        
    end
    
    % Adjust the distributed load data
	nset_new = nset_new(nset_new ~= 0);
    dsload(i).nset = unique(nset_new);
end
%-------------------------------------------------------------------------%
end