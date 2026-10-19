%% Boundary_conditions
% Boundary_conditions are added to the FE-system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 12.02.2014 (University of Magdeburg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% boundary_conditions:		Data structure containing the boundary 
% 							condition data (read from Abaqus *.inp file)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% DOF_C:	Vector of fixed degrees of freedom
% NOI:		Nodes of interest (results for transient simulation are only
%           saved at these nodes)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DOF_C,NOI] = DirichletBCs_homogeneous(boundary_conditions)

% Number of different boundary conditions that are applied to the system
number_of_BCs = size(boundary_conditions,1);

% Degrees of freedom
DOF_C = cell(number_of_BCs,1);
NOI = cell(number_of_BCs,1); 

% Implement BCs
for i = 1:number_of_BCs
    % All mechanical degrees of freedom are fixed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(boundary_conditions(i).type,'fixed') || ...
       strcmp(boundary_conditions(i).type,'UxUy 0')
        nset = boundary_conditions(i).nset;
        if size(nset,1) == 3
            nset1 = (nset(1):nset(3):nset(2))';
            if nset1(end) == nset(2)
                nset = nset1;
            end
        end

        [DOF] = CreateBoundaryConditions_Node(nset,'xy');
        
    % Symmetry boundary conditions in x-direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(boundary_conditions(i).type,'XSYMM') || ...
           strcmp(boundary_conditions(i).type,'Ux 0')
        nset = boundary_conditions(i).nset;
        if size(nset,1) == 3
            nset1 = (nset(1):nset(3):nset(2))';
            if nset1(end) == nset(2)
                nset = nset1;
            end
        end
        
        [DOF] = CreateBoundaryConditions_Node(nset,'x');
        
    % Symmetry boundary conditions in y-direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif strcmp(boundary_conditions(i).type,'YSYMM') || ...
           strcmp(boundary_conditions(i).type,'Uy 0')
        nset = boundary_conditions(i).nset;
        if size(nset,1) == 3
            nset1 = (nset(1):nset(3):nset(2))';
            if nset1(end) == nset(2)
                nset = nset1;
            end
        end
       
    [DOF] = CreateBoundaryConditions_Node(nset,'y');

    % Node of interest
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % boundary condition for the rotation about the z-axis means a history 
    % output for this node set is requested
    elseif strcmp(boundary_conditions(i).type,'noi')
        % Nodes of interest
         noi = boundary_conditions(i).nset;
         noi = double(noi);
    end
    
    % Save fixed (constraint) degrees of freedom
    if exist('DOF','var')
        DOF_C{i} = double(DOF);
        % Clear variables
        clear DOF
    else
        DOF_C{i} = double.empty(0,1);
    end
    if exist('noi','var')
        NOI{i} = double(noi);
        % Clear variables
        clear noi
    else
        NOI{i} = double.empty(0,1);
    end
    
end

% Convert cell array into matrix
DOF_C = cell2mat(DOF_C);
NOI = cell2mat(NOI);

% Ordering of the nodes of interest
NOI = fliplr(NOI);

% Delete duplicate entries
DOF_C = unique(DOF_C);
NOI = unique(NOI,'stable');

end