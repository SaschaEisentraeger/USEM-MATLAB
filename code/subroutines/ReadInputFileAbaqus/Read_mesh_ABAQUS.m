%% Read finite element mesh (nodes and elements) from Abaqus-Input file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 06.02.2014 (University of Magdeburg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% inp_file:     name of the file to be read
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% nodes:        [n x 2] matrix containing the x-, y-coordinates of all
%               nodes in the model; the row number corresponds to the node
%               number
% elements:     [m x 4/8] matrix containing the element topology; the row
%               number corresponds to the element number#
% element_type: Name of the element used in Abaqus
% elem_mech:    vector containing the element numbers of all elements that
%               have merely mechanical degrees of freedom
% elem_elec:    vector containing the element numbers of all elements that
%               have electric degrees of freedom
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes,elements,element_type,elem_mech,elem_elec] = ...
          Read_mesh_ABAQUS(inp_file)

%% Open File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(inp_file,'r');

% File length (eof: end of file)
fseek(fid,0,'eof');
end_position = ftell(fid);
    
% Reset seeker position (bof: beginning of file)
fseek(fid,0,'bof');
TEMP = fgetl(fid);

% Find the key word *Node
while ~contains(TEMP,'*Node')
    % Get line
    TEMP = fgetl(fid);
end

%% Read nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
index = 1;
% position of the key word *Nodes in the input file
position_nodes = ftell(fid);
TEMP = fgetl(fid);

% Find the key word *Element
while ~contains(TEMP,'*Element') &&  ftell(fid)~=end_position
    C = textscan(TEMP,'%f %f %f','delimiter', ',');
    % number of nodes in the model
    max_node_number = C{1,1};
    % read next line
    TEMP = fgetl(fid);
end

% Initialize the nodal matrix
nodes = zeros(max_node_number,2);

% Reset seeker position to nodes
fseek(fid,position_nodes,'bof');
TEMP = fgetl(fid);

% Find the key word element
while ~contains(TEMP,'*Element') && ftell(fid)~=end_position
	C = textscan(TEMP,'%f %f %f','delimiter', ',');
    % The node labels are used as indices: they must be 1, 2, 3, ...
    if C{1,1} ~= index
        fclose(fid);
        error('Read_mesh_ABAQUS:nodeLabels', ...
             ['Node labels must be numbered consecutively from 1 ' ...
              '(found label %d at position %d).'], C{1,1}, index);
    end
	nodes(index,:) = [C{1,2} C{1,3}];
    index = index + 1;
    TEMP = fgetl(fid);
end

%% Read element topology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial element type (name)
counter = -1;
% the name of element types is not longer than 10 symbols
for i = 1:10
	if(strcmp(TEMP(1,end-i+1), '=') == 0)
        counter = counter + 1;
    else
        elem_type = TEMP(1,end-counter:end);
        break
	end
end
elem_type_initial = elem_type;

% number of element types
number_of_element_types = 1;
% position of the key word *Element in the input file
position_element = ftell(fid);

% number of elements with mechanical degrees of freedom
counter_mech = 0;
% number of element with electrical degrees of freedom
counter_elec = 0;

% until one of the key words *End Instance, *End Part, *Nset is read the
% end of the element topology is not reached
while ~contains(TEMP,'*End Instance') && ...
      ~contains(TEMP,'*End Part') && ...
      ~contains(TEMP,'*Nset') && ...
      ftell(fid)~=end_position
    
    TEMP = fgetl(fid);
  
	if contains(TEMP,'*Element')
        % first element type (name)
        counter = -1;
        % the name of element types is not longer than 10 symbols
        for i = 1:10
            if(strcmp(TEMP(1,end-i+1), '=') == 0)
                counter = counter + 1;
            else
                elem_type = TEMP(1,end-counter:end);
                break
            end
        end
        number_of_element_types = number_of_element_types + 1;
        TEMP = fgetl(fid);
    else
        % finite elements with piezoelectric properties
        if contains(elem_type(end),'E')  && ...
           ~contains(TEMP,'*End Instance') && ...
           ~contains(TEMP,'*End Part')     && ...
           ~contains(TEMP,'*Nset')
       
            counter_elec = counter_elec + 1;
            
        % structural finite elements
        elseif contains(elem_type(end),'4')  && ...
               ~contains(TEMP,'*End Instance') && ...
               ~contains(TEMP,'*End Part')     && ...
               ~contains(TEMP,'*Nset')
           
            counter_mech = counter_mech + 1;
            
        elseif contains(elem_type(end),'8')  && ...
               ~contains(TEMP,'*End Instance') && ...
               ~contains(TEMP,'*End Part')     && ...
               ~contains(TEMP,'*Nset')
           
            counter_mech = counter_mech + 1;
        
        elseif contains(elem_type(end),'R')  && ...
               ~contains(TEMP,'*End Instance') && ...
               ~contains(TEMP,'*End Part')     && ...
               ~contains(TEMP,'*Nset')
           
            counter_mech = counter_mech + 1;
            
        end
	end
end

max_elem = counter_mech + counter_elec;

% initialize element type array
element_type = cell(number_of_element_types,1);
% first element type
element_type{1,1} = elem_type_initial;

% initialize topology matrix
if isempty(strfind(element_type{1},'8'))
    % linear geometry description
    topology_matrix = zeros(max_elem,4);
    flag_quadratic = 0;
else
    % quadratic geometry description
    topology_matrix = zeros(max_elem,8);
    flag_quadratic = 1;
end

% Reset seeker position to elements
fseek(fid,position_element,'bof');
TEMP = fgetl(fid);

% elements with mechanical degrees of freedom
if counter_mech > 0
    elem_mech = zeros(counter_mech,1);
end
index_mech = 1;
% elements with electric degrees of freedom
if counter_elec > 0
    elem_elec = zeros(counter_elec,1);
end
index_elec = 1;

% index of the element type (first element type)
num_elem_type = 1;

% element number
index = 1;

% read the topology
while ~contains(TEMP,'*End Instance') && ...
      ~contains(TEMP,'*End Part') && ...
      ~contains(TEMP,'*Nset') && ...
      ftell(fid)~=end_position
    
    % first element type has electric degrees of freedom
    if strcmp(element_type{num_elem_type}(end),'E')
        % quadratic geometry description
        if flag_quadratic == 1
            C = textscan(TEMP,'%d %d %d %d %d %d %d %d %d','delimiter', ',');
        	topology_matrix(index,:) = [C{1,2}  C{1,3}  C{1,4}  C{1,5}  ...
                                        C{1,6}  C{1,7}  C{1,8}  C{1,9}];
            % in case the elements of the same type are not successively
            % ordered
            elem_elec(index_elec,:) = C{1,1};
            index_elec = index_elec + 1;
        % linear geometry description
        else
            C = textscan(TEMP,'%d %d %d %d %d','delimiter', ',');
            topology_matrix(index,:) = [C{1,2} C{1,3} C{1,4} C{1,5}];
            % in case the elements of the same type are not successively
            % ordered
            elem_elec(index_elec,:) = C{1,1};
            index_elec = index_elec + 1;
        end
    else
        % first element type has mechanical degrees of freedom
        % quadratic geometry description
        if flag_quadratic == 1
            C = textscan(TEMP,'%d %d %d %d %d %d %d %d %d','delimiter', ',');
            topology_matrix(index,:) = [C{1,2}  C{1,3}  C{1,4}  C{1,5} ...
                                        C{1,6}  C{1,7}  C{1,8}  C{1,9}];
            % in case the elements of the same type are not successively
            % ordered
            elem_mech(index_mech,:) = C{1,1};
            index_mech = index_mech + 1;
       	% linear geometry description
        else
            C = textscan(TEMP,'%d %d %d %d %d','delimiter', ',');
            topology_matrix(index,:) = [C{1,2} C{1,3} C{1,4} C{1,5}];
            % in case the elements of the same type are not successively
            % ordered
            elem_mech(index_mech,:) = C{1,1};
            index_mech = index_mech + 1;
        end
    end
  
    TEMP = fgetl(fid);
    index = index + 1; % element number
  
	if ~contains(TEMP,'*Element')
    else
        num_elem_type = num_elem_type + 1;
        % first element type (name)
        counter = -1;
        % the name of element types is not longer than 10 symbols
        for i = 1:10
            if(strcmp(TEMP(1,end-i+1), '=') == 0)
                counter = counter + 1;
            else
                element_type{num_elem_type,1} = TEMP(1,end-counter:end);
                break
            end
        end
        TEMP = fgetl(fid);
	end
end

%% Reorder the element topology matrix
% elements are ordered by element type

if counter_mech ~= 0 && counter_elec ~= 0
    if strcmp(element_type{1}(end),'E')
        % the first element type has electric degrees of freedom
        elements(elem_elec,:) = topology_matrix(1:size(elem_elec,1),:);
        % the second element type has only mechanical degrees of freedom
        elements(elem_mech,:) = topology_matrix(size(elem_elec,1)+1:end,:);
    else
        % the first element type has only mechanical degrees of freedom
        elements(elem_mech,:) = topology_matrix(1:size(elem_mech,1),:);
        % the second element type has electric degrees of freedom
        elements(elem_elec,:) = topology_matrix(size(elem_mech,1)+1:end,:);
    end
else
    elements = topology_matrix;
end

if ~exist('elem_elec','var')
    elem_elec = [];
elseif ~exist('elem_mech','var')
    elem_mech = [];
end
    
fclose('all');

end