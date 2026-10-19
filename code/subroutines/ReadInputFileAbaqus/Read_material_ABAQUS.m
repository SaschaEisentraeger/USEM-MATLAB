%% Read material properties from Abaqus-Input file
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
% material_names:       Names of the deployed material
% material_elset:       Element sets belonging to the different materials
% material_properties:  Material properties of the different materials
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [material_names, material_elset, material_properties] = ...
          Read_material_ABAQUS(inp_file)

%% Open File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(inp_file,'r');

% File length (eof: end of file)
fseek(fid,0,'eof');
end_position = ftell(fid);
    
% Reset seeker position (bof: beginning of file)
fseek(fid,0,'bof');
TEMP = fgetl(fid);

% number of different materials in the model
number_of_materials = 0;

% Find the key word *End Instance
while ~contains(TEMP,'*End Instance') && ...
      ftell(fid)~=end_position
    if number_of_materials == 0 && ...
       contains(TEMP,'** Section')
        % position within the file where the first material is mentioned
        position = ftell(fid);
    end
    % Find the key word *Solid Section
	if ~contains(TEMP,'*Solid Section')
        % Get line
        TEMP = fgetl(fid);
    else
        number_of_materials = number_of_materials + 1;
        TEMP = fgetl(fid);
	end
end

% Reset seeker position
fseek(fid,position,'bof');
TEMP = fgetl(fid);

%% Read material data (material name, element set name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter = 0;
% cell array of the material names
material_names = cell(number_of_materials,1);
% names of the element sets for the different materials
material_elem_sets_name = cell(number_of_materials,1);

while counter ~= number_of_materials && ...
      ftell(fid)~=end_position
    if ~contains(TEMP,'*Solid Section')
        TEMP = fgetl(fid);
    else
        % find commas (,) in the string
        indices1 = strfind(TEMP,',');
        % find equal signs (=) in the string
        indices2 = strfind(TEMP,'=');
        counter = counter + 1;
        
        % Names
        material_names{counter} = TEMP(indices2(end)+1:end);
        % Element sets
        material_elem_sets_name{counter} = TEMP(indices2(1)+1:indices1(2)-1);
        
        TEMP = fgetl(fid);
    end
end

%% Read material element sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reset seeker position (bof: beginning of file)
fseek(fid,0,'bof');
TEMP = fgetl(fid);

% Initialize cell array to save the element sets
material_elset = cell(number_of_materials,1);

counter = 1;
% name of the element set belonging to the first material
set_name = ['*Elset, elset=' material_elem_sets_name{1}];
while counter <= number_of_materials && ...
      ftell(fid)~=end_position
    if ~isSetHeader(TEMP,set_name)
    else
        while true
            if contains(TEMP,'generate')
                TEMP = fgetl(fid);
                C = textscan(TEMP,'%d','Delimiter',',');
                C{1,1} = (C{1,1}(1):C{1,1}(3):C{1,1}(2))';
            else
                TEMP = fgetl(fid);
                C = textscan(TEMP,'%d','Delimiter',',');
            end
        	if isempty(C{1,1})
             	break
        	end
            material_elset{counter}= [material_elset{counter} ;C{1,1}];
        end
        if counter == number_of_materials
            break
        end
        set_name = ['*Elset, elset=' material_elem_sets_name{counter+1}];
        counter = counter + 1;
        
        % Reset seeker
        fseek(fid,0,'bof');
    end
    TEMP = fgetl(fid);
end

%% Read material properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the key word ** MATERIALS
while ~contains(TEMP,'** MATERIALS') && ...
      ftell(fid)~=end_position
	% Get line
	TEMP = fgetl(fid);
end

% Find the key word *Material,
while ~contains(TEMP,'*Material,') && ...
      ftell(fid)~=end_position
	% Get line
	TEMP = fgetl(fid);
end

idx = strfind(TEMP,'=');
mat_name_current = cell(number_of_materials,1);
mat_name_current{1} = TEMP(idx + 1:end);

counter = 1;
TEMP = fgetl(fid);

% Initialize cell array with material properties
% the columns are saved in the following order
% column #1 Density (rho)
% column #2 Permittivity (kappa)
% column #3 Elastic constants (C)
% column #4 Piezoelectric coupling constants (e)
% column #5 Damping ratio (zeta)
material_properties = cell(number_of_materials, 5);

% flag if the next line should be read
flag_line = 1;

while counter <= number_of_materials && ...
      ftell(fid)~=end_position
    % Read mass density of the material
    if contains(TEMP,'*Density')
        TEMP = fgetl(fid);
        Mat = textscan(TEMP,'%f','Delimiter',',');
        material_properties{counter,1} = Mat{1,1}(1);
    % Read permittivity of the material
    elseif contains(TEMP,'*Dielectric')
        if ~contains(TEMP,'ORTHO') && ...
           ~contains(TEMP,'ENGINEERING CONSTANTS')
            % Permittivity matrix is isotropic
            TEMP = fgetl(fid);
            Mat = textscan(TEMP,'%f','Delimiter',',');
          	material_properties{counter,2} = Mat{1,1}(1);
        else
            % Permittivity matrix is orthotropic
            TEMP = fgetl(fid);
            Mat = textscan(TEMP,'%f %f %f','Delimiter',',');
            material_properties{counter,2} = [Mat{1,1}(1) ...
                                              Mat{1,2}(1) ...
                                              Mat{1,3}(1)];
        end
    % Read elastic constants of the material    
    elseif contains(TEMP,'*Elastic')
        if ~contains(TEMP,'ORTHO') && ...
           ~contains(TEMP,'ENGINEERING CONSTANTS')
           	TEMP = fgetl(fid);
        	Mat = textscan(TEMP,'%f','Delimiter',',');
            % Young's modulus
            material_properties{counter,3}(1) = Mat{1,1}(1);
            % Poisson's ratio
           	material_properties{counter,3}(2) = Mat{1,1}(2);
        elseif contains(TEMP,'ENGINEERING CONSTANTS')
            TEMP = fgetl(fid);
            Mat = textscan(TEMP,'%f','Delimiter',',');
            material_properties{counter,3} = Mat{1,1};
            while true
                TEMP = fgetl(fid);
                if contains(TEMP,'*')
                    % do not read the next line until the next loop
                    flag_line = 0;
                 	break;
                end
                Mat = textscan(TEMP,'%f','Delimiter',',');
              	material_properties{counter,3} = ...
               [material_properties{counter,3}; Mat{1,1}];
            end
 
                % Read material data
            	E1   = material_properties{counter,3}(1);
            	E2   = material_properties{counter,3}(2);
                E3   = material_properties{counter,3}(3);
             	nu12 = material_properties{counter,3}(4);
              	nu13 = material_properties{counter,3}(5);
               	nu23 = material_properties{counter,3}(6);
             	G12  = material_properties{counter,3}(7);
              	G13  = material_properties{counter,3}(8);
              	G23  = material_properties{counter,3}(9);
                        
             	% due to symmetry of C nu_ij/E_i =  nu_ji/E_j
               	nu21 = nu12*E2/E1;
              	nu31 = nu13*E3/E1;
               	nu32 = nu23*E3/E2;

               	% Compliance matrix
              	Snew = zeros(6,6);
               	Snew(1,1) = 0.5*1/E1;
              	Snew(2,2) = 0.5*1/E2;
               	Snew(3,3) = 0.5*1/E3;
               	Snew(4,4) = 0.5*1/G12;
               	Snew(5,5) = 0.5*1/G13;
               	Snew(6,6) = 0.5*1/G23;
               	Snew(1,2) = -nu21/E2;
               	Snew(1,3) = -nu31/E3;
              	Snew(2,3) = -nu32/E3;

                Snew = Snew + Snew';

             	% Stiffness matrix
              	Cnew = inv(Snew);
                        
             	material_properties{counter,3} = [Cnew(1,1); ...
               	Cnew(1,2); Cnew(2,2); Cnew(1,3); Cnew(2,3); ...
              	Cnew(3,3); Cnew(4,4); Cnew(5,5); Cnew(6,6)];
      	elseif contains(TEMP,'ORTHO')
            TEMP = fgetl(fid);
            Mat = textscan(TEMP,'%f','Delimiter',',');
            material_properties{counter,3} = Mat{1,1};
            while true
                TEMP = fgetl(fid);
                if contains(TEMP,'*')
                    % do not read the next line until the next loop
                    flag_line = 0;
                 	break;
                end
                Mat = textscan(TEMP,'%f','Delimiter',',');
              	material_properties{counter,3} = ...
               [material_properties{counter,3}; Mat{1,1}];
            end
        end
    % Read piezoelectric coupling constants of the material 
    elseif contains(TEMP,'*Piezoelectric')
        TEMP = fgetl(fid);
       	Mat = textscan(TEMP,'%f','Delimiter',',');
      	material_properties{counter,4} = Mat{1,1};
       	while true
            TEMP = fgetl(fid);
         	if contains(TEMP,'*')
                % do not read the next line until the next loop
                flag_line = 0;
                break;
         	end
            Mat = textscan(TEMP,'%f','Delimiter',',');
         	material_properties{counter,4} = ...
           [material_properties{counter,4}; Mat{1,1}];
        end
	% Read damping ratio of the material
    elseif contains(TEMP,'*Damping')
		idx = strfind(TEMP,'=');
		material_properties{counter,5} = str2double(TEMP(idx + 1:end));
    end
    
    if contains(TEMP,'*Material,')
        counter = counter + 1;
        idx = strfind(TEMP,'=');
        mat_name_current{counter} = TEMP(idx + 1:end);
    end
    
    if flag_line == 1
        % Get line
        TEMP = fgetl(fid);
    else
        % do read the next line in the next loop
        flag_line = 1;
    end
    
    % end of the material definition block
    if contains(TEMP,'** --------------------------')
        break;
    end
end

% If the order of the element types and material properties does not
% coincide
mat_props = cell(size(material_properties));
for i = 1:number_of_materials
    for j = 1:number_of_materials
        if strcmp(material_names{i}, mat_name_current{j})
            for ijk = 1:size(material_properties,2)
                mat_props{i,ijk} = material_properties{j,ijk};
            end
            break;
        else
        end
    end
end
material_properties = mat_props;

fclose('all');

end


function tf = isSetHeader(line, header)
% True if line is the header of the set definition "header", e.g.
% '*Nset, nset=Set-1' matches '*Nset, nset=Set-1, instance=Part-1' but not
% '*Nset, nset=Set-10'.
if ~ischar(line)
    tf = false;
    return
end
line = strtrim(line);
n = numel(header);
tf = strncmp(line, header, n) && (numel(line) == n || line(n+1) == ',');
end
