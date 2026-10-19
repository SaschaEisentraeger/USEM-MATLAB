%% Read load conditions from Abaqus-Input file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 06.02.2014 (University of Magdeburg)
%
% The applied Neumann boundary conditions are imported from an Abaqus input
% file (*.inp). For concentrated loads and body/volume/gravity loads only
% the nodal and element sets need to be imported. In case of a distributed
% load (line load: pressure or surface traction) the procedure is
% differently. For a pressure an additional set with the Pressure has to be
% defined. This set has to be identical to the boundary where the pressure
% boundary condition is applied. When a surface traction is applied an
% additional set with the name Traction has to be defined. If multiple
% distributed loads are defined in a single model just number the
% additional sets accordingly: Pressure1, Pressure2, ..., PressureN or
% Traction1, Traction2, ..., TractionM. See also the comments in the code
% below.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% inp_file:     name of the file to be read
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% cload:            concentrated force data
% dsload:           distributed load data
% body_force:       volume load data
% gravity:          gravity data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cload, dsload, body_force, gravity] = Read_load_ABAQUS(inp_file)

%% Open File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(inp_file,'r');

% File length (eof: end of file)
fseek(fid,0,'eof');
end_position = ftell(fid);
    
% Reset seeker position (bof: beginning of file)
fseek(fid,0,'bof');
TEMP = fgetl(fid);

%% Read load data (load set_name, load name, amplitude, ...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the key word ** STEP:
while ~contains(TEMP,'** STEP:') && ...
      ftell(fid)~=end_position
    TEMP = fgetl(fid);
    position = ftell(fid);
end

for i = 1:2
    TEMP = fgetl(fid);
end

% number of different loads
counter_cload = 0;
counter_dsload = 0;
counter_body_force = 0;
counter_gravity = 0;

% number of load conditions
while ftell(fid)~=end_position
    % concentrated load (point load) - mechanical
    if contains(TEMP,'*Cload')
        counter_cload = counter_cload + 1;
    % distributed load (pressure, surface load) - mechanical
    elseif contains(TEMP,'*Dsload')
        counter_dsload = counter_dsload + 1;
    % body load (volume load) - mechanical
    elseif contains(TEMP,'Body force')
        counter_body_force = counter_body_force + 1;
    % gravitation - mechanical
    elseif contains(TEMP,'GRAV')
        counter_gravity = counter_gravity + 1;
    end
    % Read next line of the input of file
	TEMP = fgetl(fid);
end

% Reset seeker position to LOADS
fseek(fid,position,'bof');
TEMP = fgetl(fid);

% set name, amplitude, type, element set, node set
if counter_cload ~= 0
cload = struct('set_name', cell(counter_cload,1), ...
               'amplitude', [],                   ...
               'type', cell(counter_cload,1),     ...
               'nset', [],                        ...
               'step_name', cell(counter_cload,1));
end
if counter_dsload ~= 0        
dsload = struct('set_name', cell(counter_dsload,1), ...
                'amplitude', [],                    ...
                'type', cell(counter_dsload,1),     ...
                'nset', [],                         ...
                'elset', [],                        ...
                'DirectionVector', [],             	...
                'step_name', cell(counter_dsload,1));
end
if counter_body_force ~= 0            
body_force = struct('set_name', cell(counter_body_force,1), ...
                    'amplitude', [],                        ...
                    'type', cell(counter_body_force,1),     ...
                    'elset', [],                            ...
                    'step_name', cell(counter_body_force,1));
end
if counter_gravity ~= 0          
gravity = struct('set_name', cell(counter_gravity,1), ...
                 'amplitude', [],                     ...
                 'type', cell(counter_gravity,1),     ...
                 'elset', [],                         ...
                 'step_name', cell(counter_gravity,1));
end

counter_c        = 0;
counter_d        = 0;
counter_bf       = 0;
counter_g        = 0;
counter_Pressure = 0;
counter_Traction = 0;

if counter_cload ~= 0      || ...
   counter_dsload ~= 0     || ...
   counter_body_force ~= 0 || ...
   counter_gravity ~= 0
while ftell(fid)~=end_position
    % step name
    if contains(TEMP,'*Step')
        ia = strfind(TEMP,'=');
      	ie = strfind(TEMP,',');
     	if size(ie,2) ~= 1
         	StepName =  TEMP(ia(1)+1:ie(2)-1);
        else
          	StepName =  TEMP(ia(1)+1:ie(1));
     	end
    end
    % concentrated load (point load) - mechanical
    if contains(TEMP,'*Cload')
        counter_c = counter_c + 1;
        % constant point force
        if ~contains(TEMP,'amplitude')
            TEMP = fgetl(fid);
            ie = strfind(TEMP,',');
            cload(counter_c).set_name = TEMP(1:ie(1)-1);
            cload(counter_c).step_name = StepName;
                
            % magnitude of the force
            cload(counter_c).amplitude = str2double(TEMP(ie(2)+1:end));
                
            % direction of the force
            % x-direction
            if contains(TEMP,', 1,')
               cload(counter_c).type = 'f_x';
            end
            % y-direction
          	if contains(TEMP,', 2,')
             	cload(counter_c).type = 'f_y';
            end
            
        % time dependent point force
        else
          	TEMP = fgetl(fid);
        	ie = strfind(TEMP,',');
           	cload(counter_c).set_name = TEMP(1:ie(1)-1);
            cload(counter_c).step_name = StepName;
                
          	% magnitude of the force
            cload(counter_c).amplitude = str2double(TEMP(ie(2)+1:end));
                
          	% direction of the force
            % x-direction
           	if contains(TEMP,', 1,')
                cload(counter_c).type = 'f_x_dt';
            end
            % y-direction
            if contains(TEMP,', 2,')
             	cload(counter_c).type = 'f_y_dt';
            end
        end
        
    % distributed load (pressure, surface load) - mechanical
    elseif contains(TEMP,'*Dsload')
        counter_d = counter_d + 1;
        % Choose pressure as load and at the same time define also a set
      	% with the name Pressure corresponding to the region where the
      	% load was applied. If you define more than one pressure load
      	% number them defined sets consecutively: Pressure1, Pressure2, ...
        % The same procedure should be applied to surface tractions. Name
        % the defined sets accordingly: Traction1, Traction2, ...
        if ~contains(TEMP,'amplitude')
            % Read the next line
            TEMP = fgetl(fid);
            % Surface traction or pressure load
            if contains(TEMP,'TRSHR') || ...
               contains(TEMP,'TRVEC')
                counter_Traction = counter_Traction + 1;
                % User defined spatial variation of the surface traction
                if contains(TEMP,'NU')
                    dsload(counter_d).type = 'Traction fx';
                else
                    % Uniform traction
                    dsload(counter_d).type = 'Traction';
                end
                % Assign set name and step name
                dsload(counter_d).set_name = ['Traction' num2str(counter_Traction)];
                dsload(counter_d).step_name = StepName;
            else
                counter_Pressure = counter_Pressure + 1;
                % User defined spatial variation of the pressure field
                if contains(TEMP,'NU')
                    dsload(counter_d).type = 'Pressure fx';
                else
                    % Uniform pressure
                    dsload(counter_d).type = 'Pressure';
                end
                % Assign set name and step name
                dsload(counter_d).set_name = ['Pressure' num2str(counter_Pressure)];
                dsload(counter_d).step_name = StepName; 
            end
        else
            % Read the next line
            TEMP = fgetl(fid);
            % Surface traction or pressure load
            if contains(TEMP,'TRSHR') || ...
               contains(TEMP,'TRVEC')
                counter_Traction = counter_Traction + 1;
                dsload(counter_d).type = 'Traction dt';
                % Assign set name and step name
                dsload(counter_d).set_name = ['Traction' num2str(counter_Traction)];
                dsload(counter_d).step_name = StepName; 
            else
                counter_Pressure = counter_Pressure + 1;
                dsload(counter_d).type = 'Pressure dt';
                % Assign set name and step name
                dsload(counter_d).set_name = ['Pressure' num2str(counter_Pressure)];
                dsload(counter_d).step_name = StepName; 
            end
        end
        
        ia = strfind(TEMP,',');
        if numel(ia) > 2
            % Assign the direction vector for the surface traction
            dsload(counter_d).DirectionVector = ...
            [str2double(TEMP(ia(3)+1:ia(4)));...
             str2double(TEMP(ia(4)+1:ia(5)));...
             str2double(TEMP(ia(5)+1:end))];
            % Amplitude/ magnitude of the pressure/surface load
            dsload(counter_d).amplitude = str2double(TEMP(ia(2)+1:ia(3)));
        else
            % Amplitude/ magnitude of the pressure/surface load
            dsload(counter_d).amplitude = str2double(TEMP(ia(2)+1:end));
        end

    % body load (volume load) - mechanical
    elseif contains(TEMP,'Body force') || contains(TEMP,'Dload') && ...
           counter_bf ~= counter_body_force
        counter_bf = counter_bf + 1;
		
		if contains(TEMP,'Body force') 
            for i = 1:2
                TEMP = fgetl(fid);
            end
        else
            TEMP = fgetl(fid);
        end
    	ie = strfind(TEMP,',');
     	body_force(counter_bf).set_name = TEMP(1:ie(1)-1);
        body_force(counter_bf).step_name = StepName;
            
      	body_force(counter_bf).type = 'Body_force';
		body_force(counter_bf).amplitude = zeros(2,1);
            
      	for i = 1:2
            ie = length(TEMP);
            if contains(TEMP,'**')
                break;
            elseif contains(TEMP,'BX')
                ia = strfind(TEMP,',');
                % magnitude in x-direction
                body_force(counter_bf).amplitude(1,1) = ...
                    str2double(TEMP(ia(2):ie(1)));
            elseif contains(TEMP,'BY')
                ia = strfind(TEMP,',');
                % magnitude in y-direction
                body_force(counter_bf).amplitude(2,1) = ...
                    str2double(TEMP(ia(2):ie(1)));
            end
            TEMP = fgetl(fid);
        end

    % gravitation - mechanical
    elseif contains(TEMP,'GRAV')
        counter_g = counter_g + 1;

        gravity(counter_g).type = 'Body_force';
        gravity(counter_g).set_name = 'whole';
        gravity(counter_g).step_name = StepName;
            
      	ia = strfind(TEMP,',');
        % gravity magnitude
      	gravity(counter_g).amplitude(1,1) = str2double(TEMP(ia(2)+1:ia(3)-1));
       	gravity(counter_g).amplitude(2,1) = str2double(TEMP(ia(3)+1:ia(4)-1));
       	gravity(counter_g).amplitude(3,1) = str2double(TEMP(ia(4)+1:end));

    end
    
    if counter_c == counter_cload       && ...
       counter_d == counter_dsload      && ...
       counter_bf == counter_body_force && ...
       counter_g == counter_gravity
        break;
    end
    
    % Read next line of the input of file
	TEMP = fgetl(fid);
end
end

%% Read load element and/or node sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of element/node sets to be saved
number_of_sets = counter_cload      + ...
                 counter_dsload     + ...
                 counter_body_force + ...
                 counter_gravity;

for i = 1:number_of_sets
    % concentrated load - element/node set
	if counter_c ~= 0
        % Reset seeker position (bof: beginning of file)
        fseek(fid,0,'bof');
        TEMP = fgetl(fid);
        % node set name
        set_name = ['*Nset, nset=' cload(counter_c).set_name];
        while ftell(fid)~=end_position
            if ~contains(TEMP,set_name)
            else
                while true
                    TEMP = fgetl(fid);
                    C = textscan(TEMP,'%d','Delimiter',',');
                    if isempty(C{1,1})
                        break
                    end
                    cload(counter_c).nset= [cload(counter_c).nset; C{1,1}];
                end
                break;
            end
            TEMP = fgetl(fid);
        end
        counter_c = counter_c - 1;

    % distributed load - element/node set
    elseif counter_d ~= 0
        % Reset seeker position (bof: beginning of file)
        fseek(fid,0,'bof');
        TEMP = fgetl(fid);
        % node set name
        set_name = ['*Nset, nset=' dsload(counter_d).set_name];
        % element set name
        set_name_el = ['*Elset, elset=' dsload(counter_d).set_name];
        while ftell(fid)~=end_position
            % node set
            if contains(TEMP,set_name)
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
                    dsload(counter_d).nset= ...
                   [dsload(counter_d).nset; C{1,1}];
                end
            end
            % element set
            if contains(TEMP,set_name_el)
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
                    dsload(counter_d).elset= ...
                   [dsload(counter_d).elset; C{1,1}];
                end
                break;
            end 
            TEMP = fgetl(fid);
        end
        counter_d = counter_d - 1;

    % body force - element/node set
    elseif counter_bf ~= 0
        % Reset seeker position (bof: beginning of file)
        fseek(fid,0,'bof');
        TEMP = fgetl(fid);
        % element set name
        set_name = ['*Elset, elset=' body_force(counter_bf).set_name];
        while ftell(fid)~=end_position
            if ~contains(TEMP,set_name)
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
                    body_force(counter_bf).elset= [body_force(counter_bf).elset; C{1,1}];
                end
                break;
            end
            TEMP = fgetl(fid);
        end
        counter_bf = counter_bf - 1;

    % gravity - element/node set
    elseif counter_g ~= 0
        % Reset seeker position (bof: beginning of file)
        fseek(fid,0,'bof');
        TEMP = fgetl(fid);
        % element set name
        set_name = ['*Elset, elset=' gravity(counter_g).set_name];
        while ftell(fid)~=end_position
            if ~contains(TEMP,set_name)
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
                    gravity(counter_g).elset= [gravity(counter_g).elset; C{1,1}];
                end
                break;
            end
            TEMP = fgetl(fid);
        end
        counter_g = counter_g - 1;
	end
end

if ~exist('cload','var')
    cload = [];
end

if ~exist('dsload','var')
    
    dsload = [];
end

if ~exist('body_force','var')
    body_force = [];
end

if ~exist('gravity','var')
    gravity = [];
end

fclose('all');

end