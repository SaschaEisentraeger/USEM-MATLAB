%% Read boundary conditions from Abaqus-Input file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 10.02.2014 (University of Magdeburg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% inp_file:     		name of the file to be read
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% boundary_conditions:	boundary condition data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [boundary_conditions] = Read_boundary_condition_ABAQUS(inp_file)

%% Open File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(inp_file,'r');

% File length (eof: end of file)
fseek(fid,0,'eof');
end_position = ftell(fid);
    
% Reset seeker position (bof: beginning of file)
fseek(fid,0,'bof');
TEMP = fgetl(fid);

%% Read boundary condition data (BC type, amplitude, ...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the key word ** STEP:
while ~contains(TEMP,'** STEP:') && ...
      ftell(fid)~=end_position
    TEMP = fgetl(fid);
    position1 = ftell(fid);
end

% Reset seeker position (bof: beginning of file)
fseek(fid,0,'bof');
TEMP = fgetl(fid);

% Find the key word ** STEP:
while ~contains(TEMP,'** BOUNDARY CONDITIONS') && ...
 	ftell(fid)~=end_position
    TEMP = fgetl(fid);
    position2 = ftell(fid);
end

position = min([position1 position2]);

% Reset seeker position
fseek(fid,position,'bof');

% number of different boundary conditions
counter_bc = 0;

while ftell(fid)~=end_position
    TEMP = fgetl(fid);
	if contains(TEMP,'*Boundary')
        TEMP = fgetl(fid);
        if contains(TEMP,'Set')
            counter_bc = counter_bc + 1;
            TEMP = fgetl(fid);
        end
	end
end

counter_bounndary_conditions = counter_bc;

% initialize boundary conditions
boundary_conditions = struct('set_name', cell(counter_bc,1), ...
                             'amplitude', [],              	 ...
                             'type', cell(counter_bc,1),     ...
                             'nset', [],                   	 ...
                             'elset', [],                    ...
                             'step_name',cell(counter_bc,1));

% Reset seeker position
fseek(fid,position,'bof');
TEMP = fgetl(fid);

while counter_bc ~= 0
    % step name
    if contains(TEMP,'*Step')
        ia = strfind(TEMP,'=');
      	ie = strfind(TEMP,',');
        if size(ie,2) ~= 1
         	StepName =  TEMP(ia(1)+1:ie(2)-1);
        else
          	StepName =  TEMP(ia(1)+1:ie(1));
        end
    elseif contains(TEMP,'Name:')
        StepName = 'Initial';
    end
    % constant boundary conditions
    if contains(TEMP,'*Boundary') && ...
        ~contains(TEMP,'*Boundary, amplitude') && ...
        ~contains(TEMP,'*Boundary, user')
        
        TEMP = fgetl(fid);
        ie = strfind(TEMP,',');
        boundary_conditions(counter_bc).set_name = TEMP(1:ie(1)-1);
        boundary_conditions(counter_bc).step_name = StepName;
        
        %  electric potential boundary conditions
     	if contains(TEMP,', 9')
         	if contains(TEMP,', 9, 9,')
                boundary_conditions(counter_bc).type = ...
                'Electric Potential';
                % Amplitude/magnitude of the electric potential
                ie = strfind(TEMP,'.');
                ia = strfind(TEMP,',');
                boundary_conditions(counter_bc).amplitude = ...
                str2double(TEMP(ia(3)+1:ie(1)-1));
            else
                boundary_conditions(counter_bc).type = ...
                'Electric Potential 0';
                boundary_conditions(counter_bc).amplitude = 0;
          	end
     	end 
        % mechanical displacement boundary conditions
        % all displacement degrees of freedom are fixed
        if contains(TEMP,'ENCASTRE') || contains(TEMP,'PINNED')
            boundary_conditions(counter_bc).type= 'fixed';
        % symmetry boundary conditions in x-direction
        elseif contains(TEMP,'XSYMM')
        	boundary_conditions(counter_bc).type = 'XSYMM';
        % symmetry boundary conditions in y-direction    
        elseif contains(TEMP,'YSYMM')
          	boundary_conditions(counter_bc).type = 'YSYMM';
        % symmetry boundary conditions in z-direction    
        elseif contains(TEMP,'ZSYMM')
            boundary_conditions(counter_bc).type = 'ZSYMM';
        % other displacement BCs
        elseif contains(TEMP,', 1,')
        	boundary_conditions(counter_bc).type = 'Ux';
            if contains(TEMP,', 1, 1,')
                % Amplitude/magnitude of the prescribed displacement
            	ie = length(TEMP);
               	ia = strfind(TEMP,',');
             	boundary_conditions(counter_bc).amplitude = ...
               	str2double(TEMP(ia(3)+1:ie(1)));
            end
            TEMP = fgetl(fid);
            while strfind(TEMP,boundary_conditions(counter_bc).set_name)
                if contains(TEMP,', 2,')
                    boundary_conditions(counter_bc).type = ...
                        [boundary_conditions(counter_bc).type 'Uy'];
                    if contains(TEMP,', 2, 2,')
                        % Amplitude/magnitude of the prescribed displacement
                        ie = length(TEMP);
                        ia = strfind(TEMP,',');
                        boundary_conditions(counter_bc).amplitude = ...
                        [boundary_conditions(counter_bc).amplitude ...
                         str2double(TEMP(ia(3)+1:ie(1)))];
                    end
                elseif contains(TEMP,', 3,')
                    boundary_conditions(counter_bc).type = ...
                        [boundary_conditions(counter_bc).type 'Uz'];
                    if contains(TEMP,', 3, 3,')
                        % Amplitude/magnitude of the prescribed displacement
                        ie = length(TEMP);
                        ia = strfind(TEMP,',');
                        boundary_conditions(counter_bc).amplitude = ...
                        [boundary_conditions(counter_bc).amplitude ...
                         str2double(TEMP(ia(3)+1:ie(1)))];
                    end
                end
                TEMP = fgetl(fid);
            end
            
            if isempty(boundary_conditions(counter_bc).amplitude)
                boundary_conditions(counter_bc).type = ...
               [boundary_conditions(counter_bc).type ' 0'];
            end
            
            counter_bc = counter_bc - 1;
            continue
        elseif contains(TEMP,', 2,')
          	boundary_conditions(counter_bc).type = 'Uy';
         	if contains(TEMP,', 2, 2,')
             	% Amplitude/magnitude of the prescribed displacement
                ie = length(TEMP);
              	ia = strfind(TEMP,',');
             	boundary_conditions(counter_bc).amplitude = ...
               	[boundary_conditions(counter_bc).amplitude ...
               	 str2double(TEMP(ia(3)+1:ie(1)))];
         	end
            TEMP = fgetl(fid);
            while strfind(TEMP,boundary_conditions(counter_bc).set_name)
                if contains(TEMP,', 3,')
                    boundary_conditions(counter_bc).type = ...
                        [boundary_conditions(counter_bc).type 'Uz'];
                    if contains(TEMP,', 3, 3,')
                        % Amplitude/magnitude of the prescribed displacement
                        ie = length(TEMP);
                        ia = strfind(TEMP,',');
                        boundary_conditions(counter_bc).amplitude = ...
                        [boundary_conditions(counter_bc).amplitude ...
                         str2double(TEMP(ia(3)+1:ie(1)))];
                    end
                end
                TEMP = fgetl(fid);
            end
            
            if isempty(boundary_conditions(counter_bc).amplitude)
                boundary_conditions(counter_bc).type = ...
               [boundary_conditions(counter_bc).type ' 0'];
            end
            
            counter_bc = counter_bc - 1;
            continue
    	elseif contains(TEMP,', 3,')
            boundary_conditions(counter_bc).type = 'Uz';
            if contains(TEMP,', 3, 3,')
             	% Amplitude/magnitude of the prescribed displacement
                ie = length(TEMP);
              	ia = strfind(TEMP,',');
             	boundary_conditions(counter_bc).amplitude = ...
               	[boundary_conditions(counter_bc).amplitude ...
               	 str2double(TEMP(ia(3)+1:ie(1)))];
            end
            
            if isempty(boundary_conditions(counter_bc).amplitude)
                boundary_conditions(counter_bc).type = ...
               [boundary_conditions(counter_bc).type ' 0'];
            end
            
            counter_bc = counter_bc - 1;
            continue
        % boundary condition for the rotation about the z-axis means
       	% a history output for this node set is requested
        elseif contains(TEMP,', 6')
        	boundary_conditions(counter_bc).type = 'noi';
        end
        
        counter_bc = counter_bc - 1;
        TEMP = fgetl(fid);
    % time-dependent boundary conditions
    elseif contains(TEMP,'*Boundary, amplitude')
        % Time dependent boundary conditions
        TEMP = fgetl(fid);
        ie = strfind(TEMP,',');
        boundary_conditions(counter_bc).set_name = TEMP(1:ie(1)-1);
        boundary_conditions(counter_bc).step_name = StepName;
        %  electric potential boundary conditions
    	if contains(TEMP,', 9, 9,')
            boundary_conditions(counter_bc).type = 'Electric Potential dt';
          	% Amplitude/magnitude of the electric potential
           	ie = length(TEMP);
          	ia = strfind(TEMP,',');
          	boundary_conditions(counter_bc).amplitude = ...
          	str2double(TEMP(ia(3)+1:ie(1)));
        % Displacement boundary conditions
        elseif contains(TEMP,', 1, 1,')
           	boundary_conditions(counter_bc).type = 'Ux dt';
           	% Amplitude/magnitude of the prescribed displacement in
           	% x-direction
        	ie = length(TEMP);
          	ia = strfind(TEMP,',');
           	boundary_conditions(counter_bc).amplitude = ...
          	str2double(TEMP(ia(3)+1:ie(1)));
        elseif contains(TEMP,', 2, 2,')
          	boundary_conditions(counter_bc).type = 'Uy dt';
          	% Amplitude/magnitude of the prescribed displacement in
          	% y-direction
          	ie = length(TEMP);
         	ia = strfind(TEMP,',');
         	boundary_conditions(counter_bc).amplitude = ...
         	str2double(TEMP(ia(3)+1:ie(1)));
        elseif contains(TEMP,', 3, 3,')
           	boundary_conditions(counter_bc).type = 'Uz dt';
          	% Amplitude/magnitude of the prescribed displacement in
          	% z-direction
            ie = length(TEMP);
           	ia = strfind(TEMP,',');
          	boundary_conditions(counter_bc).amplitude = ...
          	str2double(TEMP(ia(3)+1:ie(1)));
    	end
        counter_bc = counter_bc - 1;
        TEMP = fgetl(fid);
	% spatially distributed boundary conditions
    elseif contains(TEMP,'*Boundary, user')
        % Time dependent boundary conditions
        TEMP = fgetl(fid);
        while contains(TEMP,'Set')
            ie = strfind(TEMP,',');
            boundary_conditions(counter_bc).set_name = TEMP(1:ie(1)-1);
            boundary_conditions(counter_bc).step_name = StepName;
            % Displacement boundary conditions
            if contains(TEMP,', 1, 1,')
                boundary_conditions(counter_bc).type = 'Ux fx';
                % Amplitude/magnitude of the prescribed displacement in
                % x-direction
                ie = length(TEMP);
                ia = strfind(TEMP,',');
                boundary_conditions(counter_bc).amplitude = ...
                str2double(TEMP(ia(3)+1:ie(1)));
            elseif contains(TEMP,', 2, 2,')
                boundary_conditions(counter_bc).type = 'Uy fx';
                % Amplitude/magnitude of the prescribed displacement in
                % y-direction
                ie = length(TEMP);
                ia = strfind(TEMP,',');
                boundary_conditions(counter_bc).amplitude = ...
                str2double(TEMP(ia(3)+1:ie(1)));
            elseif contains(TEMP,', 3, 3,')
                boundary_conditions(counter_bc).type = 'Uz fx';
                % Amplitude/magnitude of the prescribed displacement in
                % z-direction
                ie = length(TEMP);
                ia = strfind(TEMP,',');
                boundary_conditions(counter_bc).amplitude = ...
                str2double(TEMP(ia(3)+1:ie(1)));
            end
            counter_bc = counter_bc - 1;
            TEMP = fgetl(fid);
        end
    else
        TEMP = fgetl(fid);
    end
end

%% Read boundary condition element and/or node sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counter_bc = 1;
while counter_bc ~= counter_bounndary_conditions + 1
    % Reset seeker position (bof: beginning of file)
	fseek(fid,0,'bof');
 	TEMP = fgetl(fid);
  	% node set name
 	set_name = ['*Nset, nset=' boundary_conditions(counter_bc).set_name];
    % element set name
	set_name_el = ['*Elset, elset=' boundary_conditions(counter_bc).set_name];
    counter_old = counter_bc;
  	while ftell(fid)~=end_position || contains(TEMP,'*End Assembly')
        % node set
        if isSetHeader(TEMP,set_name)
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
             	boundary_conditions(counter_bc).nset= ...
               [boundary_conditions(counter_bc).nset; C{1,1}];
            end
            counter_bc = counter_bc + 1;
        end
        % element set
        if isSetHeader(TEMP,set_name_el)
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
             	boundary_conditions(counter_bc-1).elset= ...
               [boundary_conditions(counter_bc-1).elset; C{1,1}];
            end
            break;
        end 
      	TEMP = fgetl(fid);
    end
    % The node set must exist (otherwise this loop would never end)
    if counter_bc == counter_old
        fclose(fid);
        error('Read_boundary_condition_ABAQUS:setNotFound', ...
             ['Node set "%s" of a boundary condition was not found. ' ...
              'Define the set in the assembly (after *End Instance).'], ...
              boundary_conditions(counter_bc).set_name);
    end
end

fclose('all');

% If a node has several Dirichlet boundary conditions (applied in one step) 
% these are split into several conditions
NumberOfBoundaryConditions = numel(boundary_conditions);
for i = 1:NumberOfBoundaryConditions
    NumBCsCurr = numel(boundary_conditions);
    AuxBC = boundary_conditions(i);
    NumberOfAmplitudes = numel(boundary_conditions(i).amplitude);
    if NumberOfAmplitudes > 1
        for j = 1:NumberOfAmplitudes
            if j == 1
                boundary_conditions(i) = AuxBC;
                boundary_conditions(i).amplitude = AuxBC.amplitude(j);
                boundary_conditions(i).type = AuxBC.type(1+2*(j-1):2*j);
            else
                boundary_conditions(NumBCsCurr + j - 1) = AuxBC;
                boundary_conditions(NumBCsCurr + j - 1).amplitude = ...
                AuxBC.amplitude(j);
                boundary_conditions(NumBCsCurr + j - 1).type = ...
                AuxBC.type(1+2*(j-1):2*j);
            end
        end
    end
end

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
