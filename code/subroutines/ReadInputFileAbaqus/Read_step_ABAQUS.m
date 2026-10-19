%% Read step definition from Abaqus-Input file
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
% step_name:        Name of the analysis step
% step_type:        Analysis type
% tend:             Simulation time for a dynamic analysis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [step_name,step_type,tend,counterSteps] = Read_step_ABAQUS(inp_file)

%% Open File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(inp_file,'r');

% File length (eof: end of file)
fseek(fid,0,'eof');
end_position = ftell(fid);
    
% Reset seeker position (bof: beginning of file)
fseek(fid,0,'bof');
TEMP = fgetl(fid);

%% Read step definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Count number of steps
numSteps = 0;
while ftell(fid)~=end_position
    % Find the key word *Step,
	if contains(TEMP,'*Step')
     	numSteps = numSteps + 1;
 	end
	TEMP = fgetl(fid);
end

% Reset seeker position (bof: beginning of file)
fseek(fid,0,'bof');
TEMP = fgetl(fid);

counterSteps = 0;
step_name = cell(numSteps,1);
step_type = cell(numSteps,1);
tend = zeros(numSteps,1);
while ftell(fid)~=end_position
    % Find the key word *Step,
	if ~contains(TEMP,'*Step')
        TEMP = fgetl(fid);
    else
        counterSteps = counterSteps + 1;
        ia = strfind(TEMP,'=');
      	ie = strfind(TEMP,',');
     	if size(ie,2) ~= 1
         	step_name{counterSteps} =  TEMP(ia(1)+1:ie(2)-1);
        else
          	step_name{counterSteps} =  TEMP(ia(1)+1:ie(1));
     	end
            
    	TEMP = fgetl(fid);
		
		% Skip the load description line
        if ~contains(TEMP,'*')
            TEMP = fgetl(fid);
        end
            
      	ia = strfind(TEMP,'*');
       	ie = strfind(TEMP,',');
      	if isempty(ie)
          	step_type{counterSteps} =  TEMP(ia(1)+1:end);
        else
          	step_type{counterSteps} =  TEMP(ia(1)+1:ie(1)-1);
        end
            
     	if strcmp(step_type{counterSteps},'Dynamic')
           	TEMP = fgetl(fid);
        	C = textscan(TEMP,'%f %f %f %f','delimiter', ',');
          	tend(counterSteps) = C{2};
        else
            tend = 0;
    	end 
        
        TEMP = fgetl(fid);
        
        if counterSteps == numSteps
            break
        end
 	end
end

fclose('all');

end