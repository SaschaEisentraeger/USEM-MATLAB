%% Loads
% Loads are added to the FE-system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 12.02.2014 (University of Magdeburg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% p:                 	Polynomial degree structure
% p.pxi:                Polynomial degree in xi
% p.peta:           	Polynomial degree in eta
% p.ndist:              Nodal distribution
% cload:              	Concentrated loads (point forces)
% dsload:             	Distributed loads (pressure, traction)
% body_force:         	Volume loads (arbitrary direction and elements)
% gravity:            	Volume load (all elements)
% mat:                	Material data
% num_DoF:            	Number of degrees of freedom
% Topology_Node:        Nodal topology (connectivity/coincidence) matrix
% xCoor:            	x-coordinates of all nodes
% yCoor:               	y-coordinates of all nodes
% num_Nodes:         	Number of nodes
% flag_mRHS:            Flag for multiple load cases; if multiple right-
%                       hand sides are available the value must be 'true'
%                       and each load is saved in a separate column of the
%                       force "matrix".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% F:                  	Load vector (time-independent)
% F_dt:               	Load vector (time-dependent)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,F_dt] = NeumannBCs(p,cload,dsload,body_force,gravity,mat,...
                    num_DoF,Topology_Node,xCoor,yCoor,num_Nodes,flag_mRHS)

% Check for multiple load cases (right-hand sides)
if ~exist('flag_mRHS','var')
    flag_mRHS = false;
    % Only one load case is considered
    nLC = 1;
end

% Number of different loads that are applied to the system
% concentrated loads
if isempty(cload)
    number_of_cloads = 0;
else
    number_of_cloads = size(cload,1);
end
% distributed loads (surface)
if isempty(dsload)
    number_of_dsloads = 0;
else
    number_of_dsloads = size(dsload,1);
end
% distributed loads (volume)
if isempty(body_force)
    number_of_body_forces = 0;
else
    number_of_body_forces = size(body_force,1);
end
% gravity
if isempty(gravity)
    number_of_gravities = 0;
else
    number_of_gravities = size(gravity,1);
end

% Number of load cases
if flag_mRHS
    nLC = number_of_cloads + number_of_dsloads + ...
          number_of_body_forces + number_of_gravities;
end

% Initialize force vector
F = zeros(num_DoF,nLC);
F_dt = zeros(num_DoF,nLC);

%% implement concentrated loads
counterLC = 0;
for i = 1:number_of_cloads
    counterLC = counterLC + 1;
    % apply a force in x-direction
    if strcmp(cload(i).type,'f_x') || strcmp(cload(i).type,'f_x_dt')
        nset = cload(i).nset;
        if size(nset,1) == 3
            nset1 = (nset(1):nset(3):nset(2))';
            if nset1(end) == nset(2)
                nset = nset1;
            end
        end
        nset = nset(nset <= num_Nodes);
        
        [F_mech] = CreateForce(nset,num_DoF,cload(i).amplitude,'x');
        
        if flag_mRHS
            if strcmp(cload(i).type,'f_x')
                F(1:num_DoF,counterLC) = F_mech;
            else
                F_dt(1:num_DoF,counterLC) = F_mech;
            end
        else
            if strcmp(cload(i).type,'f_x')
                F(1:num_DoF,1) = F(1:num_DoF,1) + F_mech;
            else
                F_dt(1:num_DoF,1) = F_dt(1:num_DoF,1) + F_mech;
            end
        end
    % apply a force in y-direction
    elseif strcmp(cload(i).type,'f_y') || strcmp(cload(i).type,'f_y_dt')
        nset = cload(i).nset;
        if size(nset,1) == 3
            nset1 = (nset(1):nset(3):nset(2))';
            if nset1(end) == nset(2)
                nset = nset1;
            end
        end
        nset = nset(nset <= num_Nodes);
        
        [F_mech] = CreateForce(nset,num_DoF,cload(i).amplitude,'y');
        
        if flag_mRHS
            if strcmp(cload(i).type,'f_y')
                F(1:num_DoF,counterLC) = F_mech;
            else
                F_dt(1:num_DoF,counterLC) = F_mech;
            end
        else
            if strcmp(cload(i).type,'f_y')
                F(1:num_DoF,1) = F(1:num_DoF,1) + F_mech;
            else
                F_dt(1:num_DoF,1) = F_dt(1:num_DoF,1) + F_mech;
            end
        end
    end
end

%% implement distributed loads (surface)
for i = 1:number_of_dsloads
    counterLC = counterLC + 1;
    % Pressure load
    if contains(dsload(i).type,'Pressure')
        nset = dsload(i).nset;
        elset = dsload(i).elset;
        if size(nset,1) == 3
            nset1 = (nset(1):nset(3):nset(2))';
            if nset1(end) == nset(2)
                nset = nset1;
            end
        end
        if size(elset,1) == 3
            elset1 = (elset(1):elset(3):elset(2))';
            if elset1(end) == elset(2)
                elset = elset1;
            end
        end
        [NSET] = BC_Edge_sub(p,Topology_Node,nset,elset);
        
        % Uniform or spatially varying pressure distribution
        if contains(dsload(i).type,'fx')
            flag_uni = 0;
        else
            flag_uni = 1;
        end
        
        % Compute energetically equivalent nodal forces (pressure)
        [F_mech] = CreatePress(p,xCoor,yCoor,num_DoF,Topology_Node,...
		           NSET,dsload(i).amplitude,elset,flag_uni);
    
        if flag_mRHS
            if strcmp(dsload(i).type,'Pressure') || ...
               strcmp(dsload(i).type,'Pressure fx')
                F(1:num_DoF,counterLC) = F_mech;
            else
                F_dt(1:num_DoF,counterLC) = F_mech;
            end
        else
            if strcmp(dsload(i).type,'Pressure') || ...
               strcmp(dsload(i).type,'Pressure fx')
                F(1:num_DoF,1) = F(1:num_DoF,1) + F_mech;
            else
                F_dt(1:num_DoF,1) = F_dt(1:num_DoF,1) + F_mech;
            end
        end
    % Surface traction
    elseif contains(dsload(i).type,'Traction')
        nset = dsload(i).nset;
        elset = dsload(i).elset;
        if size(nset,1) == 3
            nset1 = (nset(1):nset(3):nset(2))';
            if nset1(end) == nset(2)
                nset = nset1;
            end
        end
        if size(elset,1) == 3
            elset1 = (elset(1):elset(3):elset(2))';
            if elset1(end) == elset(2)
                elset = elset1;
            end
        end
        [NSET] = BC_Edge_sub(p,Topology_Node,nset,elset);

        % Uniform or spatially varying surface traction
        if contains(dsload(i).type,'fx')
            flag_uni = 0;
        else
            flag_uni = 1;
        end
        
        % Compute energetically equivalent nodal forces (surface traction)
        [F_mech] = CreateTraction(p,xCoor,yCoor,num_DoF,Topology_Node,...
                   NSET,dsload(i).amplitude,...
                   dsload(i).DirectionVector(1:2),elset,flag_uni);
        
        if flag_mRHS
            if strcmp(dsload(i).type,'Traction') || ...
               strcmp(dsload(i).type,'Traction fx')
                F(1:num_DoF,counterLC) = F_mech;
            else
                F_dt(1:num_DoF,counterLC) = F_mech;
            end
        else
            if strcmp(dsload(i).type,'Traction') || ...
               strcmp(dsload(i).type,'Traction fx')
                F(1:num_DoF,1) = F(1:num_DoF,1) + F_mech;
            else
                F_dt(1:num_DoF,1) = F_dt(1:num_DoF,1) + F_mech;
            end
        end
    end
end

%% implement distributed loads (volume)
for i = 1:number_of_body_forces
    counterLC = counterLC + 1;
    [F_mech] = CreateBload(p,xCoor,yCoor,num_DoF,Topology_Node,...
	           body_force);
    
	if flag_mRHS
		if strcmp(body_force(i).type,'Body_force')
			F(1:num_DoF,counterLC) = F_mech;
		else
			F_dt(1:num_DoF,counterLC) = F_mech;
		end
	else
		if strcmp(body_force(i).type,'Body_force')
			F(1:num_DoF,1) = F(1:num_DoF,1) + F_mech;
		else
			F_dt(1:num_DoF,1) = F_dt(1:num_DoF,1) + F_mech;
		end
	end
end

%% implement gravity
for i = 1:number_of_gravities
    counterLC = counterLC + 1;
    [F_mech] = CreateGravity(p,xCoor,yCoor,num_DoF,Topology_Node,...
	           gravity,mat);
	if flag_mRHS
		if strcmp(body_force(i).type,'Body_force')
			F(1:num_DoF,counterLC) = F_mech;
		else
			F_dt(1:num_DoF,counterLC) = F_mech;
		end
	else
		if strcmp(body_force(i).type,'Body_force')
			F(1:num_DoF,1) = F(1:num_DoF,1) + F_mech;
		else
			F_dt(1:num_DoF,1) = F_dt(1:num_DoF,1) + F_mech;
		end
	end
end

% %% Save memory
% F = sparse(F);
% F_dt = sparse(F_dt);

end