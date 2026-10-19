%% Implementation of non-homogeneous Dirichlet BCs - Dynamics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 10.04.2018 (University of Magdeburg)
%
% Nonhomogeneous Dirichlet boundary conditions are implemented. Hence, the
% stiffness matrix and the force vector are adjusted in order to account
% for prescribed displacements at the nodes of the spectral element.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% K:                  	Global stiffness matrix
% C:                  	Global damping matrix
% M:                   	Global mass matrix
% Topology_Node:      	Nodal topology matrix
% boundary_conditions:  Data for all Dirichlet boundary conditions
% xCoor:                x-coordinates of all nodes
% yCoor:                y-coordinates of all nodes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% DOF_C:                Prescribed degrees of freedom and amplitudes of the
%                   	nonhomogeneous Dirichlet boundary conditions
%                   	(time-independent)
% DOF_C_Dyn:            Prescribed degrees of freedom and amplitudes of the
%                   	nonhomogeneous Dirichlet boundary conditions
%                   	(time-dependent)
% F:                	Modified force vector (nonhomogeneous Dirichlet
%                       boundary conditions; time independent)
% Fdt:                  Modified force vector (nonhomogeneous Dirichlet
%                       boundary conditions; time dependent)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [DOF_C,DOF_C_Dyn,F,Fdt] = ...
          DirichletBCs_nonhomogeneous_dynamic(K,C,M,Topology_Node,...
          boundary_conditions,xCoor,yCoor)
     
% Degrees of freedom topology matrix
[Topology] = Topology_DoF_2D(Topology_Node);

% Number of degrees of freedom
NoDoF = max(max(Topology));

% Nonhomogeneous Dirichlet boundary conditions
%-------------------------------------------------------------------------%
% In the current implementation it is assumed that only nonhomogeneous
% (time-dependent) displacement boundary conditions are prescribed. If for
% any reason velocity or acceleration BCs should be needed it is possible
% to integrate these functions and to provide the adequate displacement BC.

% Number of different boundary conditions that are applied to the system
number_of_BCs = size(boundary_conditions,1);

% Initialize the results
BC  = zeros(NoDoF,1);
Amp = zeros(NoDoF,1);

% Indices
idx_start_BC = 1;

%% Degrees of freedom with nonhomogeneous Dirichlet BCs - Time-dependent
% First we need to determine the degrees of freedom with nonhomogeneous
% Dirichlet boundary conditions. --> Note that there is a distinct
% difference between spectral or h- and p-elements. In the case of spectral
% or h-elements all degrees of freedom retain a physical meaning (nodal
% shape functions), whereas for p-elements only the degrees of freedom
% connected to the corner nodes hold a physical interpretation as nodal
% displacements while all other degrees of freedom are unknown parameters
% of the chosen ansatz (modal shape functions).
for i = 1:number_of_BCs
    % Nodal set with Dirichlet boundary conditions
 	nset = boundary_conditions(i).nset;
    
    % Dirichlet boundary conditions for mechanical degrees of freedom
    if strcmp(boundary_conditions(i).type,'UxUy dt')
        % Prescribed displacement in x- and y-direction (uniform across all
        % nodes)
        flag_uniform = 1;
        % Direction of the prescribed displacement
        flag_direction = 12;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = [2*nset-1; 2*nset];
    elseif strcmp(boundary_conditions(i).type,'Ux dt')
       	% Prescribed displacement in x-direction (uniform across all nodes)
        flag_uniform = 1;
        % Direction of the prescribed displacement
        flag_direction = 1;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = 2*nset-1;
    elseif strcmp(boundary_conditions(i).type,'Uy dt')
        % Prescribed displacement in y-direction (uniform across all nodes)
        flag_uniform = 1;
        % Direction of the prescribed displacement
        flag_direction = 2;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = 2*nset;
    elseif strcmp(boundary_conditions(i).type,'Ux dt fx')
        % Prescribed displacement in x-direction; spatially varying
        % displacement function
        flag_uniform = 0;
        % Direction of the prescribed displacement
        flag_direction = 1;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = 2*nset-1;
    elseif strcmp(boundary_conditions(i).type,'Uy dt fx')
        % Prescribed displacement in y-direction; spatially varying
        % displacement function
        flag_uniform = 0;
        % Direction of the prescribed displacement
        flag_direction = 2;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = 2*nset;
    else
        continue;
    end
    
    % Save the variables (start and end index for saving the degrees of
    % freedom on a vector)
   	idx_end_BC = idx_start_BC - 1 + numel(dofsBC);
    % Save the degrees of freedom in a vector
  	BC(idx_start_BC:idx_end_BC) = dofsBC;
    if flag_uniform == 1
        % Save the amplitude in a vector
        Amp(idx_start_BC:idx_end_BC) = boundary_conditions(i).amplitude;
    else
        % Amplitude of the displacement function
        A(1) = boundary_conditions(i).amplitude;
        A(2) = flag_direction;
        % Coordinates of the nodes
        Points = [xCoor(nset), yCoor(nset)];
        % Compute the displacement values
        [DispVal] = SpatialDisplacementFunction(A,Points);
        % Save the amplitude in a vector
        Amp(idx_start_BC:idx_end_BC) = DispVal;
    end
    % Update the indices
	idx_start_BC  = idx_end_BC + 1;
end 

% Delete all duplicate entries in BC and adjust Amp accordingly
[BC,idx,~] = unique(BC,'stable');
Amp = Amp(idx);
% Delete all zero entries in BC and adjust Amp accordingly
idx = BC ~= 0;
BC = BC(idx);
Amp = Amp(idx);

% Dynamic
DOF_C_Dyn.BC  = BC;
DOF_C_Dyn.Amp = Amp;

% In this implementation it is assumed that all prescribed displacements
% vary with the same function in time. It is only possible to apply
% nonhomogeneous displacement boundary conditions (note that velocity and
% acceleration can be integrated to obtain the adequate displacement BC).
% Force vectors due to nonhomogeneous, time-dependent Dirichlet BCs
% Adjust the force vector - Stiffness matrix
if ~isempty(Amp)
    if (istriu(K) || istril(K)) && ~isscalar(K)
        Fdt.K = -(K(:,BC)*Amp + K(BC,:)'*Amp);
    else
        Fdt.K = -K(:,BC)*Amp;
    end
else
    Fdt.K = zeros(NoDoF,1);
end

% Adjust the force vector - Damping matrix
if ~isempty(Amp)
    if istriu(C) || istril(C)
        Fdt.C = -(C(:,BC)*Amp + C(BC,:)'*Amp);
    else
        Fdt.C = -C(:,BC)*Amp;
    end
else
    Fdt.C = zeros(NoDoF,1);
end

% Adjust the force vector - Mass matrix
if ~isempty(Amp)
    if istriu(M) || istril(M)
        Fdt.M = -(M(:,BC)*Amp + M(BC,:)'*Amp);
    else
        Fdt.M = -M(:,BC)*Amp;
    end
else
    Fdt.M = zeros(NoDoF,1);
end

%% Degrees of freedom with nonhomogeneous Dirichlet BCs - Time-independent
% Initialize the results
BC  = zeros(NoDoF,1);
Amp = zeros(NoDoF,1);

% Indices
idx_start_BC = 1;

for i = 1:number_of_BCs
    % Nodal set with Dirichlet boundary conditions
 	nset = boundary_conditions(i).nset;
    
    % Dirichlet boundary conditions for mechanical degrees of freedom
    if strcmp(boundary_conditions(i).type,'UxUy')
        % Prescribed displacement in x- and y-direction (uniform across all
        % nodes)
        flag_uniform = 1;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = [2*nset-1; 2*nset];
    elseif strcmp(boundary_conditions(i).type,'Ux')
       	% Prescribed displacement in x-direction (uniform across all nodes)
        flag_uniform = 1;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = 2*nset-1;
    elseif strcmp(boundary_conditions(i).type,'Uy')
        % Prescribed displacement in y-direction (uniform across all nodes)
        flag_uniform = 1;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = 2*nset;
    elseif strcmp(boundary_conditions(i).type,'Ux fx')
        % Prescribed displacement in x-direction; spatially varying
        % displacement function
        flag_uniform = 0;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = 2*nset-1;
    elseif strcmp(boundary_conditions(i).type,'Uy fx')
        % Prescribed displacement in y-direction; spatially varying
        % displacement function
        flag_uniform = 0;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = 2*nset;
    else
        continue;
    end
    
    % Save the variables (start and end index for saving the degrees of
    % freedom on a vector)
   	idx_end_BC = idx_start_BC - 1 + numel(dofsBC);
    % Save the degrees of freedom in a vector
  	BC(idx_start_BC:idx_end_BC) = dofsBC;
    if flag_uniform == 1
        % Save the amplitude in a vector
        Amp(idx_start_BC:idx_end_BC) = boundary_conditions(i).amplitude;
    else
        % Amplitude of the displacement function
        A = boundary_conditions(i).amplitude;
        % Coordinates of the nodes
        Points = [xCoor(nset), yCoor(nset)];
        % Compute the displacement values
        [DispVal] = SpatialDisplacementFunction(A,Points);
        % Save the amplitude in a vector
        Amp(idx_start_BC:idx_end_BC) = DispVal;
    end
    % Update the indices
	idx_start_BC  = idx_end_BC + 1;
end 

% Delete all duplicate entries in BC and adjust Amp accordingly
[BC,idx,~] = unique(BC,'stable');
Amp = Amp(idx);
% Delete all zero entries in BC and adjust Amp accordingly
idx = BC ~= 0;
BC = BC(idx);
Amp = Amp(idx);

% Force vectors due to nonhomogeneous, time-dependent Dirichlet BCs
% Adjust the force vector - Stiffness matrix
if ~isempty(Amp)
    if (istriu(K) || istril(K)) && ~isscalar(K)
        F = -(K(:,BC)*Amp + K(BC,:)'*Amp);
    else
        F = -K(:,BC)*Amp;
    end
else
    F = zeros(NoDoF,1);
end

% Constant
DOF_C.BC  = BC;
DOF_C.Amp = Amp;

end