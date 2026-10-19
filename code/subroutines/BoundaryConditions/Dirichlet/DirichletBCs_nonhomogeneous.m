%% Implementation of non-homogeneous Dirichlet boundary conditions
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
% K:					Stiffness matrix
% F:					Load vector
% NoDoF:                Total number of degrees of freedom
% boundary_conditions:  Data for all Dirichlet boundary conditions
% xCoor:                x-coordinates of all nodes
% yCoor:                y-coordinates of all nodes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% K:					Stiffness matrix
% F:					Load vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,F] = DirichletBCs_nonhomogeneous(K,F,NoDoF,...
                 boundary_conditions,xCoor,yCoor)

% Number of different boundary conditions that are applied to the system
number_of_BCs = size(boundary_conditions,1);

% Initialize the results
BC  = zeros(NoDoF,1);
Amp = zeros(NoDoF,1);

% Indices
idx_start_BC = 1;

%% Degrees of freedom with nonhomogeneous Dirichlet boundary conditions
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
    if strcmp(boundary_conditions(i).type,'UxUy')
        % Prescribed displacement in x- and y-direction (uniform across all
        % nodes)
        flag_uniform = 1;
        % Direction of the prescribed displacement
        flag_direction = 12;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = [2*nset-1; 2*nset];
    elseif strcmp(boundary_conditions(i).type,'Ux')
       	% Prescribed displacement in x-direction (uniform across all nodes)
        flag_uniform = 1;
        % Direction of the prescribed displacement
        flag_direction = 1;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = 2*nset-1;
    elseif strcmp(boundary_conditions(i).type,'Uy')
        % Prescribed displacement in y-direction (uniform across all nodes)
        flag_uniform = 1;
        % Direction of the prescribed displacement
        flag_direction = 2;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = 2*nset;
    elseif strcmp(boundary_conditions(i).type,'Ux fx')
        % Prescribed displacement in x-direction; spatially varying
        % displacement function
        flag_uniform = 0;
        % Direction of the prescribed displacement
        flag_direction = 1;
        % Degrees of freedom corresponding to the nodal set
        dofsBC = 2*nset-1;
    elseif strcmp(boundary_conditions(i).type,'Uy fx')
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

if ~isempty(BC)
    %% Adjust the stiffness matrix and the force vector
    % 1) Compute the force vector due to the boundary conditions (multiply
    % the column of the stiffness matrix corresponding a degree of freedom
    % a with prescribed displacement amplitude with this value) --> needs
    % to be done for all prescribed degrees of freedom. Note that the
    % stiffness matrix has to be the original one.
    if (istriu(K) || istril(K)) && ~isscalar(K)
        F_DBC = K(:,BC)*Amp + K(BC ,:)'*Amp;
    else
        F_DBC = K(:,BC)*Amp;
    end

    % 2) Save the components of the main diagonal corresponding to the
    % prescribed degrees of freedom
    Kd = diag(K(BC,BC));

    % 3) Adjust the stiffness matrix: All components of the stiffness
    % matrix in the rows and columns corresponding to the prescribed
    % degrees of freedom are set to zero.
    K(BC,:) = 0;
    K(:,BC) = 0;

    % Retrieve the diagonal components from the original stiffness matrix
    K(BC,BC) = diag(Kd);

    % 4) Compute the adjusted force vector
    F = F - F_DBC;

    % 5) Adjust the prescribed degrees of freedom in the force vector
    if (istriu(K) || istril(K)) && ~isscalar(K)
        F(BC) = 2*Kd.*Amp;
    else
        F(BC) = Kd.*Amp;
    end
end

end