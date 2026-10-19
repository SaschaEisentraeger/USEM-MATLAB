%% Concentrated (point) force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 03.12.2010 (University of Magdeburg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% Nodes:    	Nodal numbers of the nodes with a point force
% num_DoF:      Number of degrees of freedom
% amp:          Magnitude of the point force
% DoF:        	Direction of the point force
%                   = x - point force in (global) x-direction
%                   = y - point force in (global) y-direction
%                   = xy - point force in (global) x- and y-direction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% F_mech:  		External load vector
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F_mech] = CreateForce(Nodes,num_DoF,amp,DoF)

% Degrees of freedom
[dof_vector] = CreateBoundaryConditions_Node(Nodes,DoF);

% External force vector
F_mech = zeros(num_DoF,1);

% Assign point forces
F_mech(dof_vector) = amp;

end