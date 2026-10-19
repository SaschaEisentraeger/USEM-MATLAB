%% Fixed degrees of freedom
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 02.04.2011 (University of Magdeburg)
%
% The degrees of freedom that are fixed are returned for every node in the
% set nset that is passed to this function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% Node:         Nodes with Dirichlet boundary conditions
% DoF:        	Fixed degree of freedom
%               x:  ux is fixed at the nodes
%               y:  uy is fixed at the nodes
%               xy: ux and uy are fixed at the nodes
% Topo_Node:	Nodal topology matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% dof_vector:  	Vector of fixed degrees of freedom
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dof_vector] = CreateBoundaryConditions_Node(Node,DoF)

% Constraint degrees of freedom
if strcmp(DoF,'x') == 1
	dof_vector = 2*Node-1;
elseif strcmp(DoF,'y') == 1
	dof_vector = 2*Node;
elseif strcmp(DoF,'xy') == 1    
  	dof_vector = [2*Node-1; 2*Node];
end

end