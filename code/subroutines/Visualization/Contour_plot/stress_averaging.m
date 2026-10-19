%% Stress/strain averaging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 27.11.2014 (University of Magdeburg)
%
% The stress and strain values are computed independently for each finite
% element. Due to the C0-continuity of the shape functions and therefore,
% also the displacement field the stresses and strains are not compatible
% at neighboring elements. For a smooth variation of these quantities it
% is possible to average the values from adjacent elements.
% Remark: The jumps in the stresses are a measure of the quality of the
% solution and therefore, a smooth stress distribution that is reached by
% averaging seems to be more accurate than it actually is.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% xCoor:	x-coordinates of all points where stresses and strains are
% 			computed
% yCoor:	y-coordinates of all points where stresses and strains are
% 			computed
% stress_x,stress_y,stress_xy,strain_x,strain_y,strain_xy:
% Stress and strain values saved for each element
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% stress_x,stress_y,stress_xy,strain_x,strain_y,strain_xy:
% Averaged stress and strain values saved for each element
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stress_x,stress_y,stress_xy,strain_x,strain_y,strain_xy] = ...
          stress_averaging(xCoor,yCoor,stress_x,stress_y,stress_xy,...
          strain_x,strain_y,strain_xy)

MAT = [xCoor yCoor];
[~,idx1,idx2] = unique(MAT,'rows');

% Number of nodes
nN = length(idx1);

% Initialize variables
ID = cell(nN,1);

parfor i = 1:nN
    % Indices of elements sharing one node
   	id = idx2 == i;
    ID{i} = find(id == 1);
end

for i = 1:nN
    % Index
    id = ID{i};
    % Number of values
    nID = length(id);
    % arithmetic averaging
    stress_x_sum  = sum(stress_x(id))/nID;
    stress_y_sum  = sum(stress_y(id))/nID;
    stress_xy_sum = sum(stress_xy(id))/nID;
    strain_x_sum  = sum(strain_x(id))/nID;
    strain_y_sum  = sum(strain_y(id))/nID;
    strain_xy_sum = sum(strain_xy(id))/nID;
    % write averaged values into the stress and strain vectors
    stress_x(id)  = stress_x_sum;
    stress_y(id)  = stress_y_sum;
    stress_xy(id) = stress_xy_sum;
    strain_x(id)  = strain_x_sum;
    strain_y(id)  = strain_y_sum;
    strain_xy(id) = strain_xy_sum;
end

end