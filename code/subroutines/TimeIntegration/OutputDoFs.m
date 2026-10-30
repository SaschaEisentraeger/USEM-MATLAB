function dof_poi = OutputDoFs(NOI,save_var,ActiveDoF,num_DoF)
%OUTPUTDOFS DOFs of the reduced system whose time history is returned.

% DOFs of the reduced system whose history is returned
if save_var == 1
    dof_poi = 1:numel(ActiveDoF);
    return
end

% Both DOFs of the nodes of interest, numbered in the reduced system
dof_poi = zeros(2*numel(NOI),1);
dof_poi(1:2:end) = 2*NOI-1;
dof_poi(2:2:end) = 2*NOI;
vecID = zeros(num_DoF,1);
vecID(ActiveDoF) = 1:length(ActiveDoF);
dof_poi = vecID(dof_poi);
dof_poi(dof_poi == 0) = [];

if save_var == 0 && isempty(dof_poi)
    % No node of interest: save the first four active DOFs
    dof_poi = [1 2 3 4];
    warning(['No nodes of interest have been defined. Therefore, the ' ...
             'degrees of freedom DOF = [' num2str(ActiveDoF(1)) ','   ...
             num2str(ActiveDoF(2)) ',' num2str(ActiveDoF(3)) ','     ...
             num2str(ActiveDoF(4)) '] are saved.'])
end
end
