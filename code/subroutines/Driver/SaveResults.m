function file = SaveResults(model,settings,results)
%SAVERESULTS Save model, settings and results of a simulation.
%
%   file = SaveResults(model, settings, results) writes the three structs
%   to a MAT file (v7.3) in settings.ResultsFolder and returns its name:
%
%     static            Static/<name>_Static/<name>_ndof<N>.mat
%     dynamic           Dynamic/<name>_<Method>/<name>_ndof<N>_<mass>.mat
%     modal             Modal/<name>_Modal/<name>_ndof<N>_<mass>.mat
%     frequency domain  Frequency_Domain/<name>_Frequency_Domain/
%                       <name>_ndof<N>_<mass>.mat
%
%   <mass> is CMM (consistent mass matrix), LMM_NQ (nodal quadrature),
%   LMM_RowSum or LMM_HRZ. A dynamic analysis that saved the nodes of
%   interest (save_var = 0) also writes <...>_Uy.mat with the time vector
%   T and the y-displacements U of these nodes, sorted by x.

name = model.name;
N = num2str(model.num_DoF);
massTag = {'_CMM','_LMM_NQ','_LMM_RowSum','_LMM_HRZ'};
massTag = massTag{settings.mass_matrix};

switch settings.analysis_type
    case 'Static_analysis'
        folder = fullfile('Static',[name '_Static']);
        base = [name '_ndof' N];
        massTag = '';
    case 'Dynamic_analysis'
        folder = fullfile('Dynamic',[name '_' num2str(settings.TimeIntegration.Method)]);
        base = [name '_ndof' N];
    case 'Modal_analysis'
        folder = fullfile('Modal',[name '_Modal']);
        base = [name '_ndof' N];
    case 'FrequencyDomain_analysis'
        folder = fullfile('Frequency_Domain',[name '_Frequency_Domain']);
        base = [name '_ndof' N];
end

folder = fullfile(settings.ResultsFolder,folder);
if ~exist(folder,'dir')
    mkdir(folder);
end
file = fullfile(folder,[base massTag '.mat']);
save(file,'model','settings','results','-v7.3');

% Dynamic analysis: y-displacements of the nodes of interest, sorted by x
% (only if both components of every node of interest have been stored)
if strcmp(settings.analysis_type,'Dynamic_analysis') && ...
   settings.save_var == 0 && size(results.y_comp,2) == 2*numel(results.NOI)
    [~,idx] = sort(model.xCoor(results.NOI));
    U = results.y_comp(:,2*idx,:);
    T = results.t_comp;
    save(fullfile(folder,[base '_Uy.mat']),'T','U','-v7.3');
end

end
