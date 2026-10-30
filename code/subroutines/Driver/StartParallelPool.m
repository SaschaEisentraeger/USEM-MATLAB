function StartParallelPool()
%STARTPARALLELPOOL Open a parallel pool if the Parallel Computing Toolbox is available.

% Open a pool with one worker per physical core if the Parallel Computing
% Toolbox is available; without it the parfor loops run serially.
if isempty(ver('parallel')) || ~license('test','Distrib_Computing_Toolbox')
    return
end
if isempty(gcp('nocreate'))
    NCores = feature('numcores');
    maxNumCompThreads(NCores);
    parpool(NCores);
end
end
