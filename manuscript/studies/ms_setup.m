function [root, dataDir] = ms_setup()
%MS_SETUP Paths for the studies of the manuscript (MATLAB or GNU Octave).
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(fullfile(root, 'code', 'subroutines')));
addpath(fullfile(root, 'tests', 'helpers'));
% vectorized element routines (cyl_Ke_PM, cyl_Ke_PP, cyl_quad), used by ms_Ke
addpath(fullfile(root, 'studies', 'thick_cylinder', 'matlab'));
% Several routines of the code base contain parfor loops. The studies run
% as several MATLAB batch sessions side by side (run_*.sh), so no session
% may start a pool of its own: the parfor loops then run serially. The
% setting is persistent; run_*.sh restore it after the last session.
if ~exist('OCTAVE_VERSION', 'builtin') && batchStartupOptionUsed ...
        && ~isempty(ver('parallel')) && isempty(gcp('nocreate'))
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;
end
if exist('OCTAVE_VERSION', 'builtin')
    addpath(fullfile(root, 'studies', 'metric_basis', 'octave'));
    warning('off', 'all');
end
dataDir = fullfile(fileparts(here), 'data');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
end
