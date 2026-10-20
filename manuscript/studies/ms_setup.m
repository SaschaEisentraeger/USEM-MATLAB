function [root, dataDir] = ms_setup()
%MS_SETUP Paths for the studies of the manuscript (MATLAB or GNU Octave).
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(fullfile(root, 'code', 'subroutines')));
addpath(fullfile(root, 'tests', 'helpers'));
if exist('OCTAVE_VERSION', 'builtin')
    addpath(fullfile(root, 'studies', 'metric_basis', 'octave'));
    warning('off', 'all');
end
dataDir = fullfile(fileparts(here), 'data');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
end
