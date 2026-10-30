function [root, dataDir] = ms_setup()
%MS_SETUP Paths for the studies of the manuscript (MATLAB).
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(genpath(fullfile(root, 'code', 'subroutines')));
addpath(fullfile(root, 'tests', 'helpers'));
% vectorized element routines (cyl_Ke_PM, cyl_Ke_PP, cyl_quad), used by ms_Ke
addpath(fullfile(root, 'studies', 'thick_cylinder', 'matlab'));
dataDir = fullfile(fileparts(here), 'data');
if ~exist(dataDir, 'dir'), mkdir(dataDir); end
end
