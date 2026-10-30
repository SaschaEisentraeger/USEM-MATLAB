function root = cyl_setup_paths()
%CYL_SETUP_PATHS Search path for the thick-cylinder study.
%   Adds the code base (code/subroutines) and the test helpers
%   (tests/helpers). Returns the repository root.

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..', '..');
addpath(genpath(fullfile(root, 'code', 'subroutines')));
addpath(fullfile(root, 'tests', 'helpers'));
addpath(here);
end
