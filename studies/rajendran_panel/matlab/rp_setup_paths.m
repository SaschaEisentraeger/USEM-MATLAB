function rp_setup_paths()
%RP_SETUP_PATHS Search path for the panel study (code base, test helpers,
%   helpers of the thick-cylinder study).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'thick_cylinder', 'matlab'));
cyl_setup_paths();
addpath(here);
end
