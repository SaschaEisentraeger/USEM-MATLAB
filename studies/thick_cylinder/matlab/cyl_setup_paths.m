function root = cyl_setup_paths()
%CYL_SETUP_PATHS Search path for the thick-cylinder study.
%   Adds the code base (code/subroutines), the test helpers (tests/helpers)
%   and, in GNU Octave, the symbolic package and the compatibility folder.
%   Returns the repository root.

here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..', '..');
addpath(genpath(fullfile(root, 'code', 'subroutines')));
addpath(fullfile(root, 'tests', 'helpers'));
addpath(here);
if exist('OCTAVE_VERSION', 'builtin')
    % Octave's unique under another name, called by octave_compat/unique.m
    orig = fullfile(tempdir, 'usem_octave_compat');
    if ~exist(fullfile(orig, 'octave_unique_orig.m'), 'file')
        if ~exist(orig, 'dir'), mkdir(orig); end
        src = fileread(which('unique'));
        src = regexprep(src, 'function\s+(\[[^\]]*\]\s*=\s*)?unique\s*\(', ...
                        'function $1octave_unique_orig (', 'once');
        fid = fopen(fullfile(orig, 'octave_unique_orig.m'), 'w');
        fprintf(fid, '%s', src);
        fclose(fid);
    end
    addpath(orig);
    addpath(fullfile(here, 'octave_compat'));
    pkg load symbolic
    warning('off', 'all');
end
end
