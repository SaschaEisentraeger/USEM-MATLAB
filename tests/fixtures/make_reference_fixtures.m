function make_reference_fixtures(resultsDir)
%MAKE_REFERENCE_FIXTURES Extract regression data from saved workspaces.
%
%   make_reference_fixtures(resultsDir) reads the workspaces that the
%   former script call_Static.m saved to code/subroutines/Results/Static
%   (Lee & Bathe patch tests, see docs/chapters/input-files.tex; loose
%   variables, whereas SaveResults now stores the structs model, settings
%   and results) and stores the variables needed to repeat the analysis
%   in tests/fixtures/reference/<case>.mat.
%
%   The fixtures are the input data AFTER the ABAQUS input file has been
%   read and the boundary conditions have been adapted to the spectral
%   elements, together with the reference solution u. They allow the
%   regression tests to run without the (not distributed) input files.
%
%   The fixtures in the repository were generated once from the results
%   provided with the code (static analyses, plane stress, GLL nodes).

if nargin < 1
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    resultsDir = fullfile(root, 'code', 'subroutines', 'Results', 'Static');
end
outDir = fullfile(fileparts(mfilename('fullpath')), 'reference');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

vars = {'filename', 'p', 'p_org', 'antype', 'nodes', 'topology', ...
        'Topology_Node', 'xCoor', 'yCoor', 'C', 'rho', 'zeta', ...
        'material_elset', 'material_properties', 'mass_matrix', ...
        'cload', 'dsload', 'body_force', 'gravity', ...
        'boundary_conditions', 'BC_DOFs', 'ActiveDoF', 'num_DoF', ...
        'num_Nodes', 'num_Elem', 'F', 'u'};

files = dir(fullfile(resultsDir, '**', '*.mat'));
for k = 1:numel(files)
    S = load(fullfile(files(k).folder, files(k).name), vars{:}); %#ok<NASGU>
    [~, name] = fileparts(files(k).name);
    save(fullfile(outDir, [name '.mat']), '-struct', 'S', '-v7');
    fprintf('%s\n', name);
end

end
