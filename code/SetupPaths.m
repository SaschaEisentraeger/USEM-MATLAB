function SetupPaths()
%SETUPPATHS Add the USEM code to the MATLAB path.
%
%   SetupPaths() adds code/subroutines (with all subfolders) and the
%   example input files code/InputFilesAbaqus to the MATLAB path. Call it
%   once per session before using the routines outside Start_Sim:
%
%       addpath('path/to/USEM/code'); SetupPaths

codeDir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(codeDir,'subroutines')), ...
        genpath(fullfile(codeDir,'InputFilesAbaqus')));

end
