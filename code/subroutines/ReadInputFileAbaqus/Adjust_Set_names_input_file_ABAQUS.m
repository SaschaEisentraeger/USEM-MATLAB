function [] = Adjust_Set_names_input_file_ABAQUS(inp_file,out_file)
%ADJUST_SET_NAMES_INPUT_FILE_ABAQUS Make the set names of an input file unique.
%
%   Adjust_Set_names_input_file_ABAQUS(inp_file) writes a copy of the
%   ABAQUS input file, <name>_adjusted.inp next to it (or in the current
%   folder if inp_file has no path), in which the set names are unique:
%   'Set-' after 'End Instance' becomes 'PickedSet-', and the sets named
%   'Traction'/'Pressure' are numbered consecutively. Since ABAQUS 6.12 a
%   set name can occur in several parts of the file, while the readers
%   (Read_*_ABAQUS) identify sets by their names.
%
%   Adjust_Set_names_input_file_ABAQUS(inp_file, out_file) writes the copy
%   to out_file instead (used by ReadModel with a temporary file).

% Read the content of the file into a string
Cont = fileread(inp_file);

% Find the keyword 'End Instance'
IdxEndInstance = strfind(Cont, 'End Instance');

% Replace the keyword 'Set-' by 'PickedSet-'
ContAux= strrep(Cont(IdxEndInstance:end),'Set-','PickedSet-');
ContNew = [Cont(1:IdxEndInstance-1) ContAux];

% Find the keywords 'Traction' or 'Pressure' and number those sets
% consecutively. Note that there is always an element set (elset) and a
% corresponding nodal set (nset).
idxTrac = strfind(ContNew, '=Traction');
idxPres = strfind(ContNew, '=Pressure');
nTrac = numel(idxTrac)/2;
nPres = numel(idxPres)/2;

for i = 1:nTrac
    if i == 1
        ContAux = strrep(ContNew(1:idxTrac(2*i)+10),'=Traction',...
                  ['=Traction' num2str(i)]);
        ContNew = [ContAux ContNew(idxTrac(2*i)+11:end)];
    else
        ContAux = strrep(ContNew(1+idxTrac(2*(i-1))+10:idxTrac(2*i)+10),...
                  '=Traction',['=Traction' num2str(i)]);
        ContNew = [ContNew(1:idxTrac(2*(i-1))+10) ContAux ...
                   ContNew(idxTrac(2*i)+11:end)];
    end
    idxTrac = idxTrac + 2;
end

for i = 1:nPres
    if i == 1
        ContAux = strrep(ContNew(1:idxPres(2*i)+10),'=Pressure',...
                  ['=Pressure' num2str(i)]);
        ContNew = [ContAux ContNew(idxPres(2*i)+11:end)];
    else
        ContAux = strrep(ContNew(1+idxPres(2*(i-1))+10:idxPres(2*i)+10),...
                  '=Pressure',['=Pressure' num2str(i)]);
        ContNew = [ContNew(1:idxPres(2*(i-1))+10) ContAux ...
                   ContNew(idxPres(2*i)+11:end)];
    end
    idxPres = idxPres + 2;
end

% Generate output file
if nargin < 2
    [pathstr, name, ext] = fileparts(inp_file);
    out_file = fullfile(pathstr,[name '_adjusted' ext]);
end
fout = fopen(out_file,'w');
if fout < 0
    error('Adjust_Set_names_input_file_ABAQUS:write', ...
          'Cannot write "%s".', out_file);
end

% '%s': the file content must not be interpreted as a format string
fprintf(fout,'%s',ContNew);

% Close the output file (only this one: other open files are not touched)
fclose(fout);

end