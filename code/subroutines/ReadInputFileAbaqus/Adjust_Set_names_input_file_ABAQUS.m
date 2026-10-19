function [] = Adjust_Set_names_input_file_ABAQUS(inp_file)

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
[pathstr, name, ext] = fileparts(inp_file);
if isempty(pathstr)
    fout = fopen([name '_adjusted' ext],'w');
else
    fout = fopen(fullfile(pathstr,[name '_adjusted' ext]),'w');
end

% '%s': the file content must not be interpreted as a format string
fprintf(fout,'%s',ContNew);

% Close File
fclose('all');

end