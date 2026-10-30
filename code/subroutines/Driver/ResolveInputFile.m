function file = ResolveInputFile(filename)
%RESOLVEINPUTFILE Full path of an input file given as a path or found on the MATLAB path.

% Full path of the input file (given as a path or found on the MATLAB path)
if isfile(filename)
    file = filename;
    if isempty(fileparts(file))
        file = fullfile(pwd,file);
    end
else
    file = which(filename);
    if isempty(file)
        error('ReadModel:fileNotFound', ...
              'Input file "%s" not found (neither a file nor on the path).', ...
              filename);
    end
end
end
