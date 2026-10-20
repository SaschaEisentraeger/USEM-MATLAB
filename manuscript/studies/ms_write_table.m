function ms_write_table(file, T)
%MS_WRITE_TABLE Struct of equally long column vectors -> whitespace table
%(header line with the field names), read by pgfplots.
names = fieldnames(T);
n = numel(T.(names{1}));
fid = fopen(file, 'w');
fprintf(fid, '%s ', names{:});
fprintf(fid, '\n');
for i = 1:n
    for k = 1:numel(names)
        fprintf(fid, '%.6e ', T.(names{k})(i));
    end
    fprintf(fid, '\n');
end
fclose(fid);
end
