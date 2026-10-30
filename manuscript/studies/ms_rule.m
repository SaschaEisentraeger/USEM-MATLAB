function [nG, conv] = ms_rule(family, form, p)
%MS_RULE Calibrated number of Gauss points per direction of the element
%matrices (ms_quadrature_calibration) for a mesh family, a formulation
%('PP', 'PM', 'PH') and a degree p. conv = false marks rules that do not
%converge (folded elements); nG is then the reference rule max(64, 8p).
persistent tables
if isempty(tables), tables = containers.Map(); end
if ~isKey(tables, family)
    [~, dataDir] = ms_setup();
    file = fullfile(dataDir, 'quadrature', [family '.dat']);
    if ~exist(file, 'file')
        error('ms_rule:missing', 'No calibrated rule for "%s" (run ms_quadrature_calibration).', family);
    end
    fid = fopen(file);
    head = strsplit(strtrim(fgetl(fid)));
    vals = fscanf(fid, '%f', [numel(head), Inf]).';
    fclose(fid);
    tables(family) = struct('head', {head}, 'vals', vals);
end
T = tables(family);
row = T.vals(:, strcmp(T.head, 'p')) == p;
if ~any(row)
    error('ms_rule:degree', 'No calibrated rule for "%s", p = %d.', family, p);
end
nG = T.vals(row, strcmp(T.head, ['nG_' form]));
conv = T.vals(row, strcmp(T.head, ['conv_' form])) == 1;
end
