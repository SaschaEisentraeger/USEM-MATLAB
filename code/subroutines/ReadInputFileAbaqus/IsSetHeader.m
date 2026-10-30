function tf = IsSetHeader(line, header)
%ISSETHEADER True if a line of an input file is the header of a given set.

% True if line is the header of the set definition "header", e.g.
% '*Nset, nset=Set-1' matches '*Nset, nset=Set-1, instance=Part-1' but not
% '*Nset, nset=Set-10'.
if ~ischar(line)
    tf = false;
    return
end
line = strtrim(line);
n = numel(header);
tf = strncmp(line, header, n) && (numel(line) == n || line(n+1) == ',');
end
