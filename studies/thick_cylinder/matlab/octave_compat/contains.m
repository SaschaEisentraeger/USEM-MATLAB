function tf = contains(str, pattern, varargin)
%CONTAINS Octave stand-in for MATLAB's contains (char/cellstr input).
ignoreCase = numel(varargin) >= 2 && strcmpi(varargin{1}, 'IgnoreCase') && varargin{2};
if ischar(pattern), pattern = {pattern}; end
if ischar(str)
    tf = any(cellfun(@(q) hit(str, q, ignoreCase), pattern));
else
    tf = cellfun(@(s) any(cellfun(@(q) hit(s, q, ignoreCase), pattern)), str);
end
end

function h = hit(s, q, ic)
if ic, s = lower(s); q = lower(q); end
h = ~isempty(strfind(s, q));
end
