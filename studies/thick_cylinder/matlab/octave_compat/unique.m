function varargout = unique(varargin)
%UNIQUE Octave stand-in: Octave 8 does not return the third output (J) of
%   unique(..., 'stable'), which BC_Edge_sub uses. All work is done by
%   Octave's own unique (made callable as octave_unique_orig by
%   cyl_setup_paths); only a missing J is computed, such that C(J) = A
%   (C(J,:) = A for 'rows').
varargout = cell(1, max(nargout, 1));
[varargout{:}] = octave_unique_orig(varargin{:});
if nargout >= 3 && isempty(varargout{3}) && ~isempty(varargin{1})
    A = varargin{1};  C = varargout{1};
    if any(strcmp(varargin(2:end), 'rows'))
        [~, J] = ismember(A, C, 'rows');
    else
        [~, J] = ismember(A(:), C(:));
    end
    varargout{3} = J(:);
end
end
