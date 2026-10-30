function varargout = ExpandToAllDoFs(NumDoF,ActiveDoFs,varargin)
%EXPANDTOALLDOFS Histories of all DOFs from those of the active DOFs.
%
%   [U, V, ...] = ExpandToAllDoFs(NumDoF, ActiveDoFs, uvec, vvec, ...)
%   returns the histories [nt x NumDoF] with zeros at the constrained DOFs
%   if the histories uvec, vvec, ... [nt x nSaved] contain every active
%   DOF (nSaved = numel(ActiveDoFs)). Otherwise (only some DOFs saved)
%   they are returned unchanged.

varargout = varargin;
for k = 1:numel(varargin)
    h = varargin{k};
    if size(h,2) == numel(ActiveDoFs)
        X = zeros(size(h,1),NumDoF);
        X(:,ActiveDoFs) = h;
        varargout{k} = X;
    end
end

end
