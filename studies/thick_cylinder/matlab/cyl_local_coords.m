function [xi, eta] = cyl_local_coords(p)
%CYL_LOCAL_COORDS Reference coordinates of the nodes of a spectral element.
%   Node order of ShapeFunctions2D_SEM / pSpectralElement_NodesTopology:
%   corners (-1,-1), (1,-1), (1,1), (-1,1); edge 1 (eta = -1), edge 2
%   (xi = 1), edge 3 (eta = 1), edge 4 (xi = -1), each with increasing
%   coordinate; interior nodes row by row (xi running fastest).
t = nodal_distribution(p.pxi);  t = t(:);
s = nodal_distribution(p.peta); s = s(:);
ti = t(2:end-1);  si = s(2:end-1);
[XI, ETA] = meshgrid(ti, si);  XI = XI.';  ETA = ETA.';
xi  = [-1; 1; 1; -1; ti; ones(size(si)); ti; -ones(size(si)); XI(:)];
eta = [-1; -1; 1; 1; -ones(size(ti)); si; ones(size(ti)); si; ETA(:)];
end
