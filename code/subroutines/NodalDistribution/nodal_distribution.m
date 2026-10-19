function [x] = nodal_distribution(p,distribution)

% returns the nodal distribution that is used to generate the shape
% functions for the spectral element method (SEM) --> supporting points for
% the Lagrange interpolation polynomials

%--------------------------------------------------------------------------
% Equidistant points
if distribution == 1
    x = equispaced_nodes(p);
%--------------------------------------------------------------------------
% Gauﬂ-Lobatto-Legendre points
elseif distribution == 2
    x = lglnodes(p);
%--------------------------------------------------------------------------
% Gauﬂ-Lobatto-Chebyshev points
elseif distribution == 3
    x = cgl_nodes(p);
else
    error('nodal_distribution:ndist', ...
          'Unknown nodal distribution %d (1: equidistant, 2: GLL, 3: GLC).', ...
          distribution);
end
%--------------------------------------------------------------------------

x = -x;

end