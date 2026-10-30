function x = nodal_distribution(p)
%NODAL_DISTRIBUTION Gauss-Lobatto-Legendre nodes on [-1,1] in ascending order.
%
%   x = nodal_distribution(p) returns the p+1 supporting points of the
%   Lagrange interpolation polynomials of the spectral element method.

x = -lglnodes(p);

end
