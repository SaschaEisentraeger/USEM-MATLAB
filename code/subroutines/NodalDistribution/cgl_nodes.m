function [x] = cgl_nodes(N)
%-------------------------------------------------------------------------%
%
% cgl_nodes.m
%
% Patera "A Spectral Element Method for Fluid Dynamics: Laminar Flow in a 
%         Channel Expansion"
% Journal of Computational Physics 1984, 54, 468-488
%
% Computes the Chebyshev-Gauss-Lobatto nodes
% N:		Polynomial degree
%
%-------------------------------------------------------------------------%

% Chebyshev-Gauss-Lobatto nodes 
x = cos(pi*(0:N)/N)';
% Set values that are less than the machine precision to 0
x(abs(x) < eps) = 0;
% Symmetric points
% round N/2 to the next larger integer
n1 = ceil((N+1)/2);
% round N/2 to the next smaller integer
n2 = floor((N+1)/2);
% Copy the first half the nodes and reverse their order
x(n1+1:N+1) = -flip(x(1:n2));

end