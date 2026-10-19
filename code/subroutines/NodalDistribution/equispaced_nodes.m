function [x] = equispaced_nodes(N)
%-------------------------------------------------------------------------%
%
% equispaced_nodes.m
%
% Returns equispaced nodes for a polynomial of order N (i.e. N+1 nodes are
% needed)
%
% N:		Polynomial degree
%
%-------------------------------------------------------------------------%

% N1: number of nodes (polynomial degree + 1)
N1 = N + 1;

% equispaced nodes 
x = linspace(1,-1,N1)';

end