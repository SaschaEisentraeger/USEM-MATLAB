function [N,Np] = LagrangePolynomials(xii)
    %
    % This function computes the Lagrange-Interpolation functions and its
    % derivatives on an equispaced grid on the domain [-1,1]. The output format
    % utilizes the Matlab anonymous function feature.
    %
    % Inputs: 
    %               xii - Node-distribution [n x 1]
    %
    % Outputs:  
    %               N   - Interpolation functions [n x 1]
    %               Np  - Derivatives of interpolation functions [n x 1]
    %
    %% Numeric calculation by means of barycentric weights

    % Polynomial degree
    n = size(xii,1)-1;

    % Computing values of PHI'
    PHIp = zeros(n+1,1);
    for i = 1:n+1
        PHIp(i) = LagrangeNodeDerivative(i,xii');
    end

    % Computing barycentric weights
    c = 1./PHIp;

    % Allocation of solution vector
    N  = cell(n+1,1);
    Np = cell(n+1,1);
    for i = 1:n+1
        N{i} = @(xi)LagrangeBarycentric(xi,i,c,xii');
        Np{i} = @(xi)LagrangeBarycentricDerivative(xi,i,c,xii');
    end
end
