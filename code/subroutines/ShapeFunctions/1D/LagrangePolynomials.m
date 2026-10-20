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
        PHIp(i) = PHIp_fun(i,xii');
    end

    % Computing barycentric weights
    c = 1./PHIp;

    % Allocation of solution vector
    N  = cell(n+1,1);
    Np = cell(n+1,1);
    for i = 1:n+1
        N{i} = @(xi)l(xi,i,c,xii');
        Np{i} = @(xi)lp(xi,i,c,xii');
    end
end

%% For explanation of the subsequent functions look in the diploma thesis of D. Schmicker.
function ret = PHI(xi,nodes)
    ret = ones(size(xi));
    for i = 1:size(nodes,2)
        ret = ret.*(xi-nodes(i));
    end
end

function ret = PHIp_fun(node,nodes)
    ret = 1;
    for i = 1:size(nodes,2)
        if node~=i
            ret = ret*(nodes(node)-nodes(i));
        end
    end 
end

% l and lp are evaluated for all points xi at once; the operations are the
% same (and in the same order) for every point as in a pointwise loop
function ret = l(xi,i,c,nodes)
    ret = c(i)*PHI(xi,nodes)./(xi-nodes(i));
    ret(xi==nodes(i)) = 1;
end

function ret = lp(xi,i,c,nodes)
    % sum_(k~=i) 1/(xi-x_k)
    s = zeros(size(xi));
    for k = 1:max(size(nodes))
        if k~=i
            s = s+1./(xi-nodes(k));
        end
    end
    ret = s.*l(xi,i,c,nodes);
    % xi at node i
    ret(xi==nodes(i)) = s(xi==nodes(i));
    % xi at another node z: c(i)/(xi-x_i)/c(z)
    for z = 1:max(size(nodes))
        if z~=i
            at = xi==nodes(z);
            ret(at) = c(i)./(xi(at)-nodes(i))/c(z);
        end
    end
end
