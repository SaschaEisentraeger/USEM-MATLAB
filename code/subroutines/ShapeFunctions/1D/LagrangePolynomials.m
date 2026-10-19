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

function ret = l(xi,i,c,nodes)
    ret = zeros(size(xi));
    for j = 1:max(size(xi))
        if xi(j)==nodes(i)
            ret(j) = 1;
        else
            ret(j) = c(i)*PHI(xi(j),nodes)./(xi(j)-nodes(i));
        end
    end
end

function ret = lp(xi,i,c,nodes)
    ret = zeros(size(xi));
    for j = 1:max(size(xi)) 
        if xi(j)==nodes(i)
            ret(j) = 0;
            for k = 1:max(size(nodes))
               if k~=i
                  ret(j) = ret(j)+1/(xi(j)-nodes(k)); 
               end
            end
        elseif ~isempty(find(nodes==xi(j), 1))
            z = nodes==xi(j);
            ret(j) = c(i)/(xi(j)-nodes(i))/c(z);
        else
            ret(j) = 0;
            for k = 1:max(size(nodes))
               if k~=i
                  ret(j) = ret(j)+1/(xi(j)-nodes(k)); 
               end
            end
            ret(j) = ret(j)*l(xi(j),i,c,nodes);
        end
    end
end