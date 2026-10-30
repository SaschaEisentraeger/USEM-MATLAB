function ret = LagrangeNodePolynomial(xi,nodes)
%LAGRANGENODEPOLYNOMIAL Product of (xi - x_k) over all nodes x_k.
%
% For explanation of the subsequent functions look in the diploma thesis of D. Schmicker.

    ret = ones(size(xi));
    for i = 1:size(nodes,2)
        ret = ret.*(xi-nodes(i));
    end
end
