function ret = LagrangeBarycentricDerivative(xi,i,c,nodes)
%LAGRANGEBARYCENTRICDERIVATIVE Derivative of the barycentric Lagrange polynomial i.

    % sum_(k~=i) 1/(xi-x_k)
    s = zeros(size(xi));
    for k = 1:max(size(nodes))
        if k~=i
            s = s+1./(xi-nodes(k));
        end
    end
    ret = s.*LagrangeBarycentric(xi,i,c,nodes);
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
