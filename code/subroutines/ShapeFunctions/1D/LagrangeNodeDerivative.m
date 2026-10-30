function ret = LagrangeNodeDerivative(node,nodes)
%LAGRANGENODEDERIVATIVE Derivative of the node polynomial at a node, prod_(k~=i) (x_i - x_k).

    ret = 1;
    for i = 1:size(nodes,2)
        if node~=i
            ret = ret*(nodes(node)-nodes(i));
        end
    end 
end
