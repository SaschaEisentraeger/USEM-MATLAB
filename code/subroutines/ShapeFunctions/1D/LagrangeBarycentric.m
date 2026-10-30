function ret = LagrangeBarycentric(xi,i,c,nodes)
%LAGRANGEBARYCENTRIC Lagrange polynomial i in barycentric form.
%
% LagrangeBarycentric and LagrangeBarycentricDerivative are evaluated for all points xi at once; the operations are the
% same (and in the same order) for every point as in a pointwise loop

    ret = c(i)*LagrangeNodePolynomial(xi,nodes)./(xi-nodes(i));
    ret(xi==nodes(i)) = 1;
end
