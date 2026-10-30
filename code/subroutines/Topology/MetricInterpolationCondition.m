function rc = MetricInterpolationCondition(xN,yN,p)
%METRICINTERPOLATIONCONDITION Reciprocal condition number of the metric
%interpolation matrix of one element.
%
%   rc = MetricInterpolationCondition(xN, yN, p) returns rcond(V) of the
%   generalized Vandermonde matrix V of the metric (trial) shape functions
%   for the element nodes (xN, yN) and the settings in p (degree, basis,
%   trial space). rc = 0 if V is singular to working precision, e.g. for
%   coincident nodes; rc close to 1 for well-shaped elements with the
%   default Lagrange basis.

try
    % Query point: node 1 (local coordinates (-1,-1), used by the hybrid
    % space instead of an inverse mapping)
    [~,~,~,V] = ShapeFunctions2D_Metric(xN(1),yN(1),xN(:),yN(:),p,[],[-1 -1]);
    rc = rcond(V);
catch err
    if strcmp(err.identifier,'ShapeFunctions2D_Metric:singular')
        rc = 0;
    else
        rethrow(err);
    end
end
if ~(rc >= 0)
    rc = 0;     % NaN (e.g. from coincident nodes)
end

end
