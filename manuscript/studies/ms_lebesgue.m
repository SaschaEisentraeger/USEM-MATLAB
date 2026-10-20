function [LN, LQ, LT] = ms_lebesgue(x, y, P, nGrid)
%MS_LEBESGUE Lebesgue constants of one element on an nGrid x nGrid grid
%of reference points: LN parametric (GLL) interpolation, LQ metric
%interpolation in Q_p(s,t) (method 1), LT hybrid space T_p (method 2).
if nargin < 4, nGrid = 6*max(P.pxi, P.peta) + 11; end
[Xi, Eta] = meshgrid(linspace(-1, 1, nGrid));
IP = [Xi(:), Eta(:)];
Xg = Local2GlobalCoordinatesQuad(P, IP, [x, y]);
LN = max(sum(abs(ShapeFunctions2D_SEM(IP(:,1), IP(:,2), P)), 1));
P.metricBasis = 'lagrange';
M = ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, P, 1);
LQ = max(sum(abs(M), 1));
M = ShapeFunctions2D_Metric(Xg(:,1), Xg(:,2), x, y, P, 2, IP);
LT = max(sum(abs(M), 1));
end
