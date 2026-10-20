function err = cyl_errors(model, u, uEx, prob, C, form, nq, ref)
%CYL_ERRORS Relative errors of a spectral element solution.
%
%   Displacement and strain field of the solution:
%     'PM' : metric trial functions (ShapeFunctions2D_Metric), i.e. the
%            trial space of the unsymmetric element ('MM': the same)
%     'PP' : parametric shape functions (conventional element)
%   Integrals over the elements use the signed det J of the geometry
%   mapping (as in the element formulation). ref = struct with the exact
%   energy (int eps' C eps dA) and L2 norm squared (int |u|^2 dA).
%
%   err.energy : relative energy norm of the strain error
%   err.L2     : relative L2 norm of the displacement error
%   err.nodal  : max nodal |u_h - u_ex| / max nodal |u_ex|
p = model.p;
if nargin < 7 || isempty(nq), nq = 3*p.pxi + 2; end
% metric functions: corotational skew coordinates (code-base default) unless
% the model requests Cartesian coordinates (model.metricMethod = 0)
metric = 1;
if isfield(model, 'metricMethod'), metric = model.metricMethod; end
T = model.Topology_Node;
Qd = cyl_quad(p, nq);
W = Qd.W;  N = Qd.N;  dNxi = Qd.dNxi;  dNeta = Qd.dNeta;
eE = 0;  eL = 0;
for e = 1:size(T, 1)
    ids = T(e, :).';
    x = model.xCoor(ids);  y = model.yCoor(ids);
    ux = u(2*ids-1);  uy = u(2*ids);
    [J11, J12, J21, J22, detJ] = cyl_jacobian(dNxi, dNeta, x, y);
    xg = N.'*x;  yg = N.'*y;
    if any(strcmp(form, {'PM', 'MM'}))
        [M, dMx, dMy] = ShapeFunctions2D_Metric(xg, yg, x, y, p, metric);
    else
        M = N;
        dMx = ( dNxi.*J22(:).' - dNeta.*J12(:).')./detJ(:).';
        dMy = (-dNxi.*J21(:).' + dNeta.*J11(:).')./detJ(:).';
    end
    duxh = prob.ux(xg, yg) - M.'*ux;
    duyh = prob.uy(xg, yg) - M.'*uy;
    d = [prob.exx(xg, yg) - dMx.'*ux, ...
         prob.eyy(xg, yg) - dMy.'*uy, ...
         prob.gxy(xg, yg) - (dMy.'*ux + dMx.'*uy)];
    w = W(:).*detJ(:);
    eE = eE + sum(sum((d*C).*d, 2).*w);
    eL = eL + sum((duxh.^2 + duyh.^2).*w);
end
err.energy = sqrt(max(eE, 0)/ref.energy);
err.L2 = sqrt(max(eL, 0)/ref.L2);
err.energySigned = eE/ref.energy;     % negative values flag signed-measure effects
dn = hypot(u(1:2:end) - uEx(1:2:end), u(2:2:end) - uEx(2:2:end));
err.nodal = max(dn)/max(hypot(uEx(1:2:end), uEx(2:2:end)));
end
