%% Inverse mapping from the global to the local domain (Newton's method)
%-------------------------------------------------------------------------%
%
% Sascha Eisenträger, Otto von Guericke University, 20.10.2026
%
% Local coordinates (xi,eta) of global points (x,y) inside an element,
% x(xi,eta) = sum_k N_k(xi,eta) x_k, by Newton's method, vectorized over
% the points. Unlike Global2LocalCoordinatesQuad, no toolbox is required
% and the result is not rounded.
%
%-------------------------------------------------------------------------%
%% Input
%
% (xQ,yQ):      Global coordinates of the point(s)
% (xN,yN):      Global coordinates of the element nodes
% p:            Polynomial degree template (pxi, peta, ndist)
%-------------------------------------------------------------------------%
%% Output
%
% xiQ:          Local coordinates [xi, eta] of the point(s), nQ x 2
%-------------------------------------------------------------------------%
function xiQ = Global2LocalNewton(xQ,yQ,xN,yN,p)

xQ = xQ(:);
yQ = yQ(:);
nQ = numel(xQ);

% Tolerance relative to the element size
h = max(max(xN) - min(xN), max(yN) - min(yN));
tol = 1e-14*h;

xi  = zeros(nQ,1);
eta = zeros(nQ,1);
converged = false;
for it = 1:50
    [N,dNxi,dNeta] = ShapeFunctions2D_SEM(xi,eta,p);
    rx = N'*xN - xQ;
    ry = N'*yN - yQ;
    if max(abs([rx; ry])) <= tol
        converged = true;
        break
    end
    % Jacobian d(x,y)/d(xi,eta) at every point
    xxi  = dNxi'*xN;
    yxi  = dNxi'*yN;
    xeta = dNeta'*xN;
    yeta = dNeta'*yN;
    detJ = xxi.*yeta - xeta.*yxi;
    xi  = xi  - ( yeta.*rx - xeta.*ry)./detJ;
    eta = eta - (-yxi.*rx  + xxi.*ry )./detJ;
end
if ~converged
    warning('Global2LocalNewton:noConvergence', ...
            'Inverse mapping did not converge (residual %.2e).', ...
            max(abs([rx; ry])));
end

xiQ = [xi, eta];

end
