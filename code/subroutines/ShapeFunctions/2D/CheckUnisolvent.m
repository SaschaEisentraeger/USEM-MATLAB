function CheckUnisolvent(V,px,py)
%CHECKUNISOLVENT Stop if the metric interpolation matrix is singular.
%
% Unisolvence check
% The metric shape functions are undefined if the nodes are not unisolvent
% for the polynomial space (e.g. coincident nodes or strongly distorted
% elements): stop instead of returning meaningless shape functions

rc = rcond(V);
if ~(rc >= eps)
    error('ShapeFunctions2D_Metric:singular', ...
         ['The interpolation matrix of the metric shape functions is ' ...
          'singular to working precision (rcond = %.2e): the element ' ...
          'nodes are (nearly) not unisolvent. The element is too ' ...
          'strongly distorted for the metric shape functions of degree ' ...
          '(%d,%d).'], rc, px, py);
end
end
