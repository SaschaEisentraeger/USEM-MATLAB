function R = ms_calibrate(elems, P, C, form, withMass, tol)
%MS_CALIBRATE Number of Gauss points per direction for which the element
%matrices of a set of elements are converged to round-off.
%
%   R = ms_calibrate(elems, P, C, form, withMass, tol)
%
%   elems    : cell array {x, y; ...} of element node coordinates
%   P        : polynomial degree structure (pxi = peta = p)
%   form     : 'PP', 'PM' or 'PH' (ms_Ke)
%   withMass : include the consistent mass matrix (ms_Me)
%   tol      : relative tolerance (default 1e-12)
%
%   As in the companion paper (axisymmetric solids), the rule is the
%   smallest nG >= p+1 for which the matrices with nG AND with nG+1 points
%   agree with a reference to a relative deviation below tol in the
%   Frobenius norm, maximized over the elements. The reference is computed
%   with nRef = max(64, 8p) points (the companion paper uses max(40, 8p);
%   on the nearly degenerate elements of this study, 40 points are not
%   sufficient). The round-off level of the element matrices is estimated
%   by the deviation of a second reference with nRef + 8 points; if it is
%   below 1e-8, i.e., if the integrals converge, the tolerance is raised to
%   max(tol, 10 x round-off level), since the metric functions of high
%   degree carry round-off errors of the order of cond(V) u (Sect. 3). The smallest nG is found by an exponential
%   search followed by bisection (assuming monotonic convergence); the
%   rule nG+1 is then checked explicitly. If no rule below nRef satisfies
%   the criterion (singular integrands on folded elements), R.conv = false
%   and R.nG = nRef.
%
%   R.nG, R.conv, R.dev2p (deviation of the former default rule, 2p points),
%   R.devG (deviation at nG), R.nRef, R.tol (tolerance used), R.floor
%   (round-off estimate), R.rejected (trial space rejected).

if nargin < 5, withMass = false; end
if nargin < 6, tol = 1e-12; end
p = P.pxi;
nRef = max(64, 8*p);
R = struct('nG', NaN, 'conv', false, 'dev2p', NaN, 'devG', NaN, 'nRef', nRef, ...
           'tol', tol, 'floor', NaN, 'rejected', false);
ne = size(elems, 1);
Kref = cell(ne, 1);  Mref = cell(ne, 1);
try
    for e = 1:ne
        Kref{e} = ms_Ke(P, C, elems{e,1}, elems{e,2}, nRef, form);
        if withMass, Mref{e} = ms_Me(P, elems{e,1}, elems{e,2}, nRef, form); end
    end
catch err
    if startsWith(err.identifier, 'ShapeFunctions2D_Metric:')
        R.rejected = true;
        return
    end
    rethrow(err);
end
cache = containers.Map('KeyType', 'double', 'ValueType', 'double');
    function d = dev(nq)
        if isKey(cache, nq), d = cache(nq); return, end
        d = 0;
        for k = 1:ne
            K = ms_Ke(P, C, elems{k,1}, elems{k,2}, nq, form);
            d = max(d, norm(K - Kref{k}, 'fro')/norm(Kref{k}, 'fro'));
            if withMass
                M = ms_Me(P, elems{k,1}, elems{k,2}, nq, form);
                d = max(d, norm(M - Mref{k}, 'fro')/norm(Mref{k}, 'fro'));
            end
        end
        cache(nq) = d;
    end
% round-off level: deviation of a second reference rule
R.floor = dev(nRef + 8);
if R.floor < 1e-8
    tol = max(tol, 10*R.floor);
end
R.tol = tol;
ok = @(nq) dev(nq) < tol && dev(nq + 1) < tol;
R.dev2p = dev(2*p);
% exponential search for a rule below nRef that satisfies the criterion
lo = p;  hi = NaN;  step = 1;  nq = p + 1;
while nq < nRef
    if ok(nq), hi = nq; break, end
    lo = nq;  nq = min(nq + step, nRef - 1);  step = 2*step;
    if lo == nRef - 1, break, end
end
if isnan(hi)
    R.nG = nRef;  R.conv = false;  R.devG = 0;
    return
end
% bisection between lo (fails) and hi (satisfies)
while hi - lo > 1
    mid = floor((lo + hi)/2);
    if ok(mid), hi = mid; else, lo = mid; end
end
R.nG = hi;  R.conv = true;  R.devG = dev(hi);
end
