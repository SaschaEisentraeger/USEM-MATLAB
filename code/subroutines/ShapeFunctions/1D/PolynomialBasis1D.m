%% One-dimensional polynomial bases for the metric shape functions
%-------------------------------------------------------------------------%
%
% Sascha Eisenträger, Otto von Guericke University, 19.10.2026
%
% Evaluates a basis {phi_0, ..., phi_n} of the polynomials of degree n and
% its first derivative. The bases span the same space, they differ only in
% the conditioning of the generalized Vandermonde matrix that is built from
% them in ShapeFunctions2D_Metric.
%
% Bases:
% 'monomial':   phi_k(s) = s^k
% 'legendre':   phi_k(s) = P_k(s), Legendre polynomials (three-term
%               recurrence)
% 'chebyshev':  phi_k(s) = T_k(s), Chebyshev polynomials of the first kind
% 'lagrange':   phi_k(s) = l_k(s), Lagrange polynomials through the nodes
%               of the reference element (nodal_distribution(n,ndist)),
%               evaluated with the product formula. In skew coordinates the
%               Vandermonde matrix of an affine element is then the
%               identity.
%
% Reference:
% [1] Gautschi, W. / Inglese, G.:
% Lower Bounds for the Condition Number of Vandermonde Matrices, 1988,
% Numerische Mathematik, Vol. 52, No. 3, Springer, p. 241-250
%
%-------------------------------------------------------------------------%
%% Input
%
% s:        Evaluation points (vector)
% n:        Polynomial degree
% basis:    'monomial', 'legendre', 'chebyshev' or 'lagrange'
% ndist:    Nodal distribution of the reference element (only 'lagrange')
%-------------------------------------------------------------------------%
%% Output
%
% Phi:      Basis functions at s, (n+1) x numel(s)
% dPhi:     Derivatives of the basis functions at s, (n+1) x numel(s)
%-------------------------------------------------------------------------%
function [Phi,dPhi] = PolynomialBasis1D(s,n,basis,ndist)

% Evaluation points as a row vector
s = reshape(s,1,[]);
nQ = numel(s);

Phi  = zeros(n+1,nQ);
dPhi = zeros(n+1,nQ);

switch lower(basis)
    case 'monomial'
        Phi(1,:) = 1;
        for k = 1:n
            Phi(k+1,:)  = s.^k;
            dPhi(k+1,:) = k*s.^(k-1);
        end

    case 'legendre'
        % (k+1) P_{k+1} = (2k+1) s P_k - k P_{k-1}
        % P'_{k+1} = P'_{k-1} + (2k+1) P_k
        Phi(1,:) = 1;
        if n >= 1
            Phi(2,:)  = s;
            dPhi(2,:) = 1;
        end
        for k = 1:n-1
            Phi(k+2,:)  = ((2*k+1)*s.*Phi(k+1,:) - k*Phi(k,:))/(k+1);
            dPhi(k+2,:) = dPhi(k,:) + (2*k+1)*Phi(k+1,:);
        end

    case 'chebyshev'
        % T_{k+1} = 2 s T_k - T_{k-1}
        % T'_{k+1} = 2 T_k + 2 s T'_k - T'_{k-1}
        Phi(1,:) = 1;
        if n >= 1
            Phi(2,:)  = s;
            dPhi(2,:) = 1;
        end
        for k = 1:n-1
            Phi(k+2,:)  = 2*s.*Phi(k+1,:) - Phi(k,:);
            dPhi(k+2,:) = 2*Phi(k+1,:) + 2*s.*dPhi(k+1,:) - dPhi(k,:);
        end

    case 'lagrange'
        % Supporting points of the reference element
        xn = nodal_distribution(n,ndist);
        xn = reshape(xn,[],1);
        for j = 1:n+1
            others = [1:j-1, j+1:n+1];
            % Factors (s - x_k)/(x_j - x_k), k ~= j
            F = (s - xn(others))./(xn(j) - xn(others));
            Phi(j,:) = prod(F,1);
            % Derivative: sum over m of the products that leave out the
            % m-th factor (prefix and suffix products, no division by
            % s - x_m, so that s may coincide with a node)
            pre = cumprod([ones(1,nQ); F(1:end-1,:)],1);
            suf = flipud(cumprod([ones(1,nQ); flipud(F(2:end,:))],1));
            dPhi(j,:) = sum((pre.*suf)./(xn(j) - xn(others)),1);
        end

    otherwise
        error('PolynomialBasis1D:basis', ...
              'Unknown basis "%s" (monomial, legendre, chebyshev, lagrange).', ...
              basis);
end

end
