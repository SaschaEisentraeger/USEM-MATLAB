function Y = RowNormalizedProduct(A, X)
%ROWNORMALIZEDPRODUCT Product A*X with the rows of A scaled to unit sum.
%
% ----------------- helpers -----------------

% Row-normalized sparse mat-vec: (D^{-1} A) * X with D(ii)=deg(i), adding self-loop
n = size(A,1);
A = A + speye(n); % include self for stability
Dinv = 1./max(sum(A,2),1);
Y = (spdiags(Dinv,0,n,n) * (A*X));
end
