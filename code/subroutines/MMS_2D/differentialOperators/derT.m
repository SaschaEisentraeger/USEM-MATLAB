function der = derT(u,n)
% Derivative with respect to t (nth-order)
syms t
der = diff(u,t,n);
end