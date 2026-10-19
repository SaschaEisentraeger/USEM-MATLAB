function der = derX(u,n)
% Derivative with respect to x (nth-order)
syms x
der = diff(u,x,n);
end