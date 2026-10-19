function der = derY(u,n)
% Derivative with respect to y (nth-order)
syms y
der = diff(u,y,n);
end