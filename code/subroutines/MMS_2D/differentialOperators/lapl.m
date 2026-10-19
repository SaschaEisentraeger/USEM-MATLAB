function laplU = lapl(u)
% Laplace operator (2D)
syms x y
if size(u,1) < size(u,2)
    u = u.';
end
laplU = diff(u,x,2)+diff(u,y,2);
end