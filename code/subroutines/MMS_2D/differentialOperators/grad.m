function gradU = grad(u)
% Gradient operator (2D)
syms x y

nu = numel(u);

if nu == 1
    gradU = [diff(u,x); diff(u,y)];
elseif nu == 2
    gradU = [diff(u(1),x), diff(u(2),x); diff(u(1),y), diff(u(2),y)];
end

end