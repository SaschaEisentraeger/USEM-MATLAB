function [xP,yP,uP,fP] = ManufacturedGrid(u,f,domain,t,nPlot)
%MANUFACTUREDGRID Values of a manufactured solution on a plot grid.
%
%   [xP, yP, uP, fP] = ManufacturedGrid(u, f, domain, t, nPlot) evaluates
%   the function handles u(x,y,t) and f(x,y,t) on nPlot x nPlot points of
%   the rectangle domain = [xmin, xmax; ymin, ymax]. uP and fP have one
%   page (third index) per component.

xP = linspace(domain(1,1),domain(1,2),nPlot);       % x-values for plotting
yP = linspace(domain(2,1),domain(2,2),nPlot);       % y-values for plotting
[x,y] = meshgrid(xP,yP);                            % array of coordinates

sizeU = numel(u(xP(1),yP(1),t));                   % output size of u and f
sizeF = numel(f(xP(1),yP(1),t));

uP = zeros([size(x),sizeU]);
fP = zeros([size(x),sizeF]);
for i = 1:size(x,1)
    for j = 1:size(x,2)
        uP(i,j,:) = u(x(i,j),y(i,j),t);
        fP(i,j,:) = f(x(i,j),y(i,j),t);
    end
end

end
