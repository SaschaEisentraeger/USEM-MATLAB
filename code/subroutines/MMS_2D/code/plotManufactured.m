function plotManufactured(u,f,domain,t,nP,varNames)
%PLOTMANUFACTURED Plot a manufactured solution and its right-hand side.
%
%   plotManufactured(u, f, domain, t, nP, varNames)
%   u:        function handle of the solution, u(x,y) or u(x,y,t)
%   f:        function handle of the right-hand side, same arguments
%   domain:   limits of the rectangular plot domain [xmin, xmax; ymin, ymax]
%   t:        time of the plot (transient problems)
%   nP:       number of plot points per direction
%   varNames: names of the components of u (axis labels), default {'u'}

if (nargin<6) || isempty(varNames)
    varNames = {'u'};
end

if nargin(u) == 2
    uPlot = @(x,y,t) u(x,y);
    fPlot = @(x,y,t) f(x,y);
else
    uPlot = u;
    fPlot = f;
end

[xP,yP,uP,fP] = ManufacturedGrid(uPlot,fPlot,domain,t,nP);
PlotManufacturedField(xP,yP,uP,varNames);
PlotManufacturedField(xP,yP,fP,{'f_x','f_y'});

end
