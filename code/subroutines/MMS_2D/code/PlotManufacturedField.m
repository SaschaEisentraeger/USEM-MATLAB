function PlotManufacturedField(x,y,z,label)
%PLOTMANUFACTUREDFIELD Surface plots of the components of a field.
%
%   PlotManufacturedField(x, y, z, label) opens one figure per page of z
%   (component) with a surface plot over the grid x, y. label holds the
%   component names; with fewer labels than components, the components
%   are numbered (label{1}_1, label{1}_2, ...).

% Diverging colormap with the end colours of the "Fast" map (Moreland);
% "Smooth Cool Warm": rgb1 = [0.230 0.299 0.754], rgb2 = [0.706 0.016 0.150]
nColors = 1024;
rgb1 = [0.05639932216773367, 0.056399092153948,   0.4700000908789252];
rgb2 = [0.5900001145322249,  0.07669636770019067, 0.11947505935767005];
map = diverging_map(linspace(0,1,nColors),rgb1,rgb2);

sizeZ = size(z,3);
sizeLabel = numel(label);
for i = 1:sizeZ
    figure('defaulttextinterpreter','latex','windowstate','maximized');
    surf(x,y,z(:,:,i),'FaceColor','interp');
    colormap(map)
    xlabel('$x$')
    ylabel('$y$')
    if sizeZ > sizeLabel
        zlabel(['$',label{1},'_',num2str(i),'$'])
    else
        zlabel(['$',label{i},'$'])
    end
    set(gca,'FontSize',24)
    set(gca,'TickLabelInterpreter','latex');
    axis square
    box on
    rotate3d on
end

end
