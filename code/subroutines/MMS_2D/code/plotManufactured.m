function plotManufactured(u,f,domain,t,nP,varNames)
% plot solution u and right-hand side f
% u: function handle defining solution
% f: function handle defining rhs
% domain: limits of rectangular domain for plotting [[xmin,xmax];[ymin,ymax]] 
% t: time value to plot in case of transient problem
% nP: number of plotting points in each direction
% varNames: names of variables in u, just for axis labels

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

[xP,yP,uP,fP] = evaluateManufactured(uPlot,fPlot,domain,t,nP);
createPlot(xP,yP,uP,varNames);
createPlot(xP,yP,fP,{'f_x','f_y'});


    function [xP,yP,uP,fP] = evaluateManufactured(u,f,domain,t,nPlot)
        % evaluate function handles at discrete points

        xP = linspace(domain(1,1),domain(1,2),nPlot);                       % x-values for plotting
        yP = linspace(domain(2,1),domain(2,2),nPlot);                       % y-values for plotting
        [x,y] = meshgrid(xP,yP);                                            % array of coordinate values

        sizeU = numel(u(xP(1),yP(1),t));                                    % query output size of u and f
        sizeF = numel(f(xP(1),yP(1),t));

        uP = zeros([size(x),sizeU]);                              			% allocate numerical values of u, f
        fP = zeros([size(x),sizeF]);
        for i = 1:size(x,1)                                              	% loop coordinates
            for j = 1:size(x,2)
                uP(i,j,:) = u(x(i,j),y(i,j),t);                           	% evaluate u, f
                fP(i,j,:) = f(x(i,j),y(i,j),t);
            end
        end

    end


    function createPlot(x,y,z,label)
    % Define a new colormap
    nColors = 1024;   % 8, 16, 32, 64, 128, 256, 512, 104
%-----------------------------------Fast----------------------------------%
    % rgb1 = [0.05639932216773367, 0.056399092153948,   0.4700000908789252];
    % rgb2 = [0.5900001145322249,  0.07669636770019067, 0.11947505935767005];
%-----------------------------Smooth Cool Warm----------------------------%
    % rgb1 = [0.22999950386952345, 0.2989989340493756,   0.754000138575591];
    % rgb2 = [0.7060001359117047,  0.015991824033980695, 0.15000007192220008];
%-------------------------------------------------------------------------%
    % [map] = diverging_map(linspace(0,1,nColors),rgb1,rgb2);
    
    [map] = Fast(nColors);
    % [map] = SCW(nColors);

        % generate figures with some standard settings
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

end