%% Finite element mesh of 2D spectral elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% p:                    Polynomial degree structure
% p.pxi:                Polynomial degree in xi
% p.peta:               Polynomial degree in eta
% p.ndist:              Nodal distribution
% Topology_Node:        Nodal topology (connectivity/coincidence) matrix
% xCoor:                x-coordinates of the nodes
% yCoor:                y-coordinates of the nodes
% cload:                Concentrated loads (point forces)
% dsload:               Distributed loads (pressure, traction)
% body_force:           Volume loads (arbitrary direction and elements)
% gravity:              Volume load (all elements)
% boundary_conditions:	Data structure containing the boundary condition
%                       data (read from Abaqus *.inp file)
% flagBCplot:           1: plot boundary conditions;
%                       0: do not plot boundary conditions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FE_Mesh(p,xCoor,yCoor,Topology_Node,cload,dsload,body_force,...
         gravity,boundary_conditions,flagBCplot)

fh = figure('Name','Finite Element Mesh');
fh.WindowState = 'maximized';
hold on
grid on
box on
xlabel('x [m]','FontName','Times','FontSize',40)
ylabel('y [m]','FontName','Times','FontSize',40)
set(gca,'FontName','Times','FontSize',24);

% Polynomial degrees in xi and eta
pxi = p.pxi;
peta = p.peta;

% Number of finite elements
NoElem = size(Topology_Node,1);

% Sort the nodes for the visualization
% Nodes on edge 1
NE1 = [1, 5:4+(pxi-1), 2];
% Nodes on edge 2
NE2 = [2, 5+(pxi-1):4+(pxi-1)+(peta-1), 3];
% Nodes on edge 3
NE3 = [4, 5+(pxi-1)+(peta-1):4+2*(pxi-1)+(peta-1), 3];
% Nodes on edge 4
NE4 = [1, 5+2*(pxi-1)+(peta-1):4+2*(pxi-1)+2*(peta-1), 4];
% Order of the nodes for the plot (counter-clockwise)
N = [NE1, NE2(2:end), flip(NE3(2:end)), flip(NE4(2:end))];

% Number of nodes on the boundary of the element
NB = 4 + 2*(pxi-1) + 2*(peta-1);

% Face color
FC = 'w';

% Collect plot data for all elements
vert = cell(NoElem,1);
fac = cell(NoElem,1);
for i = 1:NoElem
    % Nodal coordinates of all nodes of spectral element i
	x = xCoor(Topology_Node(i,:));
    y = yCoor(Topology_Node(i,:));
    
    % Save data for each element
    vert{i} = [x(N), y(N)];
    fac{i} = (1:NB+1);
end

% Plot all elements
% Convert cell array to mat
vert = cell2mat(vert);
fac = cell2mat(fac);
add = repmat([0;cumsum(max(fac,[],2))],1,NB+1);
fac = fac + add(1:NoElem,:);
faccolor = [100/255 100/255 100/255];
patch('Faces',fac,'Vertices',vert,'FaceVertexCData',faccolor,...
      'FaceColor',FC,'EdgeColor','k','LineWidth',1)

% Plot marker for the Dirichlet and Neumann boundary conditions
if flagBCplot == 1
    % Number of concentrated loads
    ncloads = numel(cload);
    for i = 1:ncloads
        x = xCoor(cload(i).nset);
        y = yCoor(cload(i).nset);
        plot(x,y,'rx','MarkerSize',5)
    end

    % Number of distributed loads
    ndsloads = numel(dsload);
    for i = 1:ndsloads
        x = xCoor(dsload(i).nset);
        y = yCoor(dsload(i).nset);
        plot(x,y,'bs','MarkerSize',5)
    end

    % Number of volume/body loads
    nbody = numel(body_force);
    for i = 1:nbody
        nset = unique(Topology_Node(body_force(i).elset,:));
        x = xCoor(nset);
        y = yCoor(nset);
        plot(x,y,'kd','MarkerSize',5)
    end

    % Number of gravity loads
    ngravity = numel(gravity);
    for i = 1:ngravity
        nset = unique(Topology_Node(gravity(i).elset,:));
        x = xCoor(nset);
        y = yCoor(nset);
        plot(x,y,'mv','MarkerSize',5)
    end

    % Number of Dirichlet boundary conditions
    nDBC = numel(boundary_conditions);
    for i = 1:nDBC
        x = xCoor(boundary_conditions(i).nset);
        y = yCoor(boundary_conditions(i).nset);
        plot(x,y,'go','MarkerSize',5)
    end
end

a = min(xCoor) - abs(max(xCoor)-min(xCoor))*0.05;
b = max(xCoor) + abs(max(xCoor)-min(xCoor))*0.05;
c = min(yCoor) - abs(max(xCoor)-min(xCoor))*0.05;
d = max(yCoor) + abs(max(xCoor)-min(xCoor))*0.05;

axis equal 
hold off

xlim([a b])
ylim([c d])

% SAVE figure
% axis off
% h = gcf;
% iptsetpref('ImshowBorder','tight')
% print(h,'-depsc','-painters','Mesh.eps')
% 
% savefig('Mesh.fig')

end