%% Visualization of the spectral element results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 27.11.2014 (University of Magdeburg)
%
% Contour plot of the primary and secondary results.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% p:                Polynomial degree template
% p.pxi:            Polynomial degree in xi
% p.peta:           Polynomial degree in eta
% p.ndist:          Nodal distribution
% xCoor:            x-coordinates of all nodes
% yCoor:          	y-coordinates of all nodes
% U:                Displacement vector
% scale:            Scaling factor of the displacements
% strg_cell:    	Result variable that is plotted: ux, uy, mag, sigma_x,
%                   sigma_y, tau_xy, Mises, epsilon_x, epsilon_y, gamma_xy
% Topology_Node:   	Nodal topology (connectivity/coincidence) matrix
% C:                Elasticity matrix
% material_elset:   Element sets for the different materials
% flag_ave:        	Stress averaging (1) yes; (0) no
% n_alph:           Alpha value for the faces (transparency)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [res] = Contour_plot(p,xCoor,yCoor,U,scale,strg_cell,...
                 Topology_Node,C,material_elset,flag_ave,n_alph)

% Plot the deformed triangular mesh
% 0: Do not plot original finite element mesh but instead plot the fine
% visualization mesh (surface triangulation)
% 1: Plot the original finite element mesh
plot_mesh = 1;

% For finer visualization mesh: nFac > 1 (integer)
nFac = 20;

% Number of results that are visualized
nVis = numel(strg_cell);

% Define a new colormap
nColors = 32;   % 128
rgb1 = [0.2298057,   0.298717966, 0.753683153];
rgb2 = [0.705673158, 0.01555616,  0.150232812];
[map] = diverging_map(linspace(0,1,nColors),rgb1,rgb2);

% Determine if stresses and strains are requested as output
c_str = char(strg_cell);
c_str = reshape(c_str',1,[]);
if contains(c_str,'sigma') || contains(c_str,'tau') || ...
   contains(c_str,'Mises') || contains(c_str,'epsilon') || ...
   contains(c_str,'gamma')
    flag_Stress = 1;
else
    flag_Stress = 0;
end

% Polynomial degrees in xi, eta
pxi   = p.pxi;
peta  = p.peta;

% Degrees of freedom topology data
Topo = Topology_DoF_2D(Topology_Node);

if plot_mesh == 1
    % Generate edge topology data
    [Topology_Edge,NumEdges] = Topology_Edge_2D(Topology_Node);
end

% Sub-division of each element for a more precise visualization of high
% order results
ndiv(1) = nFac*pxi+1;  	% Sub-division of the element in xi-direction
ndiv(2) = nFac*peta+1; 	% Sub-division of the element in eta-direction

% Number of elements
NumElem = size(Topo,1);

% Coordinate matrices
Xplot = cell(NumElem,1);
Yplot = cell(NumElem,1);
% Displacement matrices
Displacement_X = cell(NumElem,1);
Displacement_Y = cell(NumElem,1);
% Stress matrices
Stress_X  = cell(NumElem,1);
Stress_Y  = cell(NumElem,1);
Stress_XY = cell(NumElem,1);
% Strain matrices
Strain_X  = cell(NumElem,1);
Strain_Y  = cell(NumElem,1);
Strain_XY = cell(NumElem,1);
% Store data of the triangulation
TRI = cell(NumElem,1);

if plot_mesh == 1
    % Plot mesh (only element edges)
    Mesh1 = cell(NumElem,1);
    Mesh2 = cell(NumElem,1);
    Mesh3 = cell(NumElem,1);
    Mesh4 = cell(NumElem,1);
end

% Elasticity matrix for each element
EMAT = cell(NumElem,1);

% Nodal coordinates for each element
xElem = cell(NumElem,1);
yElem = cell(NumElem,1);

% Nodal displacement values for each element
u = cell(NumElem,1);

% Number of material element sets
num_mat = size(material_elset,1);

% Number of degrees of freedom per finite element
ndof = zeros(NumElem,1);

% Point counter
NP = zeros(NumElem,1);

for i = 1:NumElem
    % Material data - Elasticity matrix
    for klm = 1:num_mat
        % Element sets are expanded by the reader (incl. generate)
        val_Elset = material_elset{klm};
        val =  (val_Elset == i);
        if sum(val) == 1
            break;
        end
    end
    EMAT{i} = cell2mat(C(klm));

    % Nodal coordinates of element i
    xElem{i} = xCoor(Topology_Node(i,:));
    yElem{i} = yCoor(Topology_Node(i,:));
    
    % Nodal displacement values for element i
    u{i} = U(Topo(i,:));
    
    % Number of degrees of freedom for element i
    ndof(i) = numel(u{i});
    
    % Point counter
    NP(i) = (ndiv(1)*ndiv(2))*(i-1);
end

% Compute displacement and stress/strain results in parallel
parfor i = 1:NumElem
% for i = 1:NumElem
    % Compute displacements
    [x,y,ux,uy,local_coor] = Calc_results(p,xElem{i},yElem{i},u{i},ndiv);

    % Save results
    Xplot{i} = x;
    Yplot{i} = y;
    Displacement_X{i} = ux;
    Displacement_Y{i} = uy;
    
    % Triangulate the surfaces of the element (triangulation based on the
    % reference element gives better results)
    DT = delaunayTriangulation(local_coor(:,1),local_coor(:,2));
    tri = DT.ConnectivityList;
    % Save triangulation
    TRI{i} = tri  + NP(i);
    
    if plot_mesh == 1
        % Element boundary
        ElemHull = boundary(local_coor(:,1),local_coor(:,2));
        % Number of nodes on the boundary (first and last node are the same
        % in the element hull --> closed polygonal line)
        nNB = numel(ElemHull) - 1;
        % Local coordinates of the boundary nodes
        CoorElemHull = local_coor(ElemHull(1:nNB),:);
        % Element edges
        % Edge 1
        idx = find(CoorElemHull(:,2) == -1);
        xl = local_coor(ElemHull(idx),1);
        [~,idx_sort] = sort(xl);
        Mesh1{i} = [x(ElemHull(idx(idx_sort)))+ux(ElemHull(idx(idx_sort)))*scale,...
                    y(ElemHull(idx(idx_sort)))+uy(ElemHull(idx(idx_sort)))*scale];
        % Edge 2
        idx = find(CoorElemHull(:,1) == 1);
        xl = local_coor(ElemHull(idx),2);
        [~,idx_sort] = sort(xl);
        Mesh2{i} = [x(ElemHull(idx(idx_sort)))+ux(ElemHull(idx(idx_sort)))*scale,...
                    y(ElemHull(idx(idx_sort)))+uy(ElemHull(idx(idx_sort)))*scale];
        % Edge 3
        idx = find(CoorElemHull(:,2) == 1);
        xl = local_coor(ElemHull(idx),1);
        [~,idx_sort] = sort(xl);
        Mesh3{i} = [x(ElemHull(idx(idx_sort)))+ux(ElemHull(idx(idx_sort)))*scale,...
                    y(ElemHull(idx(idx_sort)))+uy(ElemHull(idx(idx_sort)))*scale];
        % Edge 4
        idx = find(CoorElemHull(:,1) == -1);
        xl = local_coor(ElemHull(idx),2);
        [~,idx_sort] = sort(xl);
        Mesh4{i} = [x(ElemHull(idx(idx_sort)))+ux(ElemHull(idx(idx_sort)))*scale,...
                    y(ElemHull(idx(idx_sort)))+uy(ElemHull(idx(idx_sort)))*scale];
    end
    
    if flag_Stress == 1
        % Calculate stress and strain results
        [stress,strain] = Results_Stress_Strain(p,xElem{i},yElem{i},...
                          u{i},EMAT{i},ndof(i),local_coor);
         % Save results
         % Stress matrices
        Stress_X{i}  = stress(1,:)';
        Stress_Y{i}  = stress(2,:)';
        Stress_XY{i} = stress(3,:)';
        % Strain matrices
        Strain_X{i}  = strain(1,:)';
        Strain_Y{i}  = strain(2,:)';
        Strain_XY{i} = strain(3,:)';
    end
end

% Save output
% Coordinates
res.X = Xplot;
res.Y = Yplot;
% Displacements
res.UX = Displacement_X;
res.UY = Displacement_Y;
% Stress values
res.SXX = Stress_X;
res.SYY = Stress_Y;
res.SXY = Stress_XY;
% Strain values
res.EXX = Strain_X;
res.EYY = Strain_Y;
res.EXY = Strain_XY;

% Save all data in vectors
Xall = cell2mat(Xplot);
Yall = cell2mat(Yplot);
UXall = cell2mat(Displacement_X);
UYall = cell2mat(Displacement_Y);
SXall = cell2mat(Stress_X);
SYall = cell2mat(Stress_Y);
SXYall = cell2mat(Stress_XY);
EXall = cell2mat(Strain_X);
EYall = cell2mat(Strain_Y);
EXYall = cell2mat(Strain_XY);
TRIall = cell2mat(TRI);

% stress and strain averaging
if flag_ave == 1 && flag_Stress == 1
	[SXall,SYall,SXYall,EXall,EYall,EXYall] = stress_averaging(Xall,...
	 Yall,SXall,SYall,SXYall,EXall,EYall,EXYall);
end

% Auxiliary variables
Zall  = zeros(size(Xall));
UZall = zeros(size(UXall));

for i = 1:nVis
    % Result variable
    strg = strg_cell{i};
    
    % Settings for the figure   
    fh = figure('Name',['Finite Element Mesh: ' strg]);
	fh.WindowState = 'maximized';
    hold on
    grid on
    box on
    
    if plot_mesh == 1
        strEdge = 'none';
    else
        strEdge = 'k';
    end
	
    % Plot results
    if strcmp(strg,'ux') == 1
        trisurf(TRIall,Xall+UXall*scale,Yall+UYall*scale,Zall+UZall*scale,UXall,...
        'FaceColor','interp','EdgeColor',strEdge,'FaceAlpha',n_alph);
		
		disp(['u_{xmax} = ' num2str(max(UXall))])
        disp(['u_{xmin} = ' num2str(min(UXall))])

    elseif strcmp(strg,'uy') == 1
        trisurf(TRIall,Xall+UXall*scale,Yall+UYall*scale,Zall+UZall*scale,UYall,...
        'FaceColor','interp','EdgeColor',strEdge,'FaceAlpha',n_alph);
		
		disp(['u_{ymax} = ' num2str(max(UYall))])
        disp(['u_{ymin} = ' num2str(min(UYall))])

    elseif strcmp(strg,'mag') == 1
        MAG = sqrt(UXall.^2 + UYall.^2 + UZall.^2);
        trisurf(TRIall,Xall+UXall*scale,Yall+UYall*scale,Zall+UZall*scale,MAG,...
        'FaceColor','interp','EdgeColor',strEdge,'FaceAlpha',n_alph);
		
		disp(['u_{mag_max} = ' num2str(max(MAG))])
        disp(['u_{mag_min} = ' num2str(min(MAG))])

    elseif strcmp(strg,'Mises') == 1
		Mises = sqrt(SXall.^2 + SYall.^2 - SXall.*SYall + 3*SXYall.^2);
        trisurf(TRIall,Xall+UXall*scale,Yall+UYall*scale,Zall+UZall*scale,Mises,...
        'FaceColor','interp','EdgeColor',strEdge,'FaceAlpha',n_alph);
		
		disp(['Mises_{max} = ' num2str(max(Mises))])
        disp(['Mises_{min} = ' num2str(min(Mises))])

    elseif strcmp(strg,'sigma_x') == 1
        trisurf(TRIall,Xall+UXall*scale,Yall+UYall*scale,Zall+UZall*scale,SXall,...
        'FaceColor','interp','EdgeColor',strEdge,'FaceAlpha',n_alph);
		
		disp(['\sigma_x_{max} = ' num2str(max(SXall))])
        disp(['\sigma_x_{min} = ' num2str(min(SXall))])

    elseif strcmp(strg,'sigma_y') == 1
        trisurf(TRIall,Xall+UXall*scale,Yall+UYall*scale,Zall+UZall*scale,SYall,...
        'FaceColor','interp','EdgeColor',strEdge,'FaceAlpha',n_alph);
		
		disp(['\sigma_y_{max} = ' num2str(max(SYall))])
        disp(['\sigma_y_{min} = ' num2str(min(SYall))])

    elseif strcmp(strg,'tau_xy') == 1
        trisurf(TRIall,Xall+UXall*scale,Yall+UYall*scale,Zall+UZall*scale,SXYall,...
        'FaceColor','interp','EdgeColor',strEdge,'FaceAlpha',n_alph);
		
		disp(['\tau_xy_{max} = ' num2str(max(SXYall))])
        disp(['\tau_xy_{min} = ' num2str(min(SXYall))])

    elseif strcmp(strg,'epsilon_x') == 1
        trisurf(TRIall,Xall+UXall*scale,Yall+UYall*scale,Zall+UZall*scale,EXall,...
        'FaceColor','interp','EdgeColor',strEdge,'FaceAlpha',n_alph);
		
		disp(['\epsilon_x_{max} = ' num2str(max(EXall))])
        disp(['\epsilon_x_{min} = ' num2str(min(EXall))])

    elseif strcmp(strg,'epsilon_y') == 1
        trisurf(TRIall,Xall+UXall*scale,Yall+UYall*scale,Zall+UZall*scale,EYall,...
        'FaceColor','interp','EdgeColor',strEdge,'FaceAlpha',n_alph);
		
		disp(['\epsilon_y_{max} = ' num2str(max(EYall))])
        disp(['\epsilon_y_{min} = ' num2str(min(EYall))])

    elseif strcmp(strg,'gamma_xy') == 1
        trisurf(TRIall,Xall+UXall*scale,Yall+UYall*scale,Zall+UZall*scale,EXYall,...
        'FaceColor','interp','EdgeColor',strEdge,'FaceAlpha',n_alph);
		
		disp(['\gamma_xy_{max} = ' num2str(max(EXYall))])
        disp(['\gamma_xy_{min} = ' num2str(min(EXYall))])

    end
    
    % Plot element edges
    if plot_mesh == 1
        % Avoid plotting edges multiple times - list of plotted edges
        idx_Edge = zeros(NumEdges,1);
        for j = 1:NumElem
            % Edge numbers
            Edge = Topology_Edge(j,:);
            % Check if these edges have been plotted before
            idx_plot = ismember(Edge,idx_Edge);
            % Edge 1
            if idx_plot(1) == 0
                plot(Mesh1{j}(:,1),Mesh1{j}(:,2),'k','LineWidth',1)
            end
            % Edge 2
            if idx_plot(2) == 0
                plot(Mesh2{j}(:,1),Mesh2{j}(:,2),'k','LineWidth',1)
            end
            % Edge 3
            if idx_plot(3) == 0
                plot(Mesh3{j}(:,1),Mesh3{j}(:,2),'k','LineWidth',1)
            end
            % Edge 4
            if idx_plot(4) == 0
                plot(Mesh4{j}(:,1),Mesh4{j}(:,2),'k','LineWidth',1)
            end
            % Update the edge list
            idx_Edge(Edge) = Edge;
        end
    end
    
    % Figure settings
    axis equal

    % Limits
    x_max = max(Xall+UXall*scale);
    x_min = min(Xall+UXall*scale);
    y_max = max(Yall+UYall*scale);
    y_min = min(Yall+UYall*scale);
    length_x = abs(x_max - x_min);
    xlim([x_min-0.05*length_x x_max+0.05*length_x])
    ylim([y_min-0.05*length_x y_max+0.05*length_x])

    view(2)

    % Axis labels
    xlabel('x','FontName','Times','FontSize',24,'FontWeight','bold','Color','k')
    ylabel('y','FontName','Times','FontSize',24,'FontWeight','bold','Color','k')

    % Get handle to current figure
    hFig =  gcf;
    % Get handle to current axis
    hAxes = gca;
    % Add the axes controls back to the figure toolbar (Matlab 2018b)
    addToolbarExplorationButtons(hFig)
    % Hide the integrated axes toolbar
    hAxes.Toolbar.Visible = 'off';
    set(hAxes,'FontName','Times','FontSize',24);
    set(hFig,'PaperPositionMode','auto')
    % Color bar with rainbow map
    colorbar
    colormap(map)
end

end