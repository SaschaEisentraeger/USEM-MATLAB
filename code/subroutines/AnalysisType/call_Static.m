%% 2D SEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Sascha Duczek 11.11.2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t_start = cputime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Static analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Material data')
disp(LineDisp)

%------------------------Material properties: FEM-------------------------%
% Number of different material sets
num_mat = size(material_elset,1);

% Cell array for material data
C    = cell(num_mat,1);
rho  = cell(num_mat,1);
zeta = cell(num_mat,1);

for ijk = 1:num_mat
	% mechanical material properties
	if size(material_properties{ijk,3},2) == 2
        % isotropic material
        Mat_C = HookeMatrix3d(material_properties{ijk,3}(1),...
                                  material_properties{ijk,3}(2));
    else
        % transversal isotropic material
      	Mat_C = full(sparse([1 1 2 1 2 3 4 6 5],[1 2 2 3 3 3 4 6 5],...
                     material_properties{ijk,3}) + ... 
                     sparse([1 1 2 1 2 3 4 6 5],[1 2 2 3 3 3 4 6 5],...
                     material_properties{ijk,3})' - ...
           diag(diag(sparse([1 1 2 1 2 3 4 6 5],[1 2 2 3 3 3 4 6 5],...
                     material_properties{ijk,3}))));
	end

    % transform 3D constitutive matrices in 2D
    % plane stress condition
    if antype == 1
        [C2d] = PlaneStressConstitutiveMatrices(Mat_C);
    % plane strain condition
    elseif antype == 2
        [C2d] = PlaneStrainConstitutiveMatrices(Mat_C);
    end

    C(ijk,1)    = {C2d};
    rho(ijk,1)  = material_properties(ijk,1);
    zeta(ijk,1) = material_properties(ijk,5); 
end
clear Mat_C C2d

% Number of all finite elements
num_Elem = size(Topology_Node,1);

% Number of nodes
num_Nodes = max(max(Topology_Node));

% Number of degrees of freedom
num_DoF = num_Nodes*2;

%% Assembly of the system matrices
disp('Assembly of the system matrices')
disp(LineDisp)

[K,M,R,dte,t_mat] = AssemblyFEM(Topology_Node,p,xCoor,yCoor,num_Elem,...
                    num_DoF,C,rho,num_mat,mass_matrix,material_elset,...
                    alpha_Rayleigh,beta_Rayleigh,zeta,0);

%%
%% Homogeneous Dirichlet boundary conditions
%%
disp('Implement homogeneous Dirichlet boundary conditions')
disp(LineDisp)

[BC_DOFs,NOI] = DirichletBCs_homogeneous(boundary_conditions);

%%
%% Loads (Neumann boundary conditions)
%%
disp('Apply loads')
disp(LineDisp)

MAT.elset = material_elset;
MAT.prop  = material_properties;
[F,F_dt] = NeumannBCs(p,cload,dsload,body_force,gravity,MAT,num_DoF,...
           Topology_Node,xCoor,yCoor,num_Nodes);
       
%%
%% Non-homogeneous Dirichlet boundary conditions
%%
disp('Implement non-homogeneous Dirichlet boundary conditions')
disp(LineDisp)
[K,F] = DirichletBCs_nonhomogeneous(K,F,num_DoF,boundary_conditions,...
        xCoor,yCoor);
       
% Active degrees of freedom
ActiveDoF = setdiff((1:num_DoF)',unique(BC_DOFs));

%% ------------------------------------------------------------------------
%% Static Analysis

disp('Static analysis - Solving the linear system of equations')
disp(LineDisp)

u = zeros(size(K,1),1);
tic;
u_RB = K(ActiveDoF,ActiveDoF)\F(ActiveDoF);
t_cpu = toc;

% Include fixed degrees of freedom
idx = setdiff((1:num_DoF),ActiveDoF);
u(idx,1) = 0;

u(ActiveDoF,1) = u_RB;

ux = u(1:2:end);
uy = u(2:2:end);

%% Save the simulation results
if saveWS
disp('Save Workspace')
disp(LineDisp)
	[~,name,~] = fileparts(filename);

	method = 'Static';

	if ~exist(['.' filesep 'subroutines' filesep 'Results' filesep 'Static' ...
			   filesep name '_' method filesep],'dir')
		mkdir(['.' filesep 'subroutines' filesep 'Results' filesep 'Static' ...
			   filesep name '_' method filesep])
	end

	save(['.' filesep 'subroutines' filesep 'Results' filesep 'Static' ...
		  filesep name '_' method filesep name '_ndof' num2str(num_DoF)...
		  '.mat'],'-v7.3')
end

%% Visualization
disp('Visualization')
disp(LineDisp)

% Plot the deformed mesh
lx = max(xCoor) - min(xCoor);
ly = max(yCoor) - min(yCoor);
ux_max_abs = max(abs(ux));
uy_max_abs = max(abs(uy));
umag_max   = max(sqrt(ux.^2 + uy.^2));
ux_min_abs = min(abs(ux));
uy_min_abs = min(abs(uy));
umag_min   = min(sqrt(ux.^2 + uy.^2));
% Determine the scaling factor for the visualization of the deformation
scale = min([round(lx/ux_max_abs) round(ly/uy_max_abs)])/10;
if isinf(scale) || isnan(scale)
    scale = 1;
end
% scale = 0;

% Results that are plotted
if plot_disp
	strg1 = {'mag','ux','uy'};
else
    strg1 = {};
end
if plot_stress
    strg2 = {'Mises','sigma_x','sigma_y','tau_xy'};
else
    strg2 = {};
end
if plot_strain
    strg3 = {'epsilon_x','epsilon_y','gamma_xy'};
else
    strg3 = {};
end
strg = [strg1(:)',strg2(:)',strg3(:)'];

% Face transparency
n_alph = 1;
% Stress averaging (1: yes; 0: no)
flag_ave = 0;
[RES] = Contour_plot(p,xCoor,yCoor,u,scale,strg,Topology_Node,C,...
material_elset,flag_ave,n_alph);