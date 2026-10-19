%% 2D SEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Sascha Duczek 11.11.2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t_start = cputime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Dynamic (transient) analysis - Time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                    alpha_Rayleigh,beta_Rayleigh,zeta,2);

%%
%% Boundary conditions (Dirichlet boundary conditions)
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
[F,FNeumann_dt] = NeumannBCs(p,cload,dsload,body_force,gravity,MAT,...
                  num_DoF,Topology_Node,xCoor,yCoor,num_Nodes,flag_mRHS);

%%
%% Boundary conditions (Dirichlet boundary conditions)
%%
disp('Implement nonhomogeneous Dirichlet boundary conditions')
disp(LineDisp)
[DOF_C,DOF_C_Dyn,FDBCdyn,FDBCdyn_dt] = ...
DirichletBCs_nonhomogeneous_dynamic(K,R,M,Topology_Node,...
boundary_conditions,xCoor,yCoor);

% Force vector
F = F + FDBCdyn;
F_dt.Neumann = FNeumann_dt;
F_dt.K       = FDBCdyn_dt.K;
F_dt.C       = FDBCdyn_dt.C;
F_dt.M       = FDBCdyn_dt.M;

%% Save all data from workspace for further post-processing
disp('Save Workspace')
disp(LineDisp)

[~,name,~] = fileparts(filename);
save(['.' filesep 'subroutines' filesep 'Temp' filesep name ...
      '_transient_analysis_ndof' num2str(num_DoF) '.mat'], '-v7.3')

save(['.' filesep 'subroutines' filesep 'Temp' filesep 'K.mat'],'K','-v7.3')
save(['.' filesep 'subroutines' filesep 'Temp' filesep 'M.mat'],'M','-v7.3')
save(['.' filesep 'subroutines' filesep 'Temp' filesep 'R.mat'],'R','-v7.3')
save(['.' filesep 'subroutines' filesep 'Temp' filesep 'F.mat'],'F','-v7.3')
save(['.' filesep 'subroutines' filesep 'Temp' filesep 'F_dt.mat'],'F_dt','-v7.3')

%% Solution of the semi-discrete equation of motion
% Active degrees of freedom
ActiveDoF = setdiff((1:num_DoF)',unique([BC_DOFs;DOF_C.BC;DOF_C_Dyn.BC]));

% Initial conditions
U0 = zeros(num_DoF,1);
V0 = zeros(num_DoF,1);

% Nodes of Interest
if isempty(NOI) == 1 && save_var == 0
    NOI = unique(round(ActiveDoF/2));
    % Results of the first three active nodes are saved if nothing is
    % prescribed
    NOI = NOI(1:min(3,numel(NOI)));
end

%% Save variables for the solver
disp('Save input data for the solver')
disp(LineDisp)

save(['.' filesep 'subroutines' filesep 'Temp' filesep ...
      'Solver_input.mat'],'t_end','Topology_Node','NOI','ActiveDoF',...
	  'num_DoF','xCoor','yCoor','save_var','SignalType','filename',...
	  'mass_matrix','dtf','dte','TimeIntegration','ExciteParam','U0',...
      'V0','flag_mRHS','-v7.3')
  
%% Call solver for the three-dimensional equations of motion 
%% Time-integration
[t_comp,y_comp,t_cpu] = Solver_Transient;

if ~isempty(DOF_C_Dyn.BC) && save_var ~= 0
    % Add the prescribed degrees of freedom (time-dependent)
    Upre = DOF_C_Dyn.Amp*ExcitationSignal(SignalType,ExciteParam,t_comp);
    y_comp(:,DOF_C_Dyn.BC) = Upre';
end
if ~isempty(DOF_C.BC) && save_var ~= 0
    % Add the prescribed degrees of freedom (time-independent)
    Upre = DOF_C.Amp;
    Upre = repmat(Upre,1,numel(t_comp));
    y_comp(:,DOF_C.BC) = Upre';
end

% Integration method number as a string
IntMeth = num2str(TimeIntegration.Method);

%% Save the simulation results
if saveWS
	if ~exist(['.' filesep 'subroutines' filesep 'Results' filesep ...
			   'Dynamic' filesep name '_' IntMeth filesep],'dir')
		mkdir(['.' filesep 'subroutines' filesep 'Results' filesep ...
			   'Dynamic' filesep name '_' IntMeth filesep])
	end

	if mass_matrix == 1
		% Consistent mass matrix (CMM)
		save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
			  filesep name '_' IntMeth filesep name '_ndof' num2str(num_DoF)...
			  '_CMM.mat'],'-v7.3')
	elseif mass_matrix == 2
		% Lumped mass matrix (LMM) obtained by the Nodal Quadrature method
		save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
			  filesep name '_' IntMeth filesep name '_ndof' num2str(num_DoF)...
			  '_LMM_NQ.mat'],'-v7.3')
	elseif mass_matrix == 3
		% Lumped mass matrix (LMM) obtained by the Row Sum method
		save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
			  filesep name '_' IntMeth filesep name '_ndof' num2str(num_DoF)...
			  '_LMM_RowSum.mat'],'-v7.3')
	elseif mass_matrix == 4
		% Lumped mass matrix (LMM) obtained by the HRZ method
		save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
			  filesep name '_' IntMeth filesep name '_ndof' num2str(num_DoF)...
			  '_LMM_HRZ.mat'],'-v7.3')
	end
end

if ~TimeIntegration.SolverSave
    % Delete temporarily saved files
    delete(['.' filesep 'subroutines' filesep 'Temp' filesep '*.mat'])
end

%% Export of the y-displacements of the nodes of interest (sorted by x)
% Only possible if the results are saved and both displacement components
% of every node of interest have been stored
if saveWS && save_var == 0 && size(y_comp,2) == 2*numel(NOI)
    % Sort according to x-coordinate
    [~,idx] = sort(xCoor(NOI));
    IDsort = zeros(2*numel(idx),1);
    IDsort(1:2:end) = 2*idx-1;
    IDsort(2:2:end) = 2*idx;
    % Save only y-components of the displacement
    U = y_comp(:,IDsort(2:2:end),:);
    T = t_comp;
    save(['.' filesep 'subroutines' filesep 'Results' filesep 'Dynamic' ...
          filesep name '_' IntMeth filesep name '_ndof' num2str(num_DoF)...
          '_Uy.mat'],'T','U','-v7.3')
end