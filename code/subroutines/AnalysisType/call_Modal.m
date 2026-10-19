%% 2D SEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Sascha Duczek 11.11.2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% t_start = cputime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Modal/Eigenvalue analysis
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
                    alpha_Rayleigh,beta_Rayleigh,zeta,1);

%%
%% Boundary conditions (Dirichlet boundary conditions)
%%
disp('Implement boundary conditions')
disp(LineDisp)

[BC_DOFs,NOI] = DirichletBCs_homogeneous(boundary_conditions);

%% ------------------------------------------------------------------------
%% Modal Analysis
disp('Modal analysis - Solving the eigenvalue problem')
disp(LineDisp)

% Active degrees of freedom
ActiveDoF = setdiff((1:num_DoF)',unique(BC_DOFs));

%% Solution of the eigenvalue problem
NoDoFActive = numel(ActiveDoF);
% The size of the Krylov subspace must be at least numEigen + 2
if NoDoFActive <= numEigen + 1
    numEigen = NoDoFActive - 2;
end
% Number of requested eigenvalues must be greater than 0
if numEigen == 0
	numEigen = NoDoFActive;
end
% Use eig to solve the eigenvalue problem
if 2*numEigen > NoDoFActive
    EVP_Solver.p = NoDoFActive;
    numEigen = NoDoFActive;
    flag_eig = true;
else
    flag_eig = false;
end

tic;
if nnz(R) == 0 || isempty(R)
%     % Sparse reverse Cuthill-McKee ordering --> smaller bandwidth
%     ro = symrcm(K);
%     K2 = K(ro,ro);
%     M2 = M(ro,ro);

%-------------------------------------------------------------------------%
    % Options for the eigenvalue solver:
    % Type of eigenvalues
    type = EVP_Solver.type;
	% Remove field from struct
	EVP_Solver = rmfield(EVP_Solver,'type');
	% Initial start vector
	if EVP_Solver.v0 == true
		% Assign random values
		EVP_Solver.v0 = randn(numel(ActiveDoF),1);
	else
		% Remove field from struct
		EVP_Solver = rmfield(EVP_Solver,'v0');
	end
	% Assign options
	opts = EVP_Solver;
%-------------------------------------------------------------------------%
    if flag_eig
        Kev = full(K(ActiveDoF,ActiveDoF));
        Mev = full(M(ActiveDoF,ActiveDoF));
        [MS,EV] = eig(Kev,Mev);
        EF = sqrt(abs(diag(EV))/(4*pi*pi));
    else
        [MS,EV] = eigs(K(ActiveDoF,ActiveDoF),M(ActiveDoF,ActiveDoF),...
		          numEigen,type,opts);
        EF = sqrt(abs(diag(EV))/(4*pi*pi));
    end

else

    % State-space method
    nDOFA = numel(ActiveDoF);
	A = [-speye(nDOFA)         , sparse(nDOFA,nDOFA); ...
          sparse(nDOFA,nDOFA)  , K(ActiveDoF,ActiveDoF)];
	B = [sparse(nDOFA,nDOFA)   , speye(nDOFA);        ...
         M(ActiveDoF,ActiveDoF), R(ActiveDoF,ActiveDoF)];
	[MS,EV] = eigs(A,-B,numEigen,'smallestabs');
    % lambda^2*M*u + lambda*R*u + K*u = 0 with the eigenvector [lambda*u; u]:
    % keep the displacement part and convert |lambda| (rad/s) to Hz
	EF = abs(diag(EV))/(2*pi);
    MS = MS(nDOFA+1:end,:);
end
% Sort eigenfrequencies and mode shapes
[EF,IdxSort] = sort(EF);
MS = MS(:,IdxSort);
t_cpu = toc;

% Mode shapes - add the fixed degrees of freedom
U = zeros(num_DoF,numEigen);
U(ActiveDoF,:) = MS;

%% Visualization
disp('Visualization')
disp(LineDisp)

% Mode shape
NumberEV = 1;
% (eigenvectors of the unsymmetric or damped problem can be complex)
u = real(U(:,NumberEV));
ux = u(1:2:end);
uy = u(2:2:end);

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

%% Save the simulation results
if saveWS
disp('Save Workspace')
disp(LineDisp)
	[~,name,~] = fileparts(filename);

	method = 'Modal';

	% Save the workspace and create a new directory if necessary
	if ~exist(['.' filesep 'subroutines' filesep 'Results' filesep 'Modal' ...
			   filesep name '_' method filesep],'dir')
		mkdir(['.' filesep 'subroutines' filesep 'Results' filesep 'Modal' ...
			   filesep name '_' method filesep])
	end

	if mass_matrix == 1
		% Consistent mass matrix (CMM)
		save(['.' filesep 'subroutines' filesep 'Results' filesep 'Modal' ...
			  filesep name '_' method filesep name '_ndof' num2str(num_DoF) ...
			  '_CMM.mat'],'-v7.3')
	elseif mass_matrix == 2
		% Lumped mass matrix (LMM) obtained by the Nodal Quadrature method
		save(['.' filesep 'subroutines' filesep 'Results' filesep 'Modal' ...
			  filesep name '_' method filesep name '_ndof' num2str(num_DoF) ...
			  '_LMM_NQ.mat'],'-v7.3')
	elseif mass_matrix == 3
		% Lumped mass matrix (LMM) obtained by the Row Sum method
		save(['.' filesep 'subroutines' filesep 'Results' filesep 'Modal' ...
			  filesep name '_' method filesep name '_ndof' num2str(num_DoF) ...
			  '_LMM_RowSum.mat'],'-v7.3')
	elseif mass_matrix == 4
		% Lumped mass matrix (LMM) obtained by the HRZ method
		save(['.' filesep 'subroutines' filesep 'Results' filesep 'Modal' ...
			  filesep name '_' method filesep name '_ndof' num2str(num_DoF) ...
			  '_LMM_HRZ.mat'],'-v7.3')
	end
end