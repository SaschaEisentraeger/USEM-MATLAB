%% Assembly of the system matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sascha Duczek 11.11.2014 (University of Magdeburg)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input
%
% Topology_Node:        Nodal topology (connectivity/coincidence) matrix
% p:                 	Polynomial degree structure
% p.pxi:                Polynomial degree in xi
% p.peta:           	Polynomial degree in eta
% p.ndist:              Nodal distribution
% xCoor:                x-coordinates of the nodes
% yCoor:                y-coordinates of the nodes
% num_Elem:             Number of elements
% num_DoF:              Number of degrees of freedom
% C:                    Elasticity matrix
% rho:                  Mass density
% num_mat:              Number of different materials
% mass_matrix:          Mass matrix flags 
%                       1: Consistent mass matrix
%                       2: Lumped mass matrix (nodal quadrature)
%                       3: Lumped mass matrix (row-sum method)
%                       4: Lumped mass matrix (HRZ-method/diagonal scaling)
% material_sets:        Element sets for different materials
% alpha_Rayleigh:       Mass specific damping parameter
% beta_Rayleigh:        Stiffness specific damping parameter
% zeta:					Damping ratio --> indicates if damping is present
% 						in the finite element (value is not considered -->
% 						alpha and beta are used to set up Rayleigh damping)
% flag_comp:            Computation flag
%                       0: Stiffness matrix
%                       1: Stiffness + mass matrix
%                       2: Stiffness + mass matrix; estimation of dt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Output
%
% K:        Global stiffness matrix
% M:        Global mass matrix
% R:        Global damping matrix
% dt:       Time increment computed from the highest numerical (element)
%           eigenfrequency
% t_cpu:   	Computational time for the assembly process
% ElemData: Element data (only computed if requested; needed by the
%           adaptive hybrid explicit method):
%           ElemData.Ke, ElemData.Me: element stiffness and mass matrices
%           (cell arrays), ElemData.Topo: DOF connectivity,
%           ElemData.damped: elements with Rayleigh damping,
%           ElemData.alpha, ElemData.beta: Rayleigh parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,M,R,dt,t_cpu,ElemData] = AssemblyFEM(Topology_Node,p,xCoor,yCoor,...
          num_Elem,num_DoF,C,rho,num_mat,mass_matrix,material_sets,...
          alpha_Rayleigh,beta_Rayleigh,zeta,flag_comp)

% Start the time measurement				   
t0 = tic;

%% Connectivity matrix for the degrees of freedom
[Topo] = Topology_DoF_2D(Topology_Node);
% Ensure that the connectivity matrix contains only integer
Topo = uint32(Topo);

% Number of degrees of freedom per finite element
ndof = numel(Topo(1,:));

% Data structure for the assembly
% Topology data
ROW = cell(num_Elem,1);
COL = cell(num_Elem,1);
% Stiffness matrix
Ke = cell(num_Elem,1);
% Mass matrix
Me = cell(num_Elem,1);
% Damping matrix
Re = cell(num_Elem,1);

% Elasticity matrix of each element is saved in a cell array
C_cell = cell(num_Elem,1);

% Mass density of each element is saved in a vector
rho_vec = zeros(num_Elem,1);

% Damping ratio of each element is saved in a vector
zeta_vec = zeros(num_Elem,1);

% Nodal coordinates of each element are saved in a cell array
x_cell = cell(num_Elem,1);
y_cell = cell(num_Elem,1);

% Time increment is saved in a vector
dt_vec = zeros(num_Elem,1);

% Options for the quadratic eigenvalue problem
opts.tol = 1e-6;

%%
disp('Prepare element data')
% Element data
for i = 1:num_Elem
     % Material properties of the finite element
    for klm = 1:num_mat
        % Element sets are expanded by the reader (incl. generate)
        val_Elset = material_sets{klm};
        val = (val_Elset == i);
        if sum(val) == 1
            break;
        end
    end
    
    % Assign material properties
    Hooke = cell2mat(C(klm));
    dens  = cell2mat(rho(klm)); if isempty(dens); dens = 1; end
	z     = cell2mat(zeta(klm)); if isempty(z); z = 0; end
    
	% Elasticity matrix
    C_cell{i} = Hooke;
	
    % Mass density
	rho_vec(i) = dens;
    
    % Damping ratio
	zeta_vec(i) = z;
    
    % Nodal coordinates
    x_cell{i} = xCoor(Topology_Node(i,:));
	y_cell{i} = yCoor(Topology_Node(i,:));
end

%%
disp('Compute element matrices in parallel')
% Compute element matrices in parallel (serially without a parallel pool;
% replace parfor by for to debug inside the loop)
parfor i = 1:num_Elem
	% Clear element matrices
    K_element = [];
    M_element = [];
    R_element = [];
	
	% Element data
    x = x_cell{i};
    y = y_cell{i};
    Hooke = C_cell{i};
    dens  = rho_vec(i);
	z     = zeta_vec(i);

	if flag_comp == 0
        % Element stiffness matrix
        [K_element] = Ke_2D(p,Hooke,x,y);

    elseif flag_comp == 1 || flag_comp == 2
        % Element stiffness and mass matrices
		[K_element,M_element] = KMe_2D(p,Hooke,dens,x,y,mass_matrix);
        % Damping matrix - Rayleigh hypothesis
        if z == 0
            R_element = sparse(ndof,ndof);
        else
            R_element = alpha_Rayleigh*M_element + beta_Rayleigh*K_element;	
        end
	end
    
    if flag_comp == 2
        % Critical time step of the central difference method from the
        % largest eigenvalue of each element, dt = 2/omega_max (a safety
        % factor of 0.99). The element estimate bounds the global maximum
        % eigenfrequency from above, so the minimum over the elements is
        % conservative. With damping, omega_max follows from the quadratic
        % eigenvalue problem in state-space form.
        if z ~= 0 && (alpha_Rayleigh ~= 0 || beta_Rayleigh~= 0)
            try
                A = [-speye(ndof)      , sparse(ndof,ndof); ...
                      sparse(ndof,ndof), K_element];
                B = [sparse(ndof,ndof), speye(ndof); ...
                     M_element        , R_element];
                omegaMax = eigs(A,-B,1,'lm',opts);
                dt_curr = 0.99*(2/abs(omegaMax));
            catch
                lambdaMax = eigs(K_element,M_element,1,'lm',opts);
                dt_curr = 0.99*(2/sqrt(abs(lambdaMax)));
            end
        else
            lambdaMax = eigs(K_element,M_element,1,'lm',opts);
            dt_curr = 0.99*(2/sqrt(abs(lambdaMax)));
        end

        % Save time increment
        dt_vec(i) = dt_curr;
    end
    
	% Save the matrices and topology data into vectors
	topo   = repmat(Topo(i,:),ndof,1);
	Ke{i}  = reshape(K_element,[],1);
	ROW{i} = reshape(topo',[],1);
	COL{i} = reshape(topo,[],1);
	if flag_comp == 1 || flag_comp == 2
		Me{i} = reshape(M_element,[],1);
		Re{i} = reshape(R_element,[],1);
	end
end

% Smallest time increment (based on highest elemental eigenfrequency)
dt = min(dt_vec);

% Element data (element-by-element time integration)
if nargout > 5
    ElemData.Ke = cellfun(@(v) reshape(v,ndof,ndof),Ke,'UniformOutput',false);
    if flag_comp == 1 || flag_comp == 2
        ElemData.Me = cellfun(@(v) reshape(v,ndof,ndof),Me,'UniformOutput',false);
    else
        ElemData.Me = {};
    end
    ElemData.Topo   = double(Topo);
    ElemData.damped = zeta_vec ~= 0;
    ElemData.alpha  = alpha_Rayleigh;
    ElemData.beta   = beta_Rayleigh;
end

%%
disp('Assembly process')
% Concatenate all the cell arrays into a matrix
k   = cell2mat(Ke);
row = cell2mat(ROW);
col = cell2mat(COL);
if flag_comp == 1 || flag_comp == 2
	m = cell2mat(Me);
	r = cell2mat(Re);
end

% Determine the Matlab version
vers = version;
num = str2double(vers(1:3));
if num < 9.8
    % The sparse command can only take integer values as input since
    % Matlab 2020a (vers. 9.8)
    if isinteger(row)
        row = double(row);
        col = double(col);
    end
end

% Save the result vectors in sparse matrices
% Stiffness matrix
K = sparse(row,col,k,num_DoF,num_DoF);
% Mass and damping matrices
if flag_comp == 1 || flag_comp == 2
    M = sparse(row,col,m,num_DoF,num_DoF); 
    if alpha_Rayleigh == 0 && beta_Rayleigh == 0
        R = sparse(num_DoF,num_DoF);
    else
        R = sparse(row,col,r,num_DoF,num_DoF);
    end
else
    M = sparse(num_DoF,num_DoF);
    R = sparse(num_DoF,num_DoF);
end

% Time increment
if flag_comp ~= 2
    dt = [];
end

% Info
info_K = whos('K');
info_M = whos('M');
info_R = whos('R');

disp('Data:')
disp('Number of degrees of freedom')
disp(['ndof = ' num2str(num_DoF)])

disp('Memory requirements in [MB]')
disp(['Stiffness matrix: ' num2str(round((info_K.bytes)/1024/1024)) 'MB'])
disp(['Mass matrix: ' num2str(round(info_M.bytes/1024/1024)) 'MB'])
disp(['Damping matrix: ' num2str(round(info_R.bytes/1024/1024)) 'MB'])

% End the time measurement	
t_cpu = toc(t0);
disp(['tcpu = ',num2str(t_cpu),'s'])

LineDisp = repmat('-',1,75);
disp(LineDisp)

end