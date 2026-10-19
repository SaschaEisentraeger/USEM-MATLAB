%% Two-Dimensional Higher Order Finite Element Simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; 
close all; clc; format long
dbstop if error     % Enter debug mode if an error is encountered

% Delete all *.mat files in the Temp directory prior to the simulation
delete(['.' filesep 'subroutines' filesep 'Temp' filesep '*.mat'])

%%
Line= '-';
LineDisp = repmat(Line,1,75);
%--------------------------------------------------------------------------
disp([repmat(Line,1,28) 'Start of Simulation!' repmat(Line,1,27)])
disp('Spectral elements')
disp(LineDisp)

disp('Loading global settings')
disp(LineDisp)
%% Load all folders into the Matlab search path
path1 = genpath(['.' filesep 'subroutines']);
path2 = genpath(['.' filesep 'InputFilesAbaqus']);
addpath(path1,path2);

%% Parallel set-up
%-------------------------------------------------------------------------%
% To adjust the preferences of the parallel computing toolbox; type 
% 'preferences' in the command window and go to Parallel Computing Toolbox
% --> Cluster Profile Manager
%-------------------------------------------------------------------------%
% Without the Parallel Computing Toolbox the parfor loops run serially
if ~isempty(ver('parallel')) && license('test','Distrib_Computing_Toolbox')
    % Check existing parallel pool
    PoolData = gcp('nocreate');
    if isempty(PoolData)
        % Number of physical cores available
        NCores = feature('numcores');
        % Control maximum number of computational threads
        maxNumCompThreads(NCores);
        % Start the parallel pool
        poolobj = parpool(NCores);
    end
end

%% Abaqus Input File
% filename = 'PatchTestLeeQuadratic_undistorted.inp';
% filename = 'PatchTestLeeQuadratic_angular.inp';
% filename = 'PatchTestLeeQuadratic_midedge.inp';
% filename = 'PatchTestLeeCubic_undistorted.inp';
% filename = 'PatchTestLeeCubic_angular.inp';
% filename = 'PatchTestLeeCubic_midedge.inp';
filename = 'PatchTestLeeCubic_Rajendran_adj.inp';

%% Set-up of the finite element model
disp('Set simulation parameters')
%% General settings
% Type of the conducted analysis
% Static_analysis
% Dynamic_analysis
% Modal_analysis
% FrequencyDomain_analysis
analysis_type = 'Static_analysis';

% Save entire workspace
saveWS = true;

% Plot settings
plot_mesh   = false;
plot_disp   = true;
plot_stress = false;
plot_strain = false;

% 2D analysis
% 1: plane stress assumption
% 2: plane strain assumption
antype = 1;

% Polynomial degree for the spectral shape functions
p.pxi  = 3;
p.peta = 3;

% Nodal distribution for the SEM
% ndist = 1: Equidistant nodes
% ndist = 2: GLL-nodes (Gauss-Lobatto-Legendre)
% ndist = 3: GLC-nodes (Gauss-Lobatto-Chebychev)
p.ndist = 2;

%% Settings for dynamic analyses
% Mass matrix formulation
% 1: Consistent mass matrix
% 2: Lumped mass matrix (nodal quadrature method)
% 3: Lumped mass matrix (row sum lumping method)
% 4: Lumped mass matrix (HRZ-method; diagonal scaling)
mass_matrix = 1;

% Set a flag if the individual loads should be considered as separate load
% cases, i.e., one right-hand side vector for each load. All loads can be
% solved simultaneously by using a force matrix instead of a vector. Use
% mRHS (multiple right-hand side) codes.
flag_mRHS = true;

% Damping matrix: Rayleigh hypothesis
% Vibration modes with prescribed damping ratio
omega_vec = [0,0];
% Damping ratio
zeta_vec = [0,0];
% Visualization of the damping ratio
vis = false;
% Damping parameter
[alpha_Rayleigh,beta_Rayleigh] = RayleighDamping(omega_vec,zeta_vec,vis);

% Options for the eigenvalue solver:
% Number of eigenvalues to be computed (modal analysis)
numEigen = 200;
% Convergence tolerance
EVP_Solver.tol = 1e-14;                 % Default value: 1e-14
% Maximum number of algorithm iterations
EVP_Solver.maxit = 500;             	% Default value: 300
% Maximum size of Krylov subspace (if p == num_DoF the command eig is
% used to solve the generalized eigenvalue problem)
EVP_Solver.p = max(2*numEigen,20);  	% Default value: max(2*numEigen,20)
% Initial start vector
EVP_Solver.v0 = false;                  % Default value: false
% Diagnostic information
EVP_Solver.disp = false;             	% Default value: false
% Type of eigenvalues
EVP_Solver.type = 'smallestabs';   	    % Default value: 'smallestabs' / 'smallestreal'

% Parameter for the dynamic analysis - force excitation
% SignalType = 0: No time-dependent signal
% SignalType = 1: Sine-burst modulated with a Hann window
% SignalType = 2: Sine-burst windowed by an exponential function
% SignalType = 3: Linear chirp
% SignalType = 4: Harmonic signal
% SignalType = 5: Piecewise linear signal
% SignalType = 6: Ricker wavelet
% SignalType = 7: Gaussian impulses
SignalType = 1; % excitation signal

% Default settings for wave propagation analyses - force excitation
tE = 2e-5;	% simulation time
if SignalType == 0
    ExciteParam = tE;
elseif SignalType == 1
    fex = 5e6;   	% excitation frequency in [Hz]
    nex = 3;     	% number of periods
    ExciteParam = [fex, nex];
elseif SignalType == 2
    fex = 0;        % excitation frequency in [Hz]
    ExciteParam = fex;
elseif SignalType == 3
    fch1 = 0;     	% minimum frequency in [Hz]
    fch2 = 0;     	% maximum frequency in [Hz]
    tch  = 0;      	% time duration of linear chirp
    ExciteParam = [fch1,fch2,tch];
elseif SignalType == 4
    a    = 0;       % amplitude of the cosine function
    b    = 0;      	% amplitude of the sine function
    fcos = 0;     	% frequency of the cosine function in [Hz]
    fsin = 0;     	% frequency of the sine function in [Hz]
    ExciteParam = [a, b, fcos, fsin];
elseif SignalType == 5
    % Amplitude tableau for a piecewise linear function
    TS = [0,0;0.005,1;0.015,1;0.02,0;0.04,0;0.045,1;0.055,1;0.06,0;...
          0.08,0;0.085,1;0.095,1;0.1,0;0.12,0;0.125,1;0.135,1;0.14,0];
    ExciteParam = [TS;tE,0];
elseif SignalType == 6
    fex  = 200;   	% excitation frequency in [Hz]
    ExciteParam = fex;
elseif SignalType == 7
    a = 1;          % Amplitude
    b = 0.01;       % Position of the center of the peak for the first
                    % impulse in [s]
    c = 0.002;      % Standard deviation (controls the width of the
                    % impulse) --> only 1% of the amplitude for t = b +-3c
    d = 0.04;       % Distance between peaks
    n = 4;          % Number of peaks
    ExciteParam = [a, b, c, d, n];
end

% Time increment estimate based on the external excitation
if strcmp(analysis_type,'Dynamic_analysis')
    % Maximum allowed time increment depending on the center frequency of
    % the excitation (at least 20 sampling points for the center frequency)
    [f_max] = FFT_Signal(ExciteParam,SignalType);
    dtf = 0.05/f_max;
end

% Save results 
% if save = 1 all variables are saved to the result file
% if save = 0 only the nodes of interest (NoI) are saved (they are marked
% in the ABAQUS input file by fixing the degree-of-freedom 6 - rotation
% about the z-axis)
save_var = 0;

% Default time integration settings
% Time integration schemes
% Explicit
% 101: Central difference method
% 102: Noh-Bathe method
% 103: Soares method
% 104: Mirbagehri method
% 105: Kim method
% 106: Zhang-Xing method 1
% 107: Zhang-Xing method 2
% 108: Explicit generalized alpha method
% 109: Adaptive hybrid explicit method - Soares
% Explicit higher order
% 201: Kim method: order = 3,4
% 202: Two-stage Kim method: order = 4
% Implicit
% 301: Newmark method
% 302: HHT-alpha method
% 303: Generalized alpha method
% 304: Bathe method
% 305: beta1-beta2-Bathe method
% 306: rho-Bathe method
% 307: Soares method
% 308: Kim method
% 309: Kim method (collocation)
% Implicit higher order
% 401: Kim method (Lagrange): order = 1,2,3,4,5
% 402: Soares method
% 403: Pade-based method: order = 2,4,6,8,10
% Highly precise time integration method
% 501: HPTIM
% Parallel explicit dynamics (High Performance Computing)
% 601: Central difference method
TimeIntegration.Method = 101;
TimeIntegration.Order = 2;      % Only applicable to high order methods
TimeIntegration.rho = 1;        % Spectral radius of the amplification
                                % matrix (only applicable to method 403)
                                % 1: no damping; 0: maximal damping
% Time step size (if a special time step should be prescribed for implicit
% methods) --> if a value is given it is used for the analysis
TimeIntegration.dtImp = 5e-9;   % []
% Interpolate the results for high-order time integration schemes
TimeIntegration.nInterp = 1;    % factor by which the number of steps is
                                % increased
% Save the solver input data
TimeIntegration.SolverSave = true;

%% Frequency domain analysis
% Harmonic excitation frequency
fHarmonic = 10;

%% Input file
%--------------------------------------------------------------------------
disp('Read ABAQUS *.inp-file')
%%
% Parts of the file name
[~,name,ext] = fileparts(filename);
% Since ABAQUS 6.12 element and node set names can be used several times in
% the different parts of the input file. This causes problems since this
% code relies on different names for different sets. Therefore the
% following routine renames the sets.
Adjust_Set_names_input_file_ABAQUS(filename);
%%
% Adjusted input file is imported and thereafter deleted.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read elements and nodes (topology)
t0 = tic;
disp('Read mesh')
[nodes, topology, element_type, ~, ~] =  ...
 Read_mesh_ABAQUS([name '_adjusted' ext]);
toc(t0)
% Read material data
t0 = tic;
disp('Read material data')
[material_names, material_elset, material_properties] = ...
 Read_material_ABAQUS([name '_adjusted' ext]);
toc(t0)
% Read step definition
t0 = tic;
disp('Read step definition')
[step_name, step_type, t_end, numSteps] = ...
 Read_step_ABAQUS([name '_adjusted' ext]);
toc(t0)
% Read loads
t0 = tic;
disp('Read loads')
[cload, dsload, body_force, gravity] = ...
 Read_load_ABAQUS([name '_adjusted' ext]);
toc(t0)
% Read boundary conditions
t0 = tic;
disp('Read boundary conditions')
[boundary_conditions] = ...
 Read_boundary_condition_ABAQUS([name '_adjusted' ext]);
toc(t0)

delete(['.' filesep filename(1:end-4) '_adjusted.inp']);
disp(LineDisp)

xCoor = nodes(:,1);
yCoor = nodes(:,2);
Topology_Node = topology;

disp(LineDisp)
% Check if the local coordinate systems have the same orientation, if that
% is not the case pxi and peta have to be assigned the same value (it might
% also be possible to rearrange the element topology)
if p.pxi ~= p.peta
    % Determine the edge topology and connectivity
    [Topology_Edge,NumEdges] = Topology_Edge_2D(Topology_Node);
    [EdgeElemList,LocalEdgeNum] = ConnectivityEdgeElement(Topology_Edge);
    for i = 1:NumEdges
        if numel(EdgeElemList{i}) == 2
            r1 = mod(LocalEdgeNum{i}(1),2);
            r2 = mod(LocalEdgeNum{i}(2),2);
            % Both edges have to have an odd or an even number only then
            % they belong to the same local coordinate direction, otherwise
            % the mismatch in the polynomial degree must be adjusted. Here,
            % the minimum polynomial degree is chosen.
            if (r1 ~= 0 && r2 == 0) || (r1 == 0 && r2 ~= 0)
                psel = min([p.pxi,p.peta]);
                p.pxi  = psel;
                p.peta = psel;
                disp([newline 'Local coordinate systems of adjacent' ...
                     ' elements are not aligned: pxi=peta=' ...
                     num2str(psel) newline])
                disp(LineDisp)
                break;
            end
        end
    end
end

% Polynomial degree of the original data/input discretization
if numel(Topology_Node(1,:)) == 4
    p_org = 1;
else
    p_org = 2;
end

%%
% Adjust topology data, nodes and boundary conditions for higher-order
% spectral elements (tensor product space)
disp('Adjust topology (connectivity) data for the high order element.')
% Adjust the nodal coordinates and the topology data for arbitrary elements
t0 = tic;
[Topology_Node,xCoor,yCoor] = pSpectralElement_NodesTopology(...
                              Topology_Node,xCoor,yCoor,p,p_org);
toc(t0)
disp('Adjust Dirichlet and Neumann boundary condition data for the high order')
disp('element.')
% Adjust the boundary conditions (Dirichlet and Neumann)
t0 = tic;
[dsload,boundary_conditions] = pSpectralElement_BCs(dsload,...
                               boundary_conditions,Topology_Node,p,p_org);
toc(t0)

%%
disp(['Loading of ' filename ' is finished!'])
disp(LineDisp)

%% Display FE-grid
%--------------------------------------------------------------------------
if plot_mesh == 1
    disp('Visualize the finite element discretization')
    t0 = tic;
    flag_BCplot = 1;    % Plot the Dirichlet and Neumann BCs
    FE_Mesh(p,xCoor,yCoor,Topology_Node,cload,dsload,body_force,gravity,...
            boundary_conditions,flag_BCplot)
    toc(t0)
end

%% Start simulation
%--------------------------------------------------------------------------

if (strcmp(analysis_type,'Static_analysis'))
	call_Static
elseif (strcmp(analysis_type,'Dynamic_analysis'))
	call_Dynamic
elseif (strcmp(analysis_type,'Modal_analysis'))
    call_Modal
elseif (strcmp(analysis_type,'FrequencyDomain_analysis'))
    call_FrequencyDomain
end

%% 
disp(LineDisp)
disp([repmat(Line,1,28) 'End of Simulation!' repmat(Line,1,29)])
disp(LineDisp)