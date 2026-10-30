function settings = DefaultSettings()
%DEFAULTSETTINGS Default parameters of a simulation.
%
%   settings = DefaultSettings() returns a struct with every parameter that
%   RunSimulation reads. Start_Sim.m overrides some of them; any script or
%   study can do the same:
%
%       settings = DefaultSettings();
%       settings.filename      = 'PatchTestLeeCubic_angular.inp';
%       settings.analysis_type = 'Modal_analysis';
%       [model, results] = RunSimulation(settings);
%
%   The parameter table of the documentation (chapter "Getting started")
%   lists all fields.

%% Model and analysis
% ABAQUS input file: name of a file on the MATLAB path (e.g. in
% code/InputFilesAbaqus) or a full path
settings.filename = 'PatchTestLeeCubic_Rajendran_adj.inp';

% Type of the analysis:
% 'Static_analysis', 'Dynamic_analysis', 'Modal_analysis',
% 'FrequencyDomain_analysis'
settings.analysis_type = 'Static_analysis';

% 2D analysis: 1 plane stress, 2 plane strain
settings.antype = 1;

%% Spectral element
% Polynomial degree of the shape functions in xi and eta. If the local
% coordinate systems of adjacent elements are not aligned, both are set to
% the smaller value (AlignPolynomialDegrees).
settings.p.pxi  = 3;
settings.p.peta = 3;

% Nodal distribution: 1 equidistant, 2 Gauss-Lobatto-Legendre (GLL),
% 3 Gauss-Lobatto-Chebyshev (GLC)
settings.p.ndist = 2;

% Basis of the generalized Vandermonde matrix of the metric (trial) shape
% functions; the shape functions are the same for every basis, only the
% round-off differs (docs: chapter "Basis of the metric Vandermonde
% matrix"): 'lagrange' (default), 'legendre', 'chebyshev', 'monomial'
settings.p.metricBasis = 'lagrange';

% Trial space of the metric shape functions:
% 1 Q_p in corotational skew coordinates (default)
% 2 hybrid space P_p(x,y) + parametric complement
settings.p.metricMethod = 1;

% Singular metric interpolation matrix (nodes not unisolvent, e.g. two
% coincident nodes on a strongly distorted element): 0 stops with an
% error; a value > 0 moves one node of such an element by this fraction of
% the element size, e.g. 1e-3 (ShiftSingularNodes)
settings.p.nodeShift = 0;

%% Mass and damping
% Mass matrix: 1 consistent, 2 lumped (GLL nodal quadrature),
% 3 lumped (row sum), 4 lumped (HRZ, diagonal scaling)
settings.mass_matrix = 1;

% Rayleigh damping R = alpha*M + beta*K, fitted to the damping ratios
% zeta_vec at the angular frequencies omega_vec [rad/s] (RayleighDamping).
% omega_vec = [0,0]: no damping. Only elements whose material has a
% damping entry are damped.
settings.omega_vec = [0,0];
settings.zeta_vec  = [0,0];
% Plot the resulting damping ratio over the frequency
settings.plot_damping = false;

%% Modal analysis
% Number of eigenvalues
settings.numEigen = 200;
% Options of eigs (see 'help eigs'). EVP_Solver.p = [] selects
% max(2*numEigen,20); if 2*numEigen exceeds the number of active DOFs, all
% eigenvalues are computed with eig.
settings.EVP_Solver.tol   = 1e-14;
settings.EVP_Solver.maxit = 500;
settings.EVP_Solver.p     = [];
settings.EVP_Solver.v0    = false;      % true: random start vector
settings.EVP_Solver.disp  = false;
settings.EVP_Solver.type  = 'smallestabs';

%% Frequency domain analysis
% Harmonic excitation frequency [Hz]
settings.fHarmonic = 10;

%% Dynamic (transient) analysis
% Time signal s(t) of the time-dependent loads (ExcitationSignal):
% 0 none, 1 sine burst (Hann window), 2 burst (exponential window),
% 3 linear chirp, 4 harmonic, 5 piecewise linear, 6 Ricker wavelet,
% 7 Gaussian impulses
settings.SignalType = 1;
% Parameters of the signal; [] selects the defaults of
% DefaultExcitationParameters(SignalType, tE)
settings.ExciteParam = [];
% Duration used by the default parameters of some signals [s]
settings.tE = 2e-5;

% One right-hand side per load case (columns of F). The transient solver
% then uses the central difference method (Method < 300) or the HHT-alpha
% method (Method >= 300), whatever Method says.
settings.flag_mRHS = true;

% Output of the transient analysis:
% 0 only the nodes of interest (fixed DOF 6 in the ABAQUS input file; the
%   first three active nodes if none are marked), 1 all DOFs
settings.save_var = 0;

% Time integration (codes: TimeIntegrationMethods, docs: chapter
% "Analyses")
settings.TimeIntegration.Method  = 101;
settings.TimeIntegration.Order   = 2;   % higher-order methods only
settings.TimeIntegration.rho     = 1;   % spectral radius (403): 1 no
                                        % damping, 0 maximal damping
settings.TimeIntegration.dtImp   = 5e-9;% prescribed time step of the
                                        % implicit methods [s]
settings.TimeIntegration.nInterp = 1;   % interpolation factor (403)
% Keep the solver input (subroutines/Temp/Solver_input.mat) to rerun
% Solver_Transient without the preprocessing
settings.TimeIntegration.SolverSave = true;

% Simulation time [s]; [] uses the step time of the ABAQUS input file
settings.t_end = [];

%% Output
% Save model, settings and results to ResultsFolder (SaveResults)
settings.saveWS = true;

% Folders for the results and the solver input of a transient analysis
% (default: code/subroutines/Results and code/subroutines/Temp)
subroutines = fileparts(fileparts(mfilename('fullpath')));
settings.ResultsFolder = fullfile(subroutines,'Results');
settings.TempFolder    = fullfile(subroutines,'Temp');

% Plots
settings.plot_mesh   = false;
settings.plot_disp   = true;
settings.plot_stress = false;
settings.plot_strain = false;

% Print progress information
settings.verbose = true;

% Open a parallel pool for the element loops (parfor) if the Parallel
% Computing Toolbox is available; without it the loops run serially
settings.parallel = true;

end
