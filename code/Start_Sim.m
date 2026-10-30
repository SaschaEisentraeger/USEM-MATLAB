%% USEM: two-dimensional unsymmetric spectral element simulation
%
% Main script. It sets the parameters of a simulation and runs it:
%   1. SetupPaths       adds the code to the MATLAB path
%   2. DefaultSettings  all parameters with their default values
%   3. (this script)    the parameters of this simulation
%   4. RunSimulation    input file -> analysis -> results and plots
%
% Only the parameters that differ from the defaults need to be set here;
% DefaultSettings documents every parameter. After the run, the workspace
% contains
%   model     mesh, materials, loads and boundary conditions (ReadModel)
%   results   results of the analysis (see 'help StaticAnalysis', ...)
%   settings  the parameters used
%
% Other scripts and studies can call RunSimulation the same way.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all; clc; format long
dbstop if error     % Enter debug mode if an error is encountered

addpath(fileparts(mfilename('fullpath')));
SetupPaths
settings = DefaultSettings();

%% ABAQUS input file (code/InputFilesAbaqus)
% settings.filename = 'PatchTestLeeQuadratic_undistorted.inp';
% settings.filename = 'PatchTestLeeQuadratic_angular.inp';
% settings.filename = 'PatchTestLeeQuadratic_midedge.inp';
% settings.filename = 'PatchTestLeeCubic_undistorted.inp';
% settings.filename = 'PatchTestLeeCubic_angular.inp';
% settings.filename = 'PatchTestLeeCubic_midedge.inp';
settings.filename = 'PatchTestLeeCubic_Rajendran_adj.inp';

%% Analysis
% 'Static_analysis', 'Dynamic_analysis', 'Modal_analysis',
% 'FrequencyDomain_analysis'
settings.analysis_type = 'Static_analysis';

% 1: plane stress, 2: plane strain
settings.antype = 1;

%% Spectral element
% Polynomial degree in xi and eta
settings.p.pxi  = 3;
settings.p.peta = 3;
% Nodal distribution: 1 equidistant, 2 GLL, 3 GLC
settings.p.ndist = 2;
% Basis of the metric Vandermonde matrix: 'lagrange' (default),
% 'legendre', 'chebyshev', 'monomial'
settings.p.metricBasis = 'lagrange';
% Trial space of the metric shape functions: 1 Q_p in skew coordinates
% (default), 2 hybrid space
settings.p.metricMethod = 1;

%% Dynamic analyses
% Mass matrix: 1 consistent, 2 GLL nodal quadrature, 3 row sum, 4 HRZ
settings.mass_matrix = 1;
% Rayleigh damping: angular frequencies [rad/s] and damping ratios
settings.omega_vec = [0,0];
settings.zeta_vec  = [0,0];

% Modal analysis: number of eigenvalues (eigs options: EVP_Solver)
settings.numEigen = 200;

% Frequency domain analysis: excitation frequency [Hz]
settings.fHarmonic = 10;

% Transient analysis: excitation signal (0 none, 1 sine burst, ... see
% DefaultSettings) and its parameters ([] = DefaultExcitationParameters)
settings.SignalType  = 1;
settings.ExciteParam = [];
% One right-hand side per load (central difference or HHT-alpha only)
settings.flag_mRHS = true;
% Saved DOFs: 0 nodes of interest, 1 all
settings.save_var = 0;
% Time integration method (codes: TimeIntegrationMethods) and options
settings.TimeIntegration.Method = 101;
settings.TimeIntegration.Order  = 2;
settings.TimeIntegration.rho    = 1;
settings.TimeIntegration.dtImp  = 5e-9;

%% Output
settings.saveWS      = true;
settings.plot_mesh   = false;
settings.plot_disp   = true;
settings.plot_stress = false;
settings.plot_strain = false;

%% Run
[model,results,settings] = RunSimulation(settings);
