function [T,U] = PadeMethod(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,...
                 ActiveDoFs,ExciteParam,SignalType,dtmax,tsim,Order,...
                 rho_inf,nInterp)
%PADEMETHOD Pade-based time integration with the common calling convention.
%
%   [T,U] = PadeMethod(K,M,C,F,F_dt,U0,V0,NumDoF,DOFsave,ActiveDoFs,...
%                      ExciteParam,SignalType,dtmax,tsim,Order,rho_inf,...
%                      nInterp)
%   adapts the input of the time integration methods (see
%   TimeIntegrationMethods) to PadeTimeIntegrator (rho_inf = 1, no
%   numerical damping) and PadeTimeIntegrator_ND (rho_inf < 1).
%
%   Order    2, 4, 6, 8, 10: diagonal Pade approximation (1,1) ... (5,5)
%            of the matrix exponential (other values: 2)
%   rho_inf  spectral radius at infinite frequency, 1 (no damping) to 0
%   nInterp  factor by which the number of output instants is increased
%            (interpolation within the time steps)
%
%   Time step: dtmax(3) (prescribed implicit step), adjusted so that tsim
%   is an integer multiple of it. The load within a time step is
%   represented by a polynomial of the order of the Pade approximation
%   (ForcePadeTimeIntegration); see below for the loads that are taken
%   into account. NumDoF and ActiveDoFs are not needed (common signature).

% Time step size: prescribed implicit step
if numel(dtmax) == 3
    dt = dtmax(3);
else
    dt = dtmax(1);
end
ns = ceil(round(tsim/dt,3,'decimal'));
dt = tsim/ns;

% Code of the Pade approximation (11, 22, 33, 44, 55)
codes = [11 22 33 44 55];
k = find(Order == [2 4 6 8 10],1);
if isempty(k)
    order = 11;
else
    order = codes(k);
end

% Load amplitudes at the interpolation points of each time step
% (column 1: time, columns 2 to end: amplitudes)
forceHistory = ForcePadeTimeIntegration(order,SignalType,ExciteParam,ns,dt);

% Matrices stored as triangular matrices are symmetric matrices of which
% only one half is stored (the unsymmetric matrices of this code are full)
if (istriu(K) || istril(K)) && ~isscalar(K)
    K = K + K';
    M = M + M';
    C = C + C';
end

% The solvers integrate a single load vector F0 times the amplitudes of
% forceHistory:
% - no signal (SignalType 0): the constant load F with amplitude 1, as in
%   the other methods (F + F_dt*s with s = 0); the amplitudes of signal 0
%   are zero and are therefore replaced by 1.
% - otherwise: the time-dependent load F_dt.Neumann times s(t). A constant
%   part F and the lifting terms F_dt.K, F_dt.C, F_dt.M of prescribed
%   displacements are not taken into account (docs: chapter "Known
%   issues").
if SignalType == 0
    F0 = F;
    forceHistory(:,2:end) = 1;
else
    F0 = F_dt.Neumann;
end

% The DOFs with homogeneous Dirichlet conditions have been removed
if rho_inf == 1
    [T,U] = PadeTimeIntegrator(order,ns,dt,forceHistory,K,M,C,F0,U0,V0,...
            [],DOFsave,nInterp);
else
    % Numerical damping: the damped solvers expect the order (1..5)
    order = order - 10;
    if exist('PadeTimeIntegrator_ND_FORTRAN','file')
        % Fortran solver (not part of the repository)
        [T,U] = PadeTimeIntegrator_ND_FORTRAN(order,rho_inf,ns,dt,...
                forceHistory,K,M,C,F0,U0,V0,[],DOFsave,nInterp);
    else
        [T,U] = PadeTimeIntegrator_ND(order,rho_inf,ns,dt,forceHistory,...
                K,M,C,F0,U0,V0,[],DOFsave,nInterp);
    end
end

end
