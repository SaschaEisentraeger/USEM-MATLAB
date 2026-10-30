function [t,dsp,vel,acc,tInterp,dspInterp,velInterp,accInterp,UAll] = ...
         PadeTimeIntegrator(porder,ns,dt,forceHistory,K0,M0,C0,F0,U0,V0,...
                            CDofs,pDOF,nInterp)
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by:
% Chongmin Song (12/2020)
%
% Edited by:
% + Sascha Eisentraeger (12/2020):
% - Add methods of orders 11 and 55
% - Add initial conditions
% - Adjust force interpolation (order of the force approx. is the 
%   identical to the order of the denominator/numerator polynomial)
% - Added interpolated results
% + Chongmin Song (01/2021):
% - LU-decomposition of the dynamic stiffness matrix for complex roots
%   (using the command lu is faster than the command decomposition for a
%   Hermitian matrix)
% + Sascha Eisentraeger (03/2021):
% - Add export of the entire displacement field for a few time steps
%-------------------------------------------------------------------------%
% Input:
% porder:       order of Pade expansion: = 11; 22; 33; 44; 55; (-1)
% ns:           number of time steps
% dt:           time steps size
% forceHistory: external loading
%               [NN x 2]: [time amplitude]; computed at a temporal grid
% 						  with a fine resolution (much less than dt); NN
% 						  values per column
%               [ns x p+1]: [time amp_1 amp_2 ... amp_(p+1)]; exact values
%                           of the force within the time steps; amp_i are
%                           computed at the GLL-points mapped to the
%                           interval [0,1]*dt
% K0:           stiffness matrix including the constrained DOFs
% M0:           mass matrix including the constrained DOFs
% C0:           damping matrix including the constrained DOFs
% F0:           prescribed nodal forces
% U0:           Initial displacements
% V0:           Initial velocities
% CDofs:        constrained degrees of freedom
% pDOF:         the DOF for output response (response at pDOF is stored for
%               output)
% nInterp:      Factor to increase the number of time steps in the output
%               by interpolation
%
% Output:
% t:            time vector
% dsp:          displacement vector
% vel:          velocity vector
% acc:          acceleration vector
% dspInterp: 	interpolated displacement vector (cubic spline interpol.)
% velInterp:   	interpolated velocity vector (cubic spline interpol.)
% accInterp:   	interpolated acceleration vector (cubic spline interpol.)
% UAll:			entire displacement field saved at a few time steps
%-------------------------------------------------------------------------%

% (De-)activate the waitbar
flag_wb = true;

% Save the full displacement field for a few time steps
nUAll = 0;

%% Initial conditions
NDof = size(K0,1);  % number of DOFs
if isempty(U0)
    U0 = zeros(NDof,1);
    nnzU0 = 0;
else
    % Number non-zero terms
    nnzU0 = nnz(U0);
end
if isempty(V0)
    V0 = zeros(NDof,1);
    nnzV0 = 0;
else
    % Number non-zero terms
    nnzV0 = nnz(V0);
end
% Number of non-zero terms in the initial conditions
nnzIC = nnzU0 + nnzV0;

%% find fixed and free DOFs
% displacement boundary condition
FDofs = (1:NDof);   % initializing unconstrained (free) DOFs
FDofs(CDofs) = [];  % remove constrained DOFs

%% exclude fixed DOFs in the analysis
dt2 = dt*dt;
% Stiffness, damping, and mass matrices
K = dt2*sparse(K0(FDofs,FDofs));
C = dt*sparse(C0(FDofs,FDofs));
M = sparse(M0(FDofs,FDofs));
% Force vector
F = dt2*sparse(F0(FDofs));
% Initial displacement
u0 = U0(FDofs);
% Initial velocity --> normalized with dt
v0 = dt*V0(FDofs);

% identify DOF for output in the reduced system
dsp = zeros(NDof,1);
dsp(FDofs) = 1:length(FDofs);
pDOF = dsp(pDOF);
pDOF(pDOF == 0) = [];

% This step is only necessary if the force history is given on a fine grid
% of points instead of the exact values
t = (0:ns)'*dt; % time at all time steps
if isempty(forceHistory)
    fn = zeros(size(t));
else
    % load factor at all time steps by linear interpolation
    fn = interp1(forceHistory(:,1),forceHistory(:,2),t,'linear',0);
end

% Polynomial variation of the force in each time step:
% Values are computed at the mid-point of the time increment
% The derivatives are computed with respect to the non-dimensional time s,
% where t = t_i + s*dt and s in [0,1] --> Df/Ds = Df/Dt*dt.
if nnz(fn) == 0
    % Set force values and derivatives to zero
    if porder == -1
        ft = fn;
    else
        ft = sparse(size(fn,1)-1,porder/11+1);
    end
else
    % Interpolate the force values and the derivatives
    switch porder
        case 11
            % Linear variation of the force (per time step)
            % GLL-points in the interval [-1,1]
            xi = lglnodes(1);
            xi = flip(xi);
            % Mapping to the interval [0,1]
            s  = 1/2*(xi + 1);
            a1 = (s-0.5).^[ 0 1];
            % Calculate the force vector of the time integration scheme:
            % The first column contains the value of the force at the
            % middle of the time step, while the second column holds the
            % first derivative of the force vector at the middle of the
            % time step.
            A1 = (a1'*a1)\a1';
            ft = (A1*[fn(1:ns),fn(2:ns+1)]')';
        case 22
            % Quadratic variation of the force (per time step)
            % GLL-points in the interval [-1,1]
            xi = lglnodes(2);
            xi = flip(xi);
            % Mapping to the interval [0,1]
            s  = 1/2*(xi + 1);
            a2 = (s-0.5).^[ 0 1 2];
            % Force vector within the time step
            if size(forceHistory,2) == 3
                % Exact force values
                f1 = forceHistory(:,3);
            else
                % Interpolate the force values within the time step
                DT1 = (s(2):1:ns)'*dt;
                f1 = interp1(forceHistory(:,1),forceHistory(:,2),DT1,'linear',0);
            end
            % Calculate the force vector of the time integration scheme:
            % The first column contains the value of the force at the
            % middle of the time step, while the second column holds the
            % first derivative of the force vector at the middle of the 
            % time step, and the third column is the second derivative of
            % the force vector at the middle of the time step. 
            % --> The forcing function is approximated by Lagrangian
            % polynomials based on GLL-points; these functions are used to
            % compute the values at the midpoint.
            A2 = (a2'*a2)\a2';
            ft = (A2*[fn(1:ns),f1(1:ns),fn(2:ns+1)]')';
        case 33
            % Cubic variation of the force (per time step)
            % GLL-points in the interval [-1,1]
            xi = lglnodes(3);
            xi = flip(xi);
            % Mapping to the interval [0,1]
            s  = 1/2*(xi + 1);
            a3 = (s-0.5).^[ 0 1 2 3];
            % Force vector within the time step
            if size(forceHistory,2) == 4
                % Exact force values
                f1 = forceHistory(:,3);
                f2 = forceHistory(:,4);
            else
                % Interpolate the force values within the time step
                DT1 = (s(2):1:ns)'*dt;
                DT2 = (s(3):1:ns)'*dt;
                f1 = interp1(forceHistory(:,1),forceHistory(:,2),DT1,'linear',0);
                f2 = interp1(forceHistory(:,1),forceHistory(:,2),DT2,'linear',0);
            end  
            % Calculate the force vector of the time integration scheme:
            % The first column contains the value of the force at the
            % middle of the time step, while the second column holds the
            % first derivative of the force vector at the middle of the 
            % time step, and the third column is the second derivative of
            % the force vector at the middle of the time step, ... 
            % --> The forcing function is approximated by Lagrangian
            % polynomials based on GLL-points; these functions are used to
            % compute the values at the midpoint.
            A3 = (a3'*a3)\a3';
            ft = (A3*[fn(1:ns),f1(1:ns),f2(1:ns),fn(2:ns+1)]')';
        case 44
            % Quartic variation of the force (per time step)
            % GLL-points in the interval [-1,1]
            xi = lglnodes(4);
            xi = flip(xi);
            % Mapping to the interval [0,1]
            s  = 1/2*(xi + 1);
            a4 = (s-0.5).^[ 0 1 2 3 4];
            % Force vector within the time step
            if size(forceHistory,2) == 5
                % Exact force values
                f1 = forceHistory(:,3);
                f2 = forceHistory(:,4);
                f3 = forceHistory(:,5);
            else
                % Interpolate the force values within the time step
                DT1 = (s(2):1:ns)'*dt;
                DT2 = (s(3):1:ns)'*dt;
                DT3 = (s(4):1:ns)'*dt;
                f1 = interp1(forceHistory(:,1),forceHistory(:,2),DT1,'linear',0);
                f2 = interp1(forceHistory(:,1),forceHistory(:,2),DT2,'linear',0);
                f3 = interp1(forceHistory(:,1),forceHistory(:,2),DT3,'linear',0);
            end  
            % Calculate the force vector of the time integration scheme:
            % The first column contains the value of the force at the
            % middle of the time step, while the second column holds the
            % first derivative of the force vector at the middle of the 
            % time step, and the third column is the second derivative of
            % the force vector at the middle of the time step, ... 
            % --> The forcing function is approximated by Lagrangian
            % polynomials based on GLL-points; these functions are used to
            % compute the values at the midpoint.
            A4 = (a4'*a4)\a4';
            ft = (A4*[fn(1:ns),f1(1:ns),f2(1:ns),f3(1:ns),fn(2:ns+1)]')';
        case 55
            % Quintic variation of the force (per time step)
            % GLL-points in the interval [-1,1]
            xi = lglnodes(5);
            xi = flip(xi);
            % Mapping to the interval [0,1]
            s  = 1/2*(xi + 1);
            a5 = (s-0.5).^[ 0 1 2 3 4 5];
            % Force vector within the time step
            if size(forceHistory,2) == 6
                % Exact force values
                f1 = forceHistory(:,3);
                f2 = forceHistory(:,4);
                f3 = forceHistory(:,5);
                f4 = forceHistory(:,6);
            else
                % Interpolate the force values within the time step
                DT1 = (s(2):1:ns)'*dt;
                DT2 = (s(3):1:ns)'*dt;
                DT3 = (s(4):1:ns)'*dt;
                DT4 = (s(5):1:ns)'*dt;
                f1 = interp1(forceHistory(:,1),forceHistory(:,2),DT1,'linear',0);
                f2 = interp1(forceHistory(:,1),forceHistory(:,2),DT2,'linear',0);
                f3 = interp1(forceHistory(:,1),forceHistory(:,2),DT3,'linear',0);
                f4 = interp1(forceHistory(:,1),forceHistory(:,2),DT4,'linear',0);
            end  
            % Calculate the force vector of the time integration scheme:
            % The first column contains the value of the force at the
            % middle of the time step, while the second column holds the
            % first derivative of the force vector at the middle of the 
            % time step, and the third column is the second derivative of
            % the force vector at the middle of the time step, ... 
            % --> The forcing function is approximated by Lagrangian
            % polynomials based on GLL-points; these functions are used to
            % compute the values at the midpoint.
            A5 = (a5'*a5)\a5';
            ft = (A5*[fn(1:ns),f1(1:ns),f2(1:ns),f3(1:ns),f4(1:ns),fn(2:ns+1)]')';
        case -1
            ft = fn;
        otherwise
            % Linear variation of the force (per time step)
            % GLL-points in the interval [-1,1]
            xi = lglnodes(1);
            xi = flip(xi);
            % Mapping to the interval [0,1]
            s  = 1/2*(xi + 1);
            a1 = (s-0.5).^[ 0 1];
            % Calculate the force vector of the time integration scheme:
            % The first column contains the value of the force at the
            % middle of the time step, while the second column holds the
            % first derivative of the force vector at the middle of the
            % time step.
            A1 = (a1'*a1)\a1';
            ft = (A1*[fn(1:ns),fn(2:ns+1)]')';
    end
end

% Command window output: Time step size
fprintf('Time step size used: %6.3e\n',dt);

switch porder
	case 11
        [dsp,vel,acc,UAll] = PadeSolver11(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
    case 22
        [dsp,vel,acc,UAll] = PadeSolver22(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
    case 33
        [dsp,vel,acc,UAll] = PadeSolver33(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
    case 44
        [dsp,vel,acc,UAll] = PadeSolver44(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
	case 55
        [dsp,vel,acc,UAll] = PadeSolver55(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
    case -1
        % Dirichlet boundary conditions have to applied otherwise the
        % analytical solution is not uniquely determined (for the
        % analytical solution the regularization of the mass matrix does
        % not help to prevent rigid body motions).
        % Exact solution based on the matrix exponential
        [dsp,vel,acc,UAll] = PadeSolverExact(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
    otherwise
        [dsp,vel,acc,UAll] = PadeSolver11(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
end

dsp = dsp';
vel = vel'/dt;
acc = acc'/(dt*dt);

if nargout > 4
    % Increase the resolution of the results by interpolation
    nInterp = ceil(nInterp);    % Ensure that nInterp is an integer value
    if nInterp > 1
        % New time vector
        nP = numel(t);
        lt = diff(t)/nInterp;
        nPNew = nInterp*(nP-1) + 1;
        tInterp = zeros(nPNew,1);
        tInterp(1:nInterp:nPNew) = t;
        for i = 2:nInterp
            tInterp(i:nInterp:nPNew) = t(1:nP-1) + (i-1)*lt;
        end
        % Interpolate displacement, velocity and acceleration vectors using
        % a cubic spline interpolation
        dspInterp = interp1(t,dsp,tInterp,'spline','extrap');
        velInterp = interp1(t,vel,tInterp,'spline','extrap');
        accInterp = interp1(t,acc,tInterp,'spline','extrap');
    else
        tInterp   = [];
        dspInterp = [];
        velInterp = [];
        accInterp = [];
    end
end

end