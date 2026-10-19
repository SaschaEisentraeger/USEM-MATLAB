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
        [dsp,vel,acc,UAll] = TimeSolver11(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
    case 22
        [dsp,vel,acc,UAll] = TimeSolver22(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
    case 33
        [dsp,vel,acc,UAll] = TimeSolver33(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
    case 44
        [dsp,vel,acc,UAll] = TimeSolver44(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
	case 55
        [dsp,vel,acc,UAll] = TimeSolver55(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
    case -1
        % Dirichlet boundary conditions have to applied otherwise the
        % analytical solution is not uniquely determined (for the
        % analytical solution the regularization of the mass matrix does
        % not help to prevent rigid body motions).
        % Exact solution based on the matrix exponential
        [dsp,vel,acc,UAll] = TimeSolverEx(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
    otherwise
        [dsp,vel,acc,UAll] = TimeSolver11(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb);
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

%% Order (1,1)

function [dsp,vel,acc,UAll] = TimeSolver11(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb)
%
% number of time steps
ns = size(ft,1) + 1;

% Save the full displacement field for a few time steps
if nUAll > 0
    nSave = ceil(ns/nUAll);
    UAll = cell(nUAll,1);
else
	UAll = [];
end

% initializing variables storing response history for output
dsp = zeros(length(pDOF),ns); % displacements
vel = dsp;  % velocities
acc = dsp;  % accelerations

dM = decomposition(M);
Fb = (dM\F);

if nnzIC ~= 0
	x = [v0; u0];
    Ax1 = SolverPadeAx(dM,K,C,x);
	
	% Initial conditions
    dsp(:,1) = u0(pDOF);
    vel(:,1) = v0(pDOF);
	
	% Initial acceleration --> normalized with dt^2
	fn = ft(1,:)*[1 (-0.5).^(1:size(ft,2)-1)]'; % force at t=0
    acc(:,1) = Ax1(pDOF) + fn*Fb(pDOF);
else
	x = zeros(2*ndof,1);
    Ax1 = x;
end

r = 2;  % Roots of the denominator polynomial Q
Kd = sparse((r*r)*M + r*C + K);
dKd = decomposition(Kd);

%%
ndof = size(K,1);

if flag_wb
    % Waitbar
    h = waitbar(0,'2nd-order Pade-based time integrator (implicit):');

    % Less updates of the waitbar (avoid slowing the simulation down)
    nUpdate = ceil(ns/1000);
end

disp('Solve the equations of motion using the Pade-based time integrator');
disp('');
disp('---------------------------------------------------------------------');
disp('');

for it = 2:ns
	is = it - 1;
	
    b1 = 4*x(1:ndof) + 2*ft(is,1)*Fb;
	b2 = 4*x(ndof+1:2*ndof);

	ftmp = -(K*b2) + r*(M*b1);
	tmp = dKd\ftmp;
	xb = [tmp ; (tmp+b2)/r];
	
	x = -x + xb;
    
% 	Ax1 = SolverPadeAx(dM,K,C,x);

	Ax1 = -Ax1 + r*xb - [b1; b2];
    
    % store responses for output
    vel(:,it) = x(pDOF);
    dsp(:,it) = x(ndof+pDOF);
	fn = ft(is,1)+0.5*ft(is,2);
    acc(:,it) = Ax1(pDOF) + fn*Fb(pDOF);
	
	% save entire displacement field
	if nUAll > 0
        if mod(it,nSave) == 0
            c = it/nSave;
            UAll{c} = x(ndof+1:2*ndof);
		elseif it == ns
            UAll{nUAll} = x(ndof+1:2*ndof);
        end
	end
    
    if flag_wb
        % Waitbar
        if mod(it,nUpdate) == 0
            waitbar(it/ns,h,['2nd-order Pade-based method (implicit): ' ...
                    num2str(it/ns*100,'%5.1f') '%'])
        end
    end
    
end

if flag_wb
    % Delete the waitbar
    close(h)
end

end

%% Order (2,2)

function [dsp,vel,acc,UAll] = TimeSolver22(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb)
%
% number of time steps
ns = size(ft,1) + 1;

% Save the full displacement field for a few time steps
if nUAll > 0
    nSave = ceil(ns/nUAll);
    UAll = cell(nUAll,1);
else
	UAll = [];
end

% initializing variables storing response history for output
dsp = zeros(length(pDOF),ns); % displacements
vel = dsp;  % velocities
acc = dsp;  % accelerations

dM = decomposition(M);
Fb = (dM\F);

if nnzIC ~= 0
	x = [v0; u0];
    Ax1 = SolverPadeAx(dM,K,C,x);
	
	% Initial conditions
    dsp(:,1) = u0(pDOF);
    vel(:,1) = v0(pDOF);
	
	% Initial acceleration --> normalized with dt^2
	fn = ft(1,:)*[1 (-0.5).^(1:size(ft,2)-1)]'; % force at t=0
    acc(:,1) = Ax1(pDOF) + fn*Fb(pDOF);
else
	x = zeros(2*ndof,1);
    Ax1 = x;
end

r = complex(3,sqrt(3)); % Roots of the denominator polynomial Q
Kd = sparse((r*r)*M + r*C + K);
[L,U,p,q] = lu(Kd,'vector');
% KdH = sparse(Kd'*Kd);	% Make the complex matrix Hermitian
% dKdH = decomposition(KdH);

%%
ndof = size(K,1);

f  = [Fb; zeros(ndof,1)];
Af1 = SolverPadeAx(dM,K,C,f);
Pf = [ 12*f -Af1 f];

npt = size(ft,2);
pt = [1 0.5.^(1:npt-1)]';

if flag_wb
    % Waitbar
    h = waitbar(0,'4th-order Pade-based time integrator (implicit):');

    % Less updates of the waitbar (avoid slowing the simulation down)
    nUpdate = ceil(ns/1000);
end

disp('Solve the equations of motion using the Pade-based time integrator');
disp('');
disp('---------------------------------------------------------------------');
disp('');

for it = 2:ns
	is = it - 1;
	
    b = 12*Ax1 + Pf(:,1:npt)*ft(is,:)';
	
    ctmp = r*(M*b(1:ndof)) - K*b(ndof+1:2*ndof);
	
	tmp = U\(L\ctmp(q));
    tmp(p) = tmp;
%     tmp = dKdH\(Kd'*ctmp);

    y = [tmp ; (tmp+b(ndof+1:2*ndof))/r];
    xb = -imag(y)/imag(r);
    
    x = x + xb;
    
% 	Ax1 = SolverPadeAx(dM,K,C,x);
% 	Ax2 = SolverPadeAx(dM,K,C,Ax1);

    Ax1b = - imag(r*y)/imag(r);
    Ax1 = Ax1 + Ax1b;
% 	Ax2b = b - 12*xb + 6*Ax1b;
% 	Ax2 = Ax2 + Ax2b;
    
    % store responses for output
    vel(:,it) = x(pDOF);
    dsp(:,it) = x(ndof+pDOF);
	fn = ft(is,:)*pt;
    acc(:,it) = Ax1(pDOF) + fn*Fb(pDOF);
    
	% save entire displacement field
	if nUAll > 0
        if mod(it,nSave) == 0
            c = it/nSave;
            UAll{c} = x(ndof+1:2*ndof);
		elseif it == ns
            UAll{nUAll} = x(ndof+1:2*ndof);
        end
	end
	
    if flag_wb
        % Waitbar
        if mod(it,nUpdate) == 0
            waitbar(it/ns,h,['4th-order Pade-based method (implicit): ' ...
                    num2str(it/ns*100,'%5.1f') '%'])
        end
    end
    
end

if flag_wb
    % Delete the waitbar
    close(h)
end

end

%% Order (3,3)

function [dsp,vel,acc,UAll] = TimeSolver33(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb)
%
% number of time steps
ns = size(ft,1) + 1;

% Save the full displacement field for a few time steps
if nUAll > 0
    nSave = ceil(ns/nUAll);
    UAll = cell(nUAll,1);
else
	UAll = [];
end

% initializing variables storing response history for output
dsp = zeros(length(pDOF),ns); % displacements
vel = dsp;  % velocities
acc = dsp;  % accelerations

dM = decomposition(M);
Fb = (dM\F);

if nnzIC ~= 0
	x = [v0; u0];
    Ax1 = SolverPadeAx(dM,K,C,x);
    Ax2 = SolverPadeAx(dM,K,C,Ax1);
	
	% Initial conditions
    dsp(:,1) = u0(pDOF);
    vel(:,1) = v0(pDOF);
	
	% Initial acceleration --> normalized with dt^2
	fn = ft(1,:)*[1 (-0.5).^(1:size(ft,2)-1)]'; % force at t=0
    acc(:,1) = Ax1(pDOF) + fn*Fb(pDOF);
else
	x = zeros(2*ndof,1);
    Ax1 = x;
    Ax2 = x;
end

r1 = 4.644370709252171; % Roots of the denominator polynomial Q
Kd1 = sparse((r1*r1)*M + r1*C + K);
dKd1 = decomposition(Kd1);

r2 = complex(3.677814645373910,3.508761919567443); % Roots of the denominator polynomial Q
Kd2 = sparse((r2*r2)*M + r2*C + K);
[L,U,p,q] = lu(Kd2,'vector');
% Kd2H = Kd2'*Kd2;    % Make the complex matrix Hermitian
% dKd2H = decomposition(Kd2H);

%%
ndof = size(K,1);

f  = [Fb; zeros(ndof,1)];
Af1 = SolverPadeAx(dM,K,C,f);
Af2 = SolverPadeAx(dM,K,C,Af1);
Pf = [ (120*f+2*Af2) -10*Af1 (10*f+Af2/2) -3*Af1/2];

npt = size(ft,2);
pt = [1 0.5.^(1:npt-1)]';

if flag_wb
    % Waitbar
    h = waitbar(0,'6th-order Pade-based time integrator (implicit):');

    % Less updates of the waitbar (avoid slowing the simulation down)
    nUpdate = ceil(ns/1000);
end

disp('Solve the equations of motion using the Pade-based time integrator');
disp('');
disp('---------------------------------------------------------------------');
disp('');

for it = 2:ns
    is = it - 1;
    
    b = 240*x + 24*Ax2 + Pf(:,1:npt)*ft(is,:)';
    
    tmp = dKd1\(M*(r1*b(1:ndof)) - K*b(ndof+1:2*ndof));
    x1 = [tmp ; (tmp+b(ndof+1:2*ndof))/r1];
    
    ctmp = r2*(M*x1(1:ndof)) - K*x1(ndof+1:2*ndof);
	
    tmp = U\(L\ctmp(q));
    tmp(p) = tmp;
%     tmp = dKd2H\(Kd2'*ctmp);

    y = [tmp ; (tmp+x1(ndof+1:2*ndof))/r2];
    xb =  -imag(y)/imag(r2);
    
    x = -x + xb;
    
% 	Ax1 = SolverPadeAx(dM,K,C,x);
% 	Ax2 = SolverPadeAx(dM,K,C,Ax1);
% 	Ax3 = SolverPadeAx(dM,K,C,Ax2);
    
    Ax1b = -imag(r2*y)/imag(r2);
    Ax1 = -Ax1 + Ax1b;
    Ax2b = x1 - real(r2*conj(r2))*xb + 2*real(r2)*Ax1b;
    Ax2 = -Ax2 +  Ax2b;
% 	Ax3b = -b + 120*xb -60*Ax1b + 12*Ax2b;
% 	Ax3 = -Ax3 +  Ax3b;
    
    % store responses for output
    vel(:,it) = x(pDOF);
    dsp(:,it) = x(ndof+pDOF);
	fn = ft(is,:)*pt;
    acc(:,it) = Ax1(pDOF) + fn*Fb(pDOF);
    
	% save entire displacement field
	if nUAll > 0
        if mod(it,nSave) == 0
            c = it/nSave;
            UAll{c} = x(ndof+1:2*ndof);
		elseif it == ns
            UAll{nUAll} = x(ndof+1:2*ndof);
        end
	end
	
    if flag_wb
        % Waitbar
        if mod(it,nUpdate) == 0
            waitbar(it/ns,h,['6th-order Pade-based method (implicit): ' ...
                    num2str(it/ns*100,'%5.1f') '%'])
        end
    end
    
end

if flag_wb
    % Delete the waitbar
    close(h)
end

end

%% Order (4,4)

function [dsp,vel,acc,UAll] = TimeSolver44(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb)
%
% number of time steps
ns = size(ft,1) + 1;

% Save the full displacement field for a few time steps
if nUAll > 0
    nSave = ceil(ns/nUAll);
    UAll = cell(nUAll,1);
else
	UAll = [];
end

% initializing variables storing response history for output
dsp = zeros(length(pDOF),ns); % displacements
vel = dsp;  % velocities
acc = dsp;  % accelerations

dM = decomposition(M);
Fb = (dM\F);

if nnzIC ~= 0
	x = [v0; u0];
    Ax1 = SolverPadeAx(dM,K,C,x);
    Ax2 = SolverPadeAx(dM,K,C,Ax1);
    Ax3 = SolverPadeAx(dM,K,C,Ax2);
	
	% Initial conditions
    dsp(:,1) = u0(pDOF);
    vel(:,1) = v0(pDOF);
	
	% Initial acceleration --> normalized with dt^2
	fn = ft(1,:)*[1 (-0.5).^(1:size(ft,2)-1)]'; % force at t=0
    acc(:,1) = Ax1(pDOF) + fn*Fb(pDOF);
else
	x = zeros(2*ndof,1);
    Ax1 = x;
    Ax2 = x;
    Ax3 = x;
end

r1 = complex(4.207578794359259, 5.314836083713504); % Roots of the denominator polynomial Q
Kd1 = sparse((r1*r1)*M + r1*C + K);
[L1,U1,p1,q1] = lu(Kd1,'vector');
% Kd1H = Kd1'*Kd1;    % Make the complex matrix Hermitian
% dKd1H = decomposition(Kd1H);

r2 = complex(5.792421205640749, 1.734468257869007); % Roots of the denominator polynomial Q
Kd2 = sparse((r2*r2)*M + r2*C + K);
[L2,U2,p2,q2] = lu(Kd2,'vector');
% Kd2H = Kd2'*Kd2;    % Make the complex matrix Hermitian
% dKd2H = decomposition(Kd2H);

%%
ndof = size(K,1);

f  = [Fb; zeros(ndof,1)];
Af1 = SolverPadeAx(dM,K,C,f);
Af2 = SolverPadeAx(dM,K,C,Af1);
Af3 = SolverPadeAx(dM,K,C,Af2);
Pf = [(1680*f+40*Af2) -(140*Af1+Af3) (140*f+8*Af2) -(21*Af1+0.25*Af3) (21*f+3*Af2/2)];

npt = size(ft,2);
pt = [1 0.5.^(1:npt-1)]';

if flag_wb
    % Waitbar
    h = waitbar(0,'8th-order Pade-based time integrator (implicit):');

    % Less updates of the waitbar (avoid slowing the simulation down)
    nUpdate = ceil(ns/1000);
end

disp('Solve the equations of motion using the Pade-based time integrator');
disp('');
disp('---------------------------------------------------------------------');
disp('');

for it = 2:ns
    is = it - 1;
    
    b = 1680*Ax1 + 40*Ax3 + Pf(:,1:npt)*ft(is,:)';
    
    ctmp = r1*(M*b(1:ndof)) - K*b(ndof+1:2*ndof);
	
	tmp = U1\(L1\ctmp(q1));
    tmp(p1) = tmp;
%     tmp = dKd1H\(Kd1'*ctmp);
	
    y1 = [tmp ; (tmp+b(ndof+1:2*ndof))/r1];
    x1 =  -imag(y1)/imag(r1);
    
    ctmp = r2*(M*x1(1:ndof)) - K*x1(ndof+1:2*ndof);
	
    tmp = U2\(L2\ctmp(q2));
    tmp(p2) = tmp;
%     tmp = dKd2H\(Kd2'*ctmp);
	
    y2 = [tmp ; (tmp+x1(ndof+1:2*ndof))/r2];
    xb =  -imag(y2)/imag(r2);
    
	x = x + xb;
    
% 	Ax1 = SolverPadeAx(dM,K,C,x);
% 	Ax2 = SolverPadeAx(dM,K,C,Ax1);
% 	Ax3 = SolverPadeAx(dM,K,C,Ax2);
% 	Ax4 = SolverPadeAx(dM,K,C,Ax3);

	Ax1b = -imag(r2*y2)/imag(r2);
	Ax1 = Ax1 + Ax1b;
	Ax2b = x1 - real(r2*conj(r2))*xb + 2*real(r2)*Ax1b;
	Ax2 = Ax2 + Ax2b;
	Ax3b = real(conj(r1)*x1 - y1 - (r2*conj(r2)*Ax1b-2*real(r2)*Ax2b));
 	Ax3 = Ax3 + Ax3b;
% 	Ax4b = b - 1680*xb+840*Ax1b-180*Ax2b+20*Ax3b;
% 	Ax4 = Ax4 + Ax4b;
    
    % store responses for output
    vel(:,it) = x(pDOF);
    dsp(:,it) = x(ndof+pDOF);
	fn = ft(is,:)*pt;
    acc(:,it) = Ax1(pDOF) + fn*Fb(pDOF);
    
	% save entire displacement field
	if nUAll > 0
        if mod(it,nSave) == 0
            c = it/nSave;
            UAll{c} = x(ndof+1:2*ndof);
		elseif it == ns
            UAll{nUAll} = x(ndof+1:2*ndof);
        end
	end
	
    if flag_wb
        % Waitbar
        if mod(it,nUpdate) == 0
            waitbar(it/ns,h,['8th-order Pade-based method (implicit): ' ...
                    num2str(it/ns*100,'%5.1f') '%'])
        end
    end
    
end

if flag_wb
    % Delete the waitbar
    close(h)
end

end

%% Order (5,5)

function [dsp,vel,acc,UAll] = TimeSolver55(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb)
%
% number of time steps
ns = size(ft,1) + 1;

% Save the full displacement field for a few time steps
if nUAll > 0
    nSave = ceil(ns/nUAll);
    UAll = cell(nUAll,1);
else
	UAll = [];
end

% initializing variables storing response history for output
dsp = zeros(length(pDOF),ns); % displacements
vel = dsp;  % velocities
acc = dsp;  % accelerations

dM = decomposition(M);
Fb = (dM\F);

if nnzIC ~= 0
	x = [v0; u0];
    Ax1 = SolverPadeAx(dM,K,C,x);
    Ax2 = SolverPadeAx(dM,K,C,Ax1);
    Ax3 = SolverPadeAx(dM,K,C,Ax2);
	Ax4 = SolverPadeAx(dM,K,C,Ax3);
	
	% Initial conditions
    dsp(:,1) = u0(pDOF);
    vel(:,1) = v0(pDOF);
	
	% Initial acceleration --> normalized with dt^2
	fn = ft(1,:)*[1 (-0.5).^(1:size(ft,2)-1)]'; % force at t=0
    acc(:,1) = Ax1(pDOF) + fn*Fb(pDOF);
else
	x = zeros(2*ndof,1);
    Ax1 = x;
    Ax2 = x;
    Ax3 = x;
	Ax4 = x;
end

r1 = 7.293477190659283; % Roots of the denominator polynomial Q
Kd1 = sparse((r1*r1)*M + r1*C + K);
dKd1 = decomposition(Kd1);

r2 = complex(4.649348606363277, 7.142045840675953); % Roots of the denominator polynomial Q
Kd2 = sparse((r2*r2)*M + r2*C + K);
[L2,U2,p2,q2] = lu(Kd2,'vector');
% Kd2H = Kd2'*Kd2;    % Make the complex matrix Hermitian
% dKd2H = decomposition(Kd2H);

r3 = complex(6.703912798307075, 3.485322832366395); % Roots of the denominator polynomial Q
Kd3 = sparse((r3*r3)*M + r3*C + K);
[L3,U3,p3,q3] = lu(Kd3,'vector');
% Kd3H = Kd3'*Kd3;    % Make the complex matrix Hermitian
% dKd3H = decomposition(Kd3H);

%%
ndof = size(K,1);

f  = [Fb; zeros(ndof,1)];
Af1 = SolverPadeAx(dM,K,C,f);
Af2 = SolverPadeAx(dM,K,C,Af1);
Af3 = SolverPadeAx(dM,K,C,Af2);
Af4 = SolverPadeAx(dM,K,C,Af3);
Pf = [(30240*f+840*Af2+2*Af4) -(2520*Af1+28*Af3) (2520*f+154*Af2+Af4/2) ...
      -(378*Af1+6*Af3) (378*f+57/2*Af2+Af4/8) -(5/4*Af3+135/2*Af1)];

npt = size(ft,2);
pt = [1 0.5.^(1:npt-1)]';

if flag_wb
    % Waitbar
    h = waitbar(0,'10th-order Pade-based time integrator (implicit):');

    % Less updates of the waitbar (avoid slowing the simulation down)
    nUpdate = ceil(ns/1000);
end

disp('Solve the equations of motion using the Pade-based time integrator');
disp('');
disp('---------------------------------------------------------------------');
disp('');

for it = 2:ns
    is = it - 1;
    
	b = 60480*x + 6720*Ax2 + 60*Ax4 + Pf(:,1:npt)*ft(is,:)';
    
	tmp = dKd1\(M*(r1*b(1:ndof)) - K*b(ndof+1:2*ndof));
    x1 = [tmp ; (tmp+b(ndof+1:2*ndof))/r1];
	
	ctmp = r2*(M*x1(1:ndof)) - K*x1(ndof+1:2*ndof);
	
    tmp = U2\(L2\ctmp(q2));
    tmp(p2) = tmp;
%     tmp = dKd2H\(Kd2'*ctmp);
	
    y2 = [tmp ; (tmp+x1(ndof+1:2*ndof))/r2];
    x2 =  -imag(y2)/imag(r2);
	
	ctmp = r3*(M*x2(1:ndof)) - K*x2(ndof+1:2*ndof);
	
    tmp = U3\(L3\ctmp(q3));
    tmp(p3) = tmp;
%     tmp = dKd3H\(Kd3'*ctmp);
	
    y3 = [tmp ; (tmp+x2(ndof+1:2*ndof))/r3];
    xb =  -imag(y3)/imag(r3);
	
	x = -x + xb;
    
% 	Ax1 = SolverPadeAx(dM,K,C,x);
% 	Ax2 = SolverPadeAx(dM,K,C,Ax1);
% 	Ax3 = SolverPadeAx(dM,K,C,Ax2);
% 	Ax4 = SolverPadeAx(dM,K,C,Ax3);
% 	Ax5 = SolverPadeAx(dM,K,C,Ax4);
    
	Ax1b = -imag(r3*y3)/imag(r3);
	Ax1 = -Ax1 + Ax1b;
	Ax2b = x2 - real(r3*conj(r3))*xb + 2*real(r3)*Ax1b;
	Ax2 = -Ax2 + Ax2b;
	Ax3b = real(conj(r2)*x2 - y2 - (r3*conj(r3)*Ax1b-2*real(r3)*Ax2b));
 	Ax3 = -Ax3 + Ax3b;
    Ax4b = real(x1 - r2*conj(r2)*x2 + (...
                      2*real(r2)*r3*conj(r3)*Ax1b - ...
           (r3*conj(r3)+4*real(r2)*real(r3))*Ax2b + ...
                     (2*real(r2)+2*real(r3))*Ax3b));
    Ax4 = -Ax4 + Ax4b;
% 	Ax5b = b + 30240*xb-15120*Ax1b+3360*Ax2b-420*Ax3b+30*Ax4b;
% 	Ax5 = Ax5 + Ax5b;
	
	% store responses for output
    vel(:,it) = x(pDOF);
    dsp(:,it) = x(ndof+pDOF);
	fn = ft(is,:)*pt;
    acc(:,it) = Ax1(pDOF) + fn*Fb(pDOF);
	
	% save entire displacement field
	if nUAll > 0
        if mod(it,nSave) == 0
            c = it/nSave;
            UAll{c} = x(ndof+1:2*ndof);
		elseif it == ns
            UAll{nUAll} = x(ndof+1:2*ndof);
        end
	end
	
    if flag_wb
        % Waitbar
        if mod(it,nUpdate) == 0
            waitbar(it/ns,h,['10th-order Pade-based method (implicit): ' ...
                    num2str(it/ns*100,'%5.1f') '%'])
        end
    end
    
end

if flag_wb
    % Delete the waitbar
    close(h)
end

end

%% Exact for a piecewise linear excitation

function [dsp,vel,acc,UAll] = TimeSolverEx(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb)
% The solution is only exact if the external excitation follows a piecewise
% linear function. Otherwise the force term has to be adjusted to the
% correct function.

% number of time steps
ns = length(ft);

% Save the full displacement field for a few time steps
if nUAll > 0
    nSave = ceil(ns/nUAll);
    UAll = cell(nUAll,1);
else
	UAll = [];
end

% initializing variables storing response history for output
dsp = zeros(length(pDOF),length(ft)); % displacements
vel = dsp;  % velocities
acc = dsp;  % accelerations

dM = decomposition(M);

if nnzIC ~= 0
	x = [v0; u0];
    Ax = A*x;
	
	% Initial conditions
    dsp(:,1) = u0(pDOF);
    vel(:,1) = v0(pDOF);
	
	% Initial acceleration --> normalized with dt^2
	fn = ft(1,:)*[1 (-0.5).^(1:size(ft,2)-1)]'; % force at t=0
    acc(:,1) = Ax1(pDOF) + fn*Fb(pDOF);
else
	x = zeros(2*ndof,1);
end

ndof = size(K,1);

I = speye(ndof);
zm = spalloc(ndof,ndof,0);

k  = dM\K;
c  = dM\C;
A = [ -c, -k; I, zm];
dA = decomposition(A);

Fb = dM\F;
f  = [Fb; zeros(ndof,1)];

%%
I = speye(2*ndof);
disp('Compute the matrix exponential');
meth = 3;
if meth == 1
    % Original Matlab function to compute the matrix exponential
    P = expm(A);
elseif meth == 2
    % Fastest version with larger differences
    P = expm2(A,1e-14,1e-14);
elseif meth == 3
    % Faster than expm but slower than expm2; more accurate than expm2
    P = expm3(A);
end

P0 = dA\(P - I);
P1 = dA\(P0 - I) - 0.5*P0;

if nnzIC == 0
	x = zeros(2*ndof,1);
else
    x = [v0; u0];
end

if flag_wb
    % Waitbar
    h = waitbar(0,'Matrix exponential time integrator (exact):');

    % Less updates of the waitbar (avoid slowing the simulation down)
    nUpdate = ceil(ns/1000);
end

disp('Solve the equations of motion using the Matrix exponential time integrator');
disp('');
disp('---------------------------------------------------------------------');
disp('');

for it = 2:ns
    
    tmp = P*x + P0*(((ft(it)+ft(it-1))/2)*f);
    tmp = tmp + P1*((ft(it)-ft(it-1))*f);
    x   = tmp;
    Ax  = A*x;
    
    % store responses for output
    vel(:,it) = x(pDOF);
    dsp(:,it) = x(ndof+pDOF);
    acc(:,it) = Ax(pDOF) + ft(it)*Fb(pDOF);
    
	% save entire displacement field
	if nUAll > 0
        if mod(it,nSave) == 0
            c = it/nSave;
            UAll{c} = x(ndof+1:2*ndof);
		elseif it == ns
            UAll{nUAll} = x(ndof+1:2*ndof);
        end
	end
	
    if flag_wb
        % Waitbar
        if mod(it,nUpdate) == 0
            waitbar(it/ns,h,['Matrix exponential time integrator (exact): ' ...
                    num2str(it/ns*100,'%5.1f') '%'])
        end
    end
    
end

if flag_wb
    % Delete the waitbar
    close(h)
end

end