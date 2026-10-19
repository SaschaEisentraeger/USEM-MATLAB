function [T,U,t_cpu] = Solver_Transient()
%-------------------------------------------------------------------------%
% Reference for different time integration methods
% Bathe, K. J. Finite-Elemente-Methoden, Chap. 9, pp. 785ff, 
% Springer Verlag, 2002
%-------------------------------------------------------------------------%
load(['.' filesep 'subroutines' filesep 'Temp' filesep ...
      'Solver_input.mat'],'ActiveDoF','dte','dtf','ExciteParam','NOI',...
	  'num_DoF','save_var','SignalType','t_end','TimeIntegration',...
	  'Topology_Node','U0','V0','flag_mRHS')

% Load matrices
load(['.' filesep 'subroutines' filesep 'Temp' filesep 'K.mat'],'K')
load(['.' filesep 'subroutines' filesep 'Temp' filesep 'M.mat'],'M')
load(['.' filesep 'subroutines' filesep 'Temp' filesep 'R.mat'],'R')
load(['.' filesep 'subroutines' filesep 'Temp' filesep 'F.mat'],'F')
load(['.' filesep 'subroutines' filesep 'Temp' filesep 'F_dt.mat'],'F_dt')

% Boundary conditions
if size(K,1) ~= numel(ActiveDoF)
    K            = K(ActiveDoF,ActiveDoF);
    M            = M(ActiveDoF,ActiveDoF);
    R            = R(ActiveDoF,ActiveDoF);
    F            = F(ActiveDoF,:);
    F_dt.Neumann = F_dt.Neumann(ActiveDoF,:);
    F_dt.K       = F_dt.K(ActiveDoF,:);
    F_dt.C       = F_dt.C(ActiveDoF,:);
    F_dt.M       = F_dt.M(ActiveDoF,:);
    U0           = U0(ActiveDoF,:);
    V0           = V0(ActiveDoF,:);
end

% Save the transient data only for the degrees of freedom corresponding to
% the selected points of interest
dof_poi = zeros(2*numel(NOI),1);
dof_poi(1:2:end) = 2*NOI-1;
dof_poi(2:2:end) = 2*NOI;

% Identify DOF for output in the reduced system
vecID = zeros(num_DoF,1);
vecID(ActiveDoF) = 1:length(ActiveDoF);
dof_poi = vecID(dof_poi);
dof_poi(dof_poi == 0) = [];

if save_var == 1
    % Save all degrees of freedom
    dof_poi = 1:numel(ActiveDoF);
elseif save_var == 0 && isempty(dof_poi)
    % If no node of interest has been specified, save DOFs 1 to 4
    dof_poi = [1 2 3 4];
    % Display a warning message in the command window
    warning(['No nodes of interest have been defined. Therefore, the ' ...
             'degrees of freedom DOF = [' num2str(ActiveDoF(1)) ','   ...
             num2str(ActiveDoF(2)) ',' num2str(ActiveDoF(3)) ','     ...
             num2str(ActiveDoF(4)) '] are saved.'])
end
 
disp('Determine the time step size');
% Time increments
% (1) dtf: 	   		based on the maximum frequency in the excitation signal
%                 	(FFT)
% (2) dte:     		based on the highest num. eigenfrequency of each
% 					element (conservative estimation)
% (3) dtImp:		user prescribed value for implicit analyses
dt = [dtf dte TimeIntegration.dtImp];

% Time interval for the transient simulation
nSTEP = 1e3;
if t_end == 0
    tsim = nSTEP*min(dt);
	% Display a warning message in the command window
    warning(['The simulation time is not specified. Therefore, ' ...
            num2str(nSTEP) ' time steps are computed by default.'])
else
    tsim = t_end;
end

% Spectral radius of the amplification matrix
rho_inf = TimeIntegration.rho;  % Applicable to method 403

%% Transient solver
% Start the measurement of the computational time
timerVal = tic;

%% Time integration
if ~flag_mRHS
%------------------Explicit methods - 2nd order accuracy------------------%
    if TimeIntegration.Method == 101 || TimeIntegration.Method == 601
        % Central difference method
        [T,U] = CentralDifferenceMethod(K,M,R,F,F_dt,U0,V0,num_DoF,...
                dof_poi,ActiveDoF,ExciteParam,SignalType,dt,tsim,true);
    elseif TimeIntegration.Method == 102
        % Noh-Bathe method
        [T,U] = NohBatheMethod(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 103
        % Soares method (explicit)
        [T,U] = SoaresMethodExplicit(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 104
        % Mirbagheri method
        [T,U] = MirbagheriMethod(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 105
        % Kim method
        [T,U] = KimMethodExplicit(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 106
        % Zhang-Xing method 1
        [T,U] = ZhangXingMethod1(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 107
        % Zhang-Xing method 2
        [T,U] = ZhangXingMethod2(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 108
        % Explicit generalized alpha method
        [T,U] = ExplicitGeneralizedAlphaMethod(K,M,R,F,F_dt,U0,V0,...
                num_DoF,dof_poi,ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 109
        % Adaptive hybrid explicit method - Soares
        % The method works element by element: the element matrices are
        % saved by call_Dynamic.m (ElemData.mat)
        fileElem = ['.' filesep 'subroutines' filesep 'Temp' filesep ...
                    'ElemData.mat'];
        if ~isfile(fileElem)
            error('Solver_Transient:method109', ...
                 ['Method 109 needs the element matrices (%s), which ' ...
                  'call_Dynamic.m saves for this method.'], fileElem);
        end
        load(fileElem,'ElemData');
        [T,U] = AdaptiveHybridExplicitMethod(K,M,R,F,F_dt,U0,V0,num_DoF,...
                dof_poi,ActiveDoF,ExciteParam,SignalType,dt,tsim,ElemData);
%----------------Explicit methods - higher order accuracy-----------------%
    elseif TimeIntegration.Method == 201
        % Kim method (3rd or 4th order)
        [T,U] = KimExplicit_ho(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim,...
                TimeIntegration.Order);
    elseif TimeIntegration.Method == 202
        % Kim method (4th order)
        [T,U] = KimMethod2Stage(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
%------------------Implicit methods - 2nd order accuracy------------------%
    elseif TimeIntegration.Method == 301
        % Newmark method
        [T,U] = NewmarkMethod(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 302
        % HHT-alpha method
        [T,U] = HHTAlphaMethod(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 303
        % Generalized-alpha method
        [T,U] = GeneralizedAlphaMethod(K,M,R,F,F_dt,U0,V0,num_DoF,...
                dof_poi,ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 304
        % Bathe-Baig method
        [T,U] = BatheBaigMethod(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 305
        % beta1-beta2-Bathe method
        [T,U] = b1b2_BatheMethod(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 306
        % rho-Bathe method
        [T,U] = rho_BatheMethod(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 307
        % Soares method (implicit)
        [T,U] = SoaresMethodImplicit(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 308
        % Kim method (implicit)
        [T,U] = KimMethodImplicit(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 309
        % Kim method (implicit)
        [T,U] = KimMethodImplicitCollocation(K,M,R,F,F_dt,U0,V0,num_DoF,...
                dof_poi,ActiveDoF,ExciteParam,SignalType,dt,tsim);
%----------------Implicit methods - higher order accuracy-----------------%
    elseif TimeIntegration.Method == 401
        % Kim method - higher order accurate time integration scheme (based
        % on Lagrange polynomials)
        [T,U] = KimLagrange_ho(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim,...
                TimeIntegration.Order);
    elseif TimeIntegration.Method == 402
        % Soares method: order depends on the presence of numerical and or
        % physical damping. No damping: order = 4; Only numerical damping:
        % order = 3; Physical damping: order = 2.
        [T,U] = SoaresMethodImplicit_ho(K,M,R,F,F_dt,U0,V0,num_DoF,...
                dof_poi,ActiveDoF,ExciteParam,SignalType,dt,tsim);
    elseif TimeIntegration.Method == 403
        % Interpolate the results
        nInterp = TimeIntegration.nInterp;
        % Pade-based high-order time integration scheme
        if numel(dt) == 3
            dtP = dt(3);
        else
            dtP = dt(1);
        end
        % Number of time steps
        ns = ceil(round(tsim/dtP,3,'decimal'));
        % Adjusted time increment
        dtP = tsim/ns;
        % Order of the Pade expansion (11, 22, 33, 44, 55)
        if TimeIntegration.Order == 2
            order = 11;
        elseif TimeIntegration.Order == 4
            order = 22;
        elseif TimeIntegration.Order == 6
            order = 33;
        elseif TimeIntegration.Order == 8
            order = 44;
        elseif TimeIntegration.Order == 10
            order = 55;
        else
            order = 11;
        end
        % Pre-compute the force vector for a polynomial variation of the
        % force within a time step
        forceHistory = ForcePadeTimeIntegration(order,SignalType,...
                       ExciteParam,ns,dtP);
        % Matrices saved as triangular matrices are converted to full
        % matrices (the unsymmetric matrices of this code are full)
        if (istriu(K) || istril(K)) && ~isscalar(K)
            K0 = K + K';
            M0 = M + M';
            C0 = R + R';
        else
            K0 = K;
            M0 = M;
            C0 = R;
        end
        % Force vector
        if SignalType == 0 && sum(F_dt.Neumann) == 0
            F0 = F;
        elseif SignalType == 0 && sum(F_dt.Neumann) ~= 0
            F0 = F + F_dt.Neumann;
        else
            F0 = F_dt.Neumann;
        end
        % High-order Pade-based implicit time integrator
        % Homogeneous Dirichlet BCs have been prescribed
    	if rho_inf == 1
            % No amplitude decay
          	[T,U] = PadeTimeIntegrator(order,ns,dtP,forceHistory,K0,M0,...
                    C0,F0,U0,V0,[],dof_poi,nInterp);
        else
            % Numerical damping
            order = order - 10;
            if exist('PadeTimeIntegrator_ND_FORTRAN','file')
                % Fortran-solver (not part of the repository)
                [T,U] = PadeTimeIntegrator_ND_FORTRAN(order,rho_inf,ns,...
                        dtP,forceHistory,K0,M0,C0,F0,U0,V0,[],dof_poi,...
                        nInterp);
            else
                % Matlab-solver
                [T,U] = PadeTimeIntegrator_ND(order,rho_inf,ns,dtP,...
                        forceHistory,K0,M0,C0,F0,U0,V0,[],dof_poi,nInterp);
            end
    	end
%-------------Highly precise time integration method (HPTIM)--------------%
    elseif TimeIntegration.Method == 501
        % HPTIM
        [T,U] = HPTIM_Newmark(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim,...
                TimeIntegration.Order);
    else
        % Newmark method --> Default method if a wrong value was chosen by
        % the user!
        [T,U] = NewmarkMethod(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    end
else
    % Deal with multiple load cases simultaneously
    if TimeIntegration.Method < 300
        % Central difference method
        [T,U] = CentralDifferenceMethod_mRHS(K,M,R,F,F_dt,U0,V0,num_DoF,...
                dof_poi,ActiveDoF,ExciteParam,SignalType,dt,tsim,true);
    else
        % HHT-alpha method
        [T,U] = HHTAlphaMethod_mRHS(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,...
                ActiveDoF,ExciteParam,SignalType,dt,tsim);
    end
end
%-------------------------------------------------------------------------%
%%

% End the measurement of the computational time
t_cpu = toc(timerVal);  
disp(['Computational time to solve the equations of motion t_{cpu} = ' ...
      num2str(t_cpu)]);

% Computational time of the FORTRAN solver  
if exist('tF','var')
	t_cpu = tF;
end

end