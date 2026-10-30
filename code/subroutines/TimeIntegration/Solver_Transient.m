function [T,U,t_cpu] = Solver_Transient(S)
%SOLVER_TRANSIENT Integrate the semi-discrete equations of motion in time.
%
%   [T, U, t_cpu] = Solver_Transient(S) integrates
%       M*a + R*v + K*u = F + F_dt.Neumann*s(t) + F_dt.K*s(t)
%                         + F_dt.C*ds/dt + F_dt.M*d2s/dt2
%   with the method S.TimeIntegration.Method (TimeIntegrationMethods) and
%   returns the time instants T, the displacement history U and the CPU
%   time of the time integration.
%
%   S is a struct (built by DynamicAnalysis) or the name of a MAT file that
%   contains its fields (Solver_input.mat written by DynamicAnalysis):
%     K, M, R, F       system matrices and constant load (all DOFs or the
%                      active DOFs only)
%     F_dt             struct with fields Neumann, K, C, M (see above)
%     U0, V0           initial displacements and velocities
%     ActiveDoF        DOFs without Dirichlet conditions
%     NOI              nodes of interest (output if save_var = 0)
%     num_DoF          number of all DOFs
%     save_var         0: DOFs of the nodes of interest, 1: all DOFs
%     t_end            simulation time; 0: 1000 time steps
%     dtf, dte         time step candidates of the excitation and of the
%                      stability limit (see DynamicAnalysis)
%     SignalType, ExciteParam
%                      excitation signal s(t) (ExcitationSignal)
%     TimeIntegration  Method, Order, rho, dtImp, nInterp
%     flag_mRHS        one column of F per load case: the central
%                      difference method (Method < 300) or the HHT-alpha
%                      method (Method >= 300) is used for every Method
%     ElemData         element matrices (method 109 only)
%
%   Reference: K. J. Bathe, Finite-Elemente-Methoden, Chap. 9, Springer,
%   2002.

if ischar(S) || isstring(S)
    S = load(S);
end

K = S.K;
M = S.M;
R = S.R;
F = S.F;
F_dt = S.F_dt;
U0 = S.U0;
V0 = S.V0;
ActiveDoF = S.ActiveDoF;
num_DoF = S.num_DoF;
TimeIntegration = S.TimeIntegration;

% Remove the DOFs with Dirichlet conditions
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

dof_poi = OutputDoFs(S.NOI,S.save_var,ActiveDoF,num_DoF);

% Time step candidates (each method chooses from them):
% dtf    based on the maximum frequency of the excitation signal
% dte    stability limit of the central difference method from the
%        highest eigenfrequency of the elements (conservative estimate)
% dtImp  prescribed time step of the implicit methods
disp('Determine the time step size');
dt = [S.dtf S.dte TimeIntegration.dtImp];

% Simulation time
nSTEP = 1e3;
if S.t_end == 0
    tsim = nSTEP*min(dt);
    warning(['The simulation time is not specified. Therefore, ' ...
             num2str(nSTEP) ' time steps are computed by default.'])
else
    tsim = S.t_end;
end

% Method and its trailing arguments
if S.flag_mRHS
    % Multiple load cases (columns of F) at once
    if TimeIntegration.Method < 300
        method = @CentralDifferenceMethod_mRHS;
        extra = {true};
    else
        method = @HHTAlphaMethod_mRHS;
        extra = {};
    end
else
    opt = TimeIntegration;
    if TimeIntegration.Method == 109
        % Adaptive hybrid explicit method: works element by element
        if ~isfield(S,'ElemData')
            error('Solver_Transient:method109', ...
                 ['Method 109 needs the element matrices (field ' ...
                  'ElemData), which DynamicAnalysis provides for this ' ...
                  'method.']);
        end
        opt.ElemData = S.ElemData;
    end
    [method,extra] = SelectTimeIntegrator(TimeIntegration.Method,opt);
end

timerVal = tic;
[T,U] = method(K,M,R,F,F_dt,U0,V0,num_DoF,dof_poi,ActiveDoF,...
        S.ExciteParam,S.SignalType,dt,tsim,extra{:});
t_cpu = toc(timerVal);
disp(['Computational time to solve the equations of motion t_{cpu} = ' ...
      num2str(t_cpu)]);

end
