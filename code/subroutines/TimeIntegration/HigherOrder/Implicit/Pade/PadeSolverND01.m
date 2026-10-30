function [dsp,vel,acc,UAll] = PadeSolverND01(rho_infinity,ft,K,M,C,F,u0,...
                              v0,pDOF,nnzIC,nUAll,flag_wb)
%PADESOLVERND01 Time stepping with the damped Pade scheme of order (0,1).
%
% Order (0,1)

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

[pcoe,qcoe,rs,cf] = TimeIntgCoeff(1,rho_infinity);
r = rs(1);
Kd = sparse((r*r)*M + r*C + K);
dKd = decomposition(Kd);

ndof = size(K,1);

f  = [Fb; zeros(ndof,1)];
Pf = f*cf';

npt = size(ft,2);
Pf = Pf(:,1:npt);
pt = [1 0.5.^(1:npt-1)]';

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

if flag_wb
    % Waitbar
    h = waitbar(0,'1st-order Pade-based time integrator (implicit):');

    % Less updates of the waitbar (avoid slowing the simulation down)
    nUpdate = ceil(ns/1000);
end

disp('Solve the equations of motion using the Pade-based time integrator');
disp('');
disp('---------------------------------------------------------------------');
disp('');

for it = 2:ns
	is = it - 1;
	
    b = [x Ax1]*pcoe' + Pf*ft(is,:)';
    
    ctmp = r*(M*b(1:ndof)) - K*b(ndof+1:2*ndof);
    
    tmp = dKd\ctmp;
    xb = [tmp ; (tmp+b(ndof+1:2*ndof))/r];
    
    x = xb;

    Ax1 = (b - qcoe(1)*x)/qcoe(2);
    
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
            waitbar(it/ns,h,['1st-order Pade-based method (implicit): ' ...
                    num2str(it/ns*100,'%5.1f') '%'])
        end
    end
    
end

if flag_wb
    % Delete the waitbar
    close(h)
end

end
