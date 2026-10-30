function [dsp,vel,acc,UAll] = PadeSolverND12(rho_infinity,ft,K,M,C,F,u0,...
                              v0,pDOF,nnzIC,nUAll,flag_wb)
%PADESOLVERND12 Time stepping with the damped Pade scheme of order (1,2).
%
% Order (1,2)

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

[pcoe,qcoe,rs,cf] = TimeIntgCoeff(2,rho_infinity);
r = rs(1);
Kd = sparse((r*r)*M + r*C + K);
[L,U,p,q] = lu(Kd,'vector');

ndof = size(K,1);

f  = [Fb; zeros(ndof,1)];
Af1 = SolverPadeAx(dM,K,C,f);
Pf = [f Af1]*cf';

npt = size(ft,2);
Pf = Pf(:,1:npt);
pt = [1 0.5.^(1:npt-1)]';

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

if flag_wb
    % Waitbar
    h = waitbar(0,'3rd-order Pade-based time integrator (implicit):');

    % Less updates of the waitbar (avoid slowing the simulation down)
    nUpdate = ceil(ns/1000);
end

disp('Solve the equations of motion using the Pade-based time integrator');
disp('');
disp('---------------------------------------------------------------------');
disp('');

for it = 2:ns
	is = it - 1;
	
    b = [x Ax1 Ax2]*pcoe' + Pf*ft(is,:)';
    
    ctmp = r*(M*b(1:ndof)) - K*b(ndof+1:2*ndof);
    
    tmp = U\(L\ctmp(q,:));
    tmp(p) = tmp;
    
    y = [tmp ; (tmp+b(ndof+1:2*ndof))/r];
    xb = -imag(y)/imag(r);
    
    x =  xb;

    Ax1b = - imag(r*y)/imag(r);
    Ax1 = Ax1b;
    Ax2 = (b - qcoe(1)*x - qcoe(2)*Ax1)/qcoe(3);
    
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
            waitbar(it/ns,h,['3rd-order Pade-based method (implicit): ' ...
                    num2str(it/ns*100,'%5.1f') '%'])
        end
    end
    
end

if flag_wb
    % Delete the waitbar
    close(h)
end

end
