function [dsp,vel,acc,UAll] = PadeSolver33(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb)
%PADESOLVER33 Time stepping with the (3,3) Pade approximation (no numerical damping).
%
% Order (3,3)

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
	x = zeros(2*size(K,1),1);
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
