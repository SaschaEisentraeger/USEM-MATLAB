function [dsp,vel,acc,UAll] = PadeSolverND45(rho_infinity,ft,K,M,C,F,u0,...
                              v0,pDOF,nnzIC,nUAll,flag_wb)
%PADESOLVERND45 Time stepping with the damped Pade scheme of order (4,5).
%
% Order (4,5)

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

[pcoe, qcoe, rs, cf] = TimeIntgCoeff(5,rho_infinity);

r1 = real(rs(1));
Kd1 = sparse((r1*r1)*M + r1*C + K);
dKd1 = decomposition(Kd1);

r2 = rs(2);
Kd2 = sparse((r2*r2)*M + r2*C + K);
[L2,U2,p2,q2] = lu(Kd2,'vector');

r3 = rs(3);
Kd3 = sparse((r3*r3)*M + r3*C + K);
[L3,U3,p3,q3] = lu(Kd3,'vector');

ndof = size(K,1);

f  = [Fb; zeros(ndof,1)];
Af1 = SolverPadeAx(dM,K,C,f);
Af2 = SolverPadeAx(dM,K,C,Af1);
Af3 = SolverPadeAx(dM,K,C,Af2);
Af4 = SolverPadeAx(dM,K,C,Af3);
Pf = [f Af1 Af2 Af3 Af4]*cf';

npt = size(ft,2);
Pf = Pf(:,1:npt);
pt = [1 0.5.^(1:npt-1)]';

if nnzIC ~= 0
	x = [v0; u0];
    Ax1 = SolverPadeAx(dM,K,C,x);
    Ax2 = SolverPadeAx(dM,K,C,Ax1);
    Ax3 = SolverPadeAx(dM,K,C,Ax2);
	Ax4 = SolverPadeAx(dM,K,C,Ax3);
	Ax5 = SolverPadeAx(dM,K,C,Ax4);
	
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
	Ax5 = x;
end

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
    
    b = [x Ax1 Ax2 Ax3 Ax4 Ax5]*pcoe' + Pf*ft(is,:)';
    
	tmp = dKd1\(M*(r1*b(1:ndof)) - K*b(ndof+1:2*ndof));
    
    x1 = [tmp ; (tmp+b(ndof+1:2*ndof))/r1];
	
	ctmp = r2*(M*x1(1:ndof)) - K*x1(ndof+1:2*ndof);
    
    tmp = U2\(L2\ctmp(q2,:));
    tmp(p2) = tmp;
    
    y2 = [tmp ; (tmp+x1(ndof+1:2*ndof))/r2];
    x2 =  -imag(y2)/imag(r2);
	
	ctmp = r3*(M*x2(1:ndof)) - K*x2(ndof+1:2*ndof);
    
    tmp = U3\(L3\ctmp(q3,:));
    tmp(p3) = tmp;
    
    y3 = [tmp ; (tmp+x2(ndof+1:2*ndof))/r3];
    xb =  -imag(y3)/imag(r3);
	
	x = xb;
    
	Ax1b = -imag(r3*y3)/imag(r3);
	Ax1 = Ax1b;
	Ax2b = x2 - real(r3*conj(r3))*xb + 2*real(r3)*Ax1b;
	Ax2 = Ax2b;
	Ax3b = real(conj(r2)*x2 - y2 - (r3*conj(r3)*Ax1b-2*real(r3)*Ax2b));
 	Ax3 = Ax3b;
    Ax4b = real(x1 - r2*conj(r2)*x2 + (...
                      2*real(r2)*r3*conj(r3)*Ax1b - ...
           (r3*conj(r3)+4*real(r2)*real(r3))*Ax2b + ...
                     (2*real(r2)+2*real(r3))*Ax3b));
    Ax4 = Ax4b;
    Ax5 = (b - [x Ax1 Ax2 Ax3 Ax4]*qcoe(1:5)')/qcoe(6);
	
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
