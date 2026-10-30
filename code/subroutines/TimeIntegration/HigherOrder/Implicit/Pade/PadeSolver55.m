function [dsp,vel,acc,UAll] = PadeSolver55(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb)
%PADESOLVER55 Time stepping with the (5,5) Pade approximation (no numerical damping).
%
% Order (5,5)

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
	x = zeros(2*size(K,1),1);
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
