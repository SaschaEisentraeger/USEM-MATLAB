function [dsp,vel,acc,UAll] = PadeSolver11(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb)
%PADESOLVER11 Time stepping with the (1,1) Pade approximation (no numerical damping).
%
% Order (1,1)

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
	x = zeros(2*size(K,1),1);
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
