function [dsp,vel,acc,UAll] = PadeSolverExact(ft,K,M,C,F,u0,v0,pDOF,nnzIC,nUAll,flag_wb)
%PADESOLVEREXACT Time stepping with the exact matrix exponential (reference solution).
%
% Exact for a piecewise linear excitation

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
	x = zeros(2*size(K,1),1);
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
	x = zeros(2*size(K,1),1);
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
