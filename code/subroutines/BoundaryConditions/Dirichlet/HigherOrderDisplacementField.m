%% Polynomial Displacement Field - Plane Stress State
clearvars; close all; clc;
% Symbolic variables
syms x y E nu
% Polynomial order (p >= 2; the equilibrium conditions are fulfilled for
% any linear displacement field independent of the choice of the
% coefficients as second derivatives of the displacement functions are
% required)
p = 8;
% Number of terms in a complete polynomial of order p
n = (p+1)*((p+1)+1)/2;
% Vector of monomials
P = sym('P',[n,1]);
inc = 0;
for i = 0:p
    px = i;
   	py = 0;
    for j = 0:i
        if px+py <= p
            inc = inc + 1;
            P(inc) = x^px*y^py;
            px = px - 1;
            py = py + 1;
        end
    end
end

% Displacement in x-direction: ux
% Unknown coefficients/symbolic variables
A = sym(zeros(1,n));
for i = 1:n
    syms(['a' num2str(i)]);
    A(i) = eval(['a' num2str(i)]);
end
ux = A*P;

% Displacement in y-direction: uy
% Unknown coefficients
B = sym(zeros(1,n));
for i = 1:n
    syms(['b' num2str(i)]);
    B(i) = eval(['b' num2str(i)]);
end
uy = B*P;

% Strains
epsilon_x = diff(ux,x);
epsilon_y = diff(uy,y);
gamma_xy  = diff(ux,y) + diff(uy,x);

% Stresses
sigma_x = E/(1-nu^2)*(epsilon_x + nu*epsilon_y);
sigma_y = E/(1-nu^2)*(epsilon_y + nu*epsilon_x);
tau_xy  = E/(2+2*nu)*gamma_xy;

% Equilibrium equations
Eq1 = diff(sigma_x,x) + diff(tau_xy,y);
Eq2 = diff(tau_xy,x)  + diff(sigma_y,y);
% Eq1 and Eq2 have to zero to fulfill the equilibrium equations. Therefore,
% we can choose a certain number of constants to fulfill the equilibrium
% conditions. Substitute values for the coefficient vectors a and b.

% Do not prescribe coefficients that are encountered twice in these terms
% (i.e. all mixed terms are not prescribed and need to be solved for or in
% other words each term that involves a product of x^i*y^j with i,j >= 1,
% should be left in the system of equations):
% Pascal's triangle:
%                                   1
%                            x              y
%                     x^2           xy          y^2
%               x^3         x^2y          xy^2         y^3
%          x^4       x^3y         x^2y^2        xy^3         y^4
%-------------------------------------------------------------------------%
%                                   1
%                           2               3
%                      4            5            6
%                7          8               9         10
%           11         12          13            14         15
%-------------------------------------------------------------------------%

% Number of coefficients that need to be prescribed
n_pre = (p-1)*((p-1)+1)/2;
% Vector of indices to mixed terms; these coefficients are not prescribed
% by an arbitrary value!
idx = zeros(1,n_pre);
% Start and end counter
count_s = 0;
count_e = 0;
% Start and end index for the coefficients
idx_s = 5;
idx_e = 5;
for i = 1:p-1
    % Update the counters
    count_s = count_e + 1;
    count_e = count_s + i - 1;
    % Indices to the coefficients
    idx(count_s:count_e) = idx_s:idx_e;
    % Update start and and index for the next loop
    idx_s = idx_e + 3;
    idx_e = idx_s + i;
end

for i = 1:n
    % Separate the equations into different terms according to the power of
    % x and y (use: CEq1 = coeffs(Eq1,[x,y]); CEq2 = coeffs(Eq2,[x,y]);).
    % The mixed terms are encountered in both Eq1 and Eq2, while all terms
    % x^i are only part of Eq1 and all terms y^i are only contained in Eq2.
	% The mixed terms are coupled through tau_xy in both equations.
    
    % Prescribe 'free' coefficients/ coefficient vectors a and b
	if ~ismember(i,idx)
        eval(['a' num2str(i) '=' num2str(i) ';']);
     	eval(['b' num2str(i) '=' num2str(i) ';']);
%         eval(['a' num2str(i) '=' num2str(1) ';']);
%      	eval(['b' num2str(i) '=' num2str(1) ';']);
	end
end
% Separate the coefficients of the equations
CEq1 = coeffs(Eq1,[x,y]);
CEq2 = coeffs(Eq2,[x,y]);

% Number of equations per equilibrium condition
nEqs1 = size(CEq1,2);
nEqs2 = size(CEq2,2);

% Substitute the prescribed coefficients in the equations that are
% separated in powers of x and y.
for i = 1:nEqs1
    % Write the formula as a string and evaluate the expression
    eval(['Eq' num2str(i) '_subs = subs(CEq1(' num2str(i) '));']);
end
for i = 1:nEqs2
    % Write the formula as a string and evaluate the expression
    eval(['Eq' num2str(nEqs1+i) '_subs = subs(CEq2(' num2str(i) '));']);
end

% Set up the system of equations
EEcell = cell(1,nEqs1+nEqs2);
if p > 1
    for i = 1:nEqs1+nEqs2
        if i == 1
            EEcell{i} = ['[Eq' num2str(i) '_subs == 0,'];
        elseif i == nEqs1+nEqs2
            EEcell{i} = ['Eq' num2str(i) '_subs == 0];'];
        else
            EEcell{i} = ['Eq' num2str(i) '_subs == 0,'];
        end
    end
    EEstr = strcat(EEcell{:});
    EE = eval(EEstr);

    % Solve the system of equations
    Scell = cell(1,nEqs1+nEqs2);
    for i = 1:nEqs1
        if i == 1
            Scell{i} = ['[a' num2str(idx(i)) ','];
        else
            Scell{i} = ['a' num2str(idx(i)) ','];
        end
    end
    for i = nEqs1+1:nEqs1+nEqs2
        if i == nEqs1+nEqs2
            Scell{i} = ['b' num2str(idx(i-nEqs1)) ']'];
        else
            Scell{i} = ['b' num2str(idx(i-nEqs1)) ','];
        end
    end
    Sstr = strcat(Scell{:});
    S = solve(EE,eval(Sstr));

    % Poisson's ratio
    nu = 0.3;
    for i = 1:nEqs1
        str = ['a' num2str(idx(i)) ' = eval(subs(S.a' num2str(idx(i)) '));'];
        eval(str);
    end
    for i = 1:nEqs2
        str = ['b' num2str(idx(i)) ' = eval(subs(S.b' num2str(idx(i)) '));'];
        eval(str);
    end
end

%% Test the accuracy
% [DispVal_x] = SpatialDisplacementFunction([1 1],[xCoor yCoor]);
% err_x = (DispVal_x - ux)./DispVal_x;
% idx_x = find(abs(err_x) > 1e-9);
% err_mean_x= mean(err_x);

% [DispVal_y] = SpatialDisplacementFunction([1 1],[xCoor yCoor]);
% err_y = (DispVal_y - ux)./DispVal_y;
% idx_y = find(abs(err_y) > 1e-9);
% err_mean_y= mean(err_y);