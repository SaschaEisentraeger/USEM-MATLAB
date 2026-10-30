%% Spectral radius of the Pade-based time integrator
% clearvars; clc; close all; format long

% Order of the Pade-expansion
M = 3;
L = M-1;

% Weighting factor [0,1]: diagonal Pade expansion (M,M) and (L,M)
a = [1,3/4,1/2,1/4,0];

% Number of different weighting parameters
na = numel(a);

% Initialize symbolic variable x
syms x
% Exponential function
func = exp(x);

% Pade series
[CoeffN_diag,CoeffD_diag] = PadeSeriesCoeff(func,M,M,false);
[CoeffN,CoeffD] = PadeSeriesCoeff(func,L,M,false);

% Mixed scheme
CoeffN_mixed = cell(1,na);
CoeffD_mixed = cell(1,na);
for k = 1:na
    CoeffN_mixed{k} = a(k)*CoeffN_diag + (1-a(k))*[CoeffN,zeros(1,M-L)];
    CoeffD_mixed{k} = a(k)*CoeffD_diag + (1-a(k))*CoeffD;
end

% Compute the spectral radius
% Argument
T = 0.01:0.01:100;
% Spectral radius
rho = cell(1,na);
for k = 1:na
    rho{k} = abs(polyval(flip(CoeffN_mixed{k}),T*1j)./...
                 polyval(flip(CoeffD_mixed{k}),T*1j));
end

% Visualization
figure('WindowState','maximized','Name','Spectral radius curve')
hold on
grid on
str_legend = strings(1,na);
for k = 1:na
    plot(T,rho{k},'LineWidth',2)
    str_legend(k) = ['Pad\''e (' num2str(M) num2str(M) '--' num2str(L) num2str(M) '): $a=' num2str(a(k)) '$'];
end
legend(str_legend,'Interpreter','latex','Location','southwest')
set(gca,'xscale','log');
ylim([0 1.2])
xlabel('$T$ [-]','Interpreter','latex')
ylabel('$\rho_\infty$ [-]','Interpreter','latex')
set(gca,'FontSize',18);