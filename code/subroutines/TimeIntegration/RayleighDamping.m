%% Rayleigh Damping
% 22/12/2020
% University of New South Wales (UNSW Sydney)
% Sascha Eisentraeger
%
% Song, Z. & Su, C.: Computation of Rayleigh Damping Coefficients for the
% Seismic Analysis of a Hydro-Powerhouse, Shock and Vibration, Hindawi
% Limited, 2017, 2017, 1-11
%
% The mass and stiffness of the upper and lower structures of a powerhouse
% are different. As such, the first two vibration modes mostly indicate the
% dynamic characteristics of the upper structure, and the precise seismic
% response of a powerhouse is difficult to obtain on the basis of Rayleigh
% damping coefficients acquired using the fundamental frequencies of this
% structure. The damping ratio of each mode is relatively accurate when the
% least square method is used, but the accuracy of the damping ratios that
% contribute substantially to seismic responses is hardly ensured. The
% error of dynamic responses may even be amplified. In this study, modes 
% that greatly influence these responses are found on the basis of mode
% participation mass, and Rayleigh damping coefficients are obtained.
% Seismic response distortion attributed to large differences in Rayleigh
% damping coefficients because of improper modal selection is avoided by
% using the proposed method,which is also simpler and more accurate than
% the least square method. Numerical experiments show that the damping
% matrix determined by using the Rayleigh damping coefficients identified
% by our method is closer to the actual value and the seismic response of
% the powerhouse is more reasonable than that revealed through the least
% square method.
%-------------------------------------------------------------------------%
% Input
%
% omega:    Vector of natural angular frequencies; typically two modes are
%           selected to determine alpha and beta; if several modes are
%           selected the least squares method is used to determine the
%           damping parameter
% zeta:     Vector of damping ratios; two different damping ratios can be
%           provided for 2 modes; if several modes are provided the first
%           value is set as damping ratio
% flag_vis: Switch visualization on or off
%
% The choice of the modes 1 and 2 is in principle arbitrary and depends on
% the a priori knowledge about the structural system. It is a good idea to
% choose the modes of vibration such that all modes that contribute
% significantly are in the range of the chosen two modes. Thus, it is at
% least ensured that the damping ratios are reasonable.
%
%-------------------------------------------------------------------------%
% Output
% alpha:    mass proportional damping parameter [1/s]
% beta:     stiffness proportional damping parameter [s]
%
%-------------------------------------------------------------------------%
function [alpha,beta] = RayleighDamping(omega,zeta,flag_vis)

% If omega is zero no damping is prescribed
if nnz(omega) == 0
    alpha = 0;
    beta  = 0;
    return
end

% Number of modes
nModes = numel(omega);
if nModes == 2
    omega1 = omega(1);
    omega2 = omega(2);
    method = 1;
else
    method = 2;
end

% Damping ratio
nDamp = numel(zeta);
if nDamp == 2
    zeta1 = zeta(1);
    zeta2 = zeta(2);
else
    zeta1 = zeta(1);
    zeta2 = zeta(1);
end

if method == 1
    % Check if omega2 > omega1
    if omega1 > omega2
        % Frequencies
        omega_aux = omega1;
        omega1 = omega2;
        omega2 = omega_aux;
        % Damping ratios
        zeta_aux = zeta1;
        zeta1 = zeta2;
        zeta2 = zeta_aux;
    end

    % Calculate the damping parameters
    if zeta1 == zeta2
        D = 2*zeta1/(omega1+omega2)*[omega1*omega2;1];
    else
        D = 2*omega1*omega2/(omega2^2 - omega1^2)*...
            [omega2,-omega1;-1/omega2,1/omega1]*[zeta1;zeta2];
    end
else
    % Weighted least squares method (weight coefficient exp(-omega_i))
    options = optimset('FunValCheck','on','Diagnostics','off', ...
    'Display','off','TolX',1e-15,'TolFun',1e-15);
    zeta = zeta1;
%     weightingFunc = @(t) exp(-t);
%     weightingFunc = @(t) 1/t;
    weightingFunc = @(t) 1;
    fun = @(x) sum(arrayfun(@(t) weightingFunc(t)*(x(1)/(2*t)+x(2)*t/2-zeta)^2,omega));
    x0 = [0,0];
    D = fminsearch(fun,x0,options);
end

% Rayleigh parameters
alpha = D(1);
beta  = D(2);

%% Visualization of the damping ratio
%-------------------------------------------------------------------------%
% ODE:  m*ddu +            c*du +       k*u = f
%         ddu + 2*zeta*omega*du + omega^2*u = f/m
%
% Mass proportional damping: c = alpha*m
%       m*ddu +      alpha*m*du +       k*u = f
%         ddu +        alpha*du + omega^2*u = f/m
% Compare the term in front of du
% 2*zeta*omega = alpha --> zeta = alpha/(2*omega)
%
% Stiffness proportional damping: c = beta*k
%       m*ddu +       beta*k*du +       k*u = f
%         ddu + beta*omega^2*du + omega^2*u = f/m
% Compare the term in front of du
% 2*zeta*omega = beta*omega^2 --> zeta = beta*omega/2
%-------------------------------------------------------------------------%
if flag_vis
    figure('WindowState','Maximized','Name','Damping ratio')
    hold on; grid on;
    omega_max = max(omega);
    omega_min = min(omega);
    om = omega_min/2:omega_max/1e3:1.25*omega_max;
    zeta_mass  = alpha./(2*om)*100;
    zeta_stiff = beta*om/2*100;
    zeta_ov    = zeta_mass + zeta_stiff;
    plot(om,zeta_ov,'k','LineWidth',2)
    plot(om,zeta_mass,'--k','LineWidth',1)
    plot(om,zeta_stiff,'-.k','LineWidth',1)
    if method == 1
        plot(omega1,zeta1*100,'rx','MarkerSize',12,'LineWidth',2);
        plot(omega2,zeta2*100,'rx','MarkerSize',12,'LineWidth',2);
        text(omega1,zeta1*102.5,...
            {['$\omega_1=$ ' num2str(omega1) ' rad/s'],...
             ['$\zeta_1=$ ' num2str(zeta1*100) ' \%']},...
             'HorizontalAlignment','left','VerticalAlignment','bottom',...
             'Interpreter','latex','FontSize',12)
        text(omega2,zeta2*102.5,...
            {['$\omega_2=$ ' num2str(omega2) ' rad/s'],...
             ['$\zeta_2=$ ' num2str(zeta2*100) ' \%']},...
             'HorizontalAlignment','right','VerticalAlignment','bottom',...
             'Interpreter','latex','FontSize',12)
    else
        plot(omega,ones(size(omega))*zeta1*100,'rx','MarkerSize',6,...
            'LineWidth',1);
    end
    set(gca,'FontSize',18);
    set(gca,'TickLabelInterpreter','latex')
    xlabel('$\omega_\mathrm{n}$ [rad/s]','Interpreter','latex')
    ylabel('$\zeta$ [\%]','Interpreter','latex')
    legend(["$\zeta$","$\zeta_\mathrm{m}$","$\zeta_\mathrm{k}$"],...
            'Interpreter','latex','Location','northeast')
end

% Round to 6 significant digits (ABAQUS)
alpha = round(alpha,6,'significant');
beta  = round(beta,6,'significant');

end