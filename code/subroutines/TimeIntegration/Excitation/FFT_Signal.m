%% Fast Fourier Transform of the excitation signal
% 20.04.2017
% University of Magdeburg
% Sascha Duczek
%-------------------------------------------------------------------------%
% Input
% ExciteParam:  Parameters describing the excitation signal
% SignalType:   Type of the excitation signal
%-------------------------------------------------------------------------%
% Output
% f_max:        Maximum frequency of interest
%-------------------------------------------------------------------------%
function [f_max] = FFT_Signal(ExciteParam,SignalType)

if SignalType == 0
    tE = ExciteParam(1);    	% simulation time
elseif SignalType == 1
    f_d = ExciteParam(1);       % excitation frequency
    n_d = ExciteParam(2);       % number of periods
elseif SignalType == 2
    f_d = ExciteParam(1);       % excitation frequency
elseif SignalType == 3
    f0  = ExciteParam(1);   	% minimum frequency
    f1  = ExciteParam(2);   	% maximum frequency
    t1  = ExciteParam(3);   	% duration of the chirp
elseif SignalType == 4
    a   = ExciteParam(1);       % amplitude of the cosine function
    b   = ExciteParam(2);       % amplitude of the sine function
    f0  = ExciteParam(3);   	% frequency of the cos-term
    f1  = ExciteParam(4);   	% frequency of the sin-term
elseif SignalType == 5
    TS  = ExciteParam;          % amplitude tableau
elseif SignalType == 6
    f_d = ExciteParam(1);       % excitation frequency
elseif SignalType == 7
    a = ExciteParam(1);         % amplitude
    b = ExciteParam(2);         % Position of the center of the peak for
                                % the first impulse in [s]
    c = ExciteParam(3);         % Standard deviation (controls the width of
                                % the % impulse) --> only 1% of the
                                % amplitude for t = b +-3c
    d = ExciteParam(4);         % distance between peaks
    n = ExciteParam(5);         % number of peaks
end

% Plot the time- and frequency-domain signal
flag_vis = false;
% Save the plots
flag_fig_save = 0;

% Sampling frequency - Fs/2 is the highest frequency contained in the fft
if SignalType == 0
    % Constant force
    Fs = 10;
elseif SignalType == 3
    % Linear chirp
    Fs = f1*10;
elseif SignalType == 4
    % Harmonic excitation
    Fs = min([f0,f1])*10;
elseif SignalType == 5
    % Piecewise linear excitation signal
    % Find peaks in the excitation signal
    [pks,loc] = findpeaks(TS(:,2));
    % Time between peaks
    if numel(pks) == 1
        DT = TS(loc,1);
        FEX = 1/(4*DT);
        Fs = FEX*50;
    elseif isempty(pks)
        DT = TS(end,1) - TS(1,1);
        FEX = 1/DT;
        Fs = FEX*1000;
    else
        DT = TS(loc(2),1) - TS(loc(1),1);
        FEX = 1/(2*DT);
        Fs = FEX*75;
    end
elseif SignalType == 7
    % Gaussian impulses
    f_d = 1/c;
    Fs = f_d*50;
else
    Fs = f_d*50;
end

% Sampling period (time increment dt)
T = 1/Fs;

% Length of signal (number of increments)
L = 1e5;

% Time vector
t = (0:L-1)*T;

if SignalType == 0
    % Constant force
    S = ones(size(t));
    S(t > tE) = 0;
elseif SignalType == 1
    % Time-dependent excitation function: Sine burst signal (modulated by a
    % Hann-window)
    S = SineBurstExcitation(t,f_d,n_d);
elseif SignalType == 2
    % Time-dependent excitation function: Burst signal
    S = BurstExcitation(t,f_d);
elseif SignalType == 3
    % Time-dependent excitation function: Linear chirp signal
    S = PeriodicChirpExcitation(t,f0,t1,f1);
elseif SignalType == 4
    % Harmonic excitation
    S = a*cos(2*pi*f0*t) + b*sin(2*pi*f1*t);
elseif SignalType == 5
    % Piecewise linear function excitation
    S = PiecewiseLinearExcitation(TS,t);
elseif SignalType == 6
    % Time-dependent excitation function: Burst signal
    S = RickerExcitation(t,f_d);
elseif SignalType == 7
    % Gaussian impulses
    S = GaussianImpulsesExcitation(t,a,b,c,d,n);
end

%-------------------------------------------------------------------------%

% Fourier transform of the signal
Y = fft(S);

% Two- and single sided spectra
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Define the frequency domain (frequency increment df = Fs/L)
f = Fs*(0:(L/2))/L;

% Determine the frequency range of interest
% Highest amplitude
[Amp_max,idx_amp] = max(P1);
% Find the maximum frequency of interest
idx_max = find(P1(idx_amp:end) > 1e-2*Amp_max,1,'last');
f_max = f(idx_amp + idx_max - 1);

%% Visualization
if flag_vis
    figure('WindowState','maximized','Name','Signal - Time domain')
    hold on; grid on; box on
    plot(t,S,'k','LineWidth',2)
    xlabel('$t$ [s]','Interpreter','latex')
    ylabel('$p(t)$','Interpreter','latex')
    axis square
    set(gca,'FontSize',18);
    set(gca,'TickLabelInterpreter','latex')
    hold off
%     xlim([0 0.25])
    if flag_fig_save == 1
        set(gca,'Units','centimeters','Position',[3,3,6.9,6.0],'FontSize',9)
        iptsetpref('ImshowBorder','tight')
        print(gcf,'-depsc','-painters','p(t).eps')
        savefig('p(t).fig')
    end

    figure('WindowState','maximized','Name','Signal - Frequency domain')
    hold on; grid on; box on
    plot(f,P1,'k','LineWidth',2)
    xlabel('$f$ [Hz]','Interpreter','latex')
    ylabel('$\hat{P}(f)$','Interpreter','latex')
    axis square
    set(gca,'FontSize',18);
    set(gca,'TickLabelInterpreter','latex')
    hold off
%     xlim([0 100])
    if flag_fig_save == 1
        set(gca,'Units','centimeters','Position',[3,3,6.9,6.0],'FontSize',9)
        iptsetpref('ImshowBorder','tight')
        print(gcf,'-depsc','-painters','P(f).eps')
        savefig('P(f).fig')
    end
end

end