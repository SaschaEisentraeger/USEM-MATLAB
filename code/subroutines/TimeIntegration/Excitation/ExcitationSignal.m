function x = ExcitationSignal(SignalType,ExciteParam,t)

if SignalType == 1
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

if SignalType == 0
	x = zeros(size(t));
elseif SignalType == 1
    % Time-dependent excitation function: Sine burst signal (modulated by a
    % Hann-window)
    x = SineBurstExcitation(t,f_d,n_d);
elseif SignalType == 2
    % Time-dependent excitation function: Burst signal
    x = BurstExcitation(t,f_d);
elseif SignalType == 3
    % Time-dependent excitation function: Linear chirp signal
    x = PeriodicChirpExcitation(t,f0,t1,f1);
elseif SignalType == 4
    % Harmonic excitation
    x = a*cos(2*pi*f0*t) + b*sin(2*pi*f1*t);
elseif SignalType == 5
    % Piecewise linear function excitation
    x = PiecewiseLinearExcitation(TS,t);
elseif SignalType == 6
    % Ricker wavelet
    x = RickerExcitation(t,f_d);
elseif SignalType == 7
    % Gaussian impulses
    x = GaussianImpulsesExcitation(t,a,b,c,d,n);
end

% x should be a column vector
if iscolumn(x)
    x = x';
end

end