function x = ExcitationSignalDerivative(SignalType,ExciteParam,t,order)

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
    x = SineBurstExcitationDerivative(t,f_d,n_d,order);
elseif SignalType == 2
    % Time-dependent excitation function: Burst signal
    x = BurstExcitationDerivative(t,f_d,order);
elseif SignalType == 3
    % Time-dependent excitation function: Linear chirp signal
    x = PeriodicChirpExcitationDerivative(t,f0,t1,f1,order);
elseif SignalType == 4
    % Harmonic excitation
    if order == 1
        x =  -2*a*f0*pi*sin(2*pi*f0*t)     +  2*b*f1*pi*cos(2*pi*f1*t);
    elseif order == 2
        x =  -4*a*f0^2*pi^2*cos(2*pi*f0*t) -  4*b*f1^2*pi^2*sin(2*pi*f1*t);
    elseif order == 3
        x =   8*a*f0^3*pi^3*sin(2*pi*f0*t) -  8*b*f1^3*pi^3*cos(2*pi*f1*t);
    elseif order == 4
        x =  16*a*f0^4*pi^4*cos(2*pi*f0*t) + 16*b*f1^4*pi^4*sin(2*pi*f1*t);
   elseif order == 5
        x = -32*a*f0^5*pi^5*sin(2*pi*f0*t) + 32*b*f1^5*pi^5*cos(2*pi*f1*t);
    end
elseif SignalType == 5
    % Piecewise linear function excitation
    x = PiecewiseLinearExcitationDerivative(TS,t,order);
elseif SignalType == 6
    % Ricker wavelet
    x = RickerExcitationDerivative(t,f_d,order);
elseif SignalType == 7
    % Gaussian impulses
    x = GaussianImpulsesExcitationDerivative(t,a,b,c,d,n,order);
end

% x should be a column vector
if iscolumn(x)
    x = x';
end

end