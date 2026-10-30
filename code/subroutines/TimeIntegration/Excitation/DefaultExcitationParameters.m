function ExciteParam = DefaultExcitationParameters(SignalType,tE)
%DEFAULTEXCITATIONPARAMETERS Example parameters of the excitation signals.
%
%   ExciteParam = DefaultExcitationParameters(SignalType, tE) returns the
%   parameter vector that ExcitationSignal expects for the given signal
%   type. tE is the simulation time used by signals 0, 3 and 5. The
%   values are examples; change them in the settings
%   (settings.ExciteParam) rather than here.

switch SignalType
    case 0
        % No time-dependent signal
        ExciteParam = tE;
    case 1
        % Sine burst modulated with a Hann window
        fex = 5e6;      % excitation frequency [Hz]
        nex = 3;        % number of periods
        ExciteParam = [fex, nex];
    case 2
        % Sine burst windowed by an exponential function
        fex = 5e6;      % excitation frequency [Hz]
        ExciteParam = fex;
    case 3
        % Linear chirp
        fch1 = 1e6;     % minimum frequency [Hz]
        fch2 = 5e6;     % maximum frequency [Hz]
        tch  = tE;      % duration of the chirp [s]
        ExciteParam = [fch1, fch2, tch];
    case 4
        % Harmonic signal a*cos(2*pi*fcos*t) + b*sin(2*pi*fsin*t)
        a    = 1;
        b    = 0;
        fcos = 5e6;     % [Hz]
        fsin = 0;       % [Hz]
        ExciteParam = [a, b, fcos, fsin];
    case 5
        % Piecewise linear signal: amplitude table [time, amplitude]
        TS = [0,0;0.005,1;0.015,1;0.02,0;0.04,0;0.045,1;0.055,1;0.06,0;...
              0.08,0;0.085,1;0.095,1;0.1,0;0.12,0;0.125,1;0.135,1;0.14,0];
        ExciteParam = [TS;tE,0];
    case 6
        % Ricker wavelet
        fex = 200;      % central frequency [Hz]
        ExciteParam = fex;
    case 7
        % Gaussian impulses
        a = 1;          % amplitude
        b = 0.01;       % centre of the first impulse [s]
        c = 0.002;      % standard deviation (1% of the amplitude at b+-3c)
        d = 0.04;       % distance between the impulses [s]
        n = 4;          % number of impulses
        ExciteParam = [a, b, c, d, n];
    otherwise
        error('DefaultExcitationParameters:SignalType', ...
              'Unknown signal type %d (valid: 0 to 7).', SignalType);
end

end
