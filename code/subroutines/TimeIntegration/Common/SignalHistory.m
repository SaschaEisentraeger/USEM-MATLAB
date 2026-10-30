function [Amp,dAmp,ddAmp] = SignalHistory(SignalType,ExciteParam,t)
%SIGNALHISTORY Excitation signal and its first two time derivatives.
%
%   [Amp, dAmp, ddAmp] = SignalHistory(SignalType, ExciteParam, t) returns
%   s(t), ds/dt and d2s/dt2 of the excitation signal at the times t (see
%   ExcitationSignal and ExcitationSignalDerivative). The load at time t
%   is then LoadAtTime(F, Fdt, Amp(i), dAmp(i), ddAmp(i)).

Amp   = ExcitationSignal(SignalType,ExciteParam,t);
dAmp  = ExcitationSignalDerivative(SignalType,ExciteParam,t,1);
ddAmp = ExcitationSignalDerivative(SignalType,ExciteParam,t,2);

end
