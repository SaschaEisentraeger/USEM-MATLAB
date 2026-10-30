function [dt,NumTimeStep] = TimeStepSize(dtmax,tsim,stabilityFactor)
%TIMESTEPSIZE Time step of a direct time integration method.
%
%   [dt, NumTimeStep] = TimeStepSize(dtmax, tsim) selects the time step of
%   an implicit method from the candidates dtmax = [dtf, dte, dtImp]
%   (excitation, stability limit of the central difference method,
%   prescribed implicit step): dtImp if the three values differ (a step
%   was prescribed), otherwise dtf.
%
%   [dt, NumTimeStep] = TimeStepSize(dtmax, tsim, stabilityFactor) selects
%   the time step of an explicit method: the smallest of dtf, dtImp and
%   stabilityFactor*dte, where stabilityFactor is the ratio of the
%   stability limit of the method to that of the central difference
%   method (1 for the central difference method itself).
%
%   In both cases the step is reduced so that tsim is an integer multiple
%   NumTimeStep of it.

if nargin < 3 || isempty(stabilityFactor)
    % Implicit method
    if numel(unique(dtmax)) == 3
        dt = dtmax(3);
    else
        dt = dtmax(1);
    end
else
    % Explicit method
    dtmax(2) = stabilityFactor*dtmax(2);
    dt = min(dtmax);
end

% Rounding to three decimals avoids an extra step if tsim/dt is an
% integer up to round-off
NumTimeStep = ceil(round(tsim/dt,3,'decimal'));
dt = tsim/NumTimeStep;

fprintf('Time step size used: %6.3e\n\n',dt);

end
