function x = PiecewiseLinearExcitationDerivative(TS,t,order)
% Number of time steps
NT = numel(t);

% Table for the values of the time series [time,amplitude]
time = TS(:,1);
amplitude = TS(:,2);

% Number of intervals
NI = numel(amplitude)-1;

% Identifier for the different times
ID = cell(NI,1);
for i = 1:NI
    if i ~= NI
        ID{i} = t <  TS(i+1) & t >= TS(i);
    else
        ID{i} = t <= TS(i+1) & t >= TS(i);
    end
end

% Generate the linear functions
Slope = zeros(NI,1);
for i = 1:NI
    % Slope of the function
    Slope(i) = (amplitude(i)-amplitude(i+1))/(time(i)-time(i+1));
end

% Compute values at t
x = zeros(NT,1);
if order == 1   % higher order derivatives are always zero
    for i = 1:NI
        Val = Slope(i)*ones(NT,1);
        x(ID{i}) = Val(ID{i});
    end
    x(abs(x) < eps) = 0;
end

end