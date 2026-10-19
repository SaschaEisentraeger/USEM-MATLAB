function x = PeriodicChirpExcitation(t,f0,t1,f1)
% Number of time steps
NT = numel(t);

% Time-dependent excitation function: Linear chirp signal
x = zeros(1,NT);

% Frequency constant k
k = (f1-f0)/t1;

for i = 1:NT
    if t(i)<t1
        % x(i) = chirp(t(i),f0,t1,f1,'linear');
        x(i) = cos(2*pi*(f0+k*t(i))*t(i));
    else
        x(i) = 0;
    end
end

end