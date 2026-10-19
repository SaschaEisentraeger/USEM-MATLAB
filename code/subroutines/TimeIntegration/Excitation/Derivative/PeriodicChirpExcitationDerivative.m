function x = PeriodicChirpExcitationDerivative(t,f0,t1,f1,order)
% Number of time steps
NT = numel(t);

% Time-dependent excitation function: Linear chirp signal
x = zeros(1,NT);

% Frequency constant k
k = (f1-f0)/t1;

for i = 1:NT
    if order == 1
        if t(i)<t1
            x(i) = -sin(2*pi*t(i)*(f0 + k*t(i)))* ...
                   (2*pi*(f0 + k*t(i)) + 2*pi*k*t(i));
        else
            x(i) = 0;
        end
    elseif order == 2
        if t(i)<t1
            x(i) = -cos(2*pi*t(i)*(f0 + k*t(i)))*         ...
                   (2*pi*(f0 + k*t(i)) + 2*pi*k*t(i))^2 - ...
                   4*k*pi*sin(2*pi*t(i)*(f0 + k*t(i)));
        else
            x(i) = 0;
        end
    elseif order == 3
        if t(i)<t1
            x(i) = sin(2*pi*t(i)*(f0 + k*t(i)))*          ...
                   (2*pi*(f0 + k*t(i)) + 2*pi*k*t(i))^3 - ...
                   12*k*pi*cos(2*pi*t(i)*(f0 + k*t(i)))*  ...
                   (2*pi*(f0 + k*t(i)) + 2*pi*k*t(i));
        else
            x(i) = 0;
        end
    elseif order == 4
    elseif order == 5
    end
end

end