function x = SineBurstExcitation(t,f_d,n_d)
% Number of time steps
NT = numel(t);

% Time-dependent excitation function: Sine burst signal (modulated by a
% Hann-window)
x = zeros(1,NT);
for i = 1:NT
    if t(i)<n_d/f_d && t(i) >= 0
        x(i) = sin(2*pi*f_d*t(i))*sin(pi*f_d*t(i)/n_d)^2;
    else
        x(i) = 0;
    end
end

end