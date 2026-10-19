function x = RickerExcitation(t,f_d)
% Number of time steps
NT = numel(t);

% Time parameter
t0 = 1/f_d;

% Amplitude
x = zeros(1,NT);

% Time-dependent excitation function: Ricker wavelet
for i = 1:NT
    a = sqrt(2)*pi;
    Psi = (a*f_d*(t(i)-t0))^2;
	x(i) = (1-Psi)*exp(-1/2*Psi);
end

x(abs(x) < eps) = 0;

end