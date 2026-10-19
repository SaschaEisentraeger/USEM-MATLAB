function x = BurstExcitation(t,f_d)
% Number of time steps
NT = numel(t);

% Time-dependent excitation function: Burst signal
mu = 4/f_d;
sigma = 1/f_d;
x = zeros(1,NT);
for i = 1:NT
	x(i) = sin(2*pi*f_d*t(i)).*exp(-0.5*((t(i)-mu)/sigma).^2);
end

x(abs(x) < eps) = 0;

end