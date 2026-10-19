function x = GaussianImpulsesExcitation(t,a,b,c,d,n)
%-------------------------------------------------------------------------%
% a: Amplitude
% b: Position of the centre of the peak for the first impulse
% c: Standard deviation (controls the width of the impulse) --> only 1% of
%    the amplitude for t = b +-3c
% d: Distance between peaks
% n: Number of impulses
%
% Area below a Guassian function: A = a*c*sqrt(2*pi)
%-------------------------------------------------------------------------%
% Number of time steps
NT = numel(t);

% Parameter
b_vec = zeros(1,n);
for i = 1:n
    b_vec(i) = b + d*(i-1);
end

% Piece-wise function
idx = cell(n+1,1);
for i = 1:n
    if i == 1
        idx{i} = t <= b_vec(i) + d/2;
    else
        idx{i} = t > b_vec(i-1) + d/2 & t <= b_vec(i) + d/2;
    end
end
idx{n+1} = t > b_vec(n) + d/2;

% Save anonymous function in cell array
f = cell(n,1);
for i = 1:n
    bf = b_vec(i);
    f{i} = @(t) a*exp(-(t-bf).^2/(2*c^2));
end
f{n+1} = @(t) zeros(1,numel(t));

% Result vector
x = zeros(NT,1);
counterS = 0;
for i = 1:n+1
    counterE = counterS + sum(idx{i});
    x(counterS+1:counterE) = f{i}(t(idx{i}));
    counterS = counterE;
end

end