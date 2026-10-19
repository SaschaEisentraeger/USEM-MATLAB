function x = GaussianImpulsesExcitationDerivative(t,a,b,c,d,n,order)
%-------------------------------------------------------------------------%
% a: Amplitude
% b: Position of the centre of the peak for the first impulse
% c: Standard deviation (controls the width of the impulse) --> only 1% of
%    the amplitude for t = b +-3c
% d: Distance between peaks
% n: Number of impulses
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
    if order == 1
        f{i} = @(t) (a*exp(-(bf - t).^2/(2*c^2)).*(2*bf - 2*t))/(2*c^2);
    elseif order == 2
        f{i} = @(t) -(a*exp(-(bf - t).^2/(2*c^2)).*...
                    (- bf^2 + 2*bf*t + c^2 - t.^2))/c^4;
    elseif order == 3
        f{i} = @(t) -(a*exp(-(bf - t).^2/(2*c^2)).*(bf - t).*...
                    (- bf^2 + 2*bf*t + 3*c^2 - t.^2))/c^6;
    elseif order == 4
        f{i} = @(t) (a*exp(-(bf - t).^2/(2*c^2)).*...
                    (bf^4 - 4*bf^3*t - 6*bf^2*c^2 + 6*bf^2*t.^2 + ...
                    12*bf*c^2*t - 4*bf*t.^3 + 3*c^4 - 6*c^2*t.^2 + ...
                    t.^4))/c^8;
    elseif order == 5
        f{i} = @(t) (a*exp(-(bf - t).^2/(2*c^2)).*(bf - t).*...
                    (bf^4 - 4*bf^3*t - 10*bf^2*c^2 + 6*bf^2*t.^2 + ...
                     20*bf*c^2*t - 4*bf*t.^3 + 15*c^4 - 10*c^2*t.^2 + ...
                     t.^4))/c^10;
    end
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