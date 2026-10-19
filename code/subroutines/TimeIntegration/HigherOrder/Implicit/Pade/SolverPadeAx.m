function y = SolverPadeAx(dM,K,C,x)
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by Chongmin Song (12/2020)
%-------------------------------------------------------------------------%

n = length(x)/2;
tmp = C*x(1:n) + K*x(n+1:end);
y = [-dM\tmp ; x(1:n)];

end

