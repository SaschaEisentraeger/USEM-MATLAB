function y = AhemSolveEffective(b,flag_diag,Einv,dE)
%AHEMSOLVEEFFECTIVE Solve with the effective matrix E of the adaptive hybrid explicit method.

% Solve E*y = b with the diagonal or the decomposed effective matrix
if flag_diag
    y = Einv.*b;
else
    y = dE\b;
end
end
