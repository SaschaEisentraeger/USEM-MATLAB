function xdd = KimStageAcceleration(solveM,F,Fdt,K,C,x,xd,SignalType,...
                                 ExciteParam,t)
%KIMSTAGEACCELERATION Acceleration at a stage of the higher-order explicit Kim method.

% Acceleration from the equation of motion at the collocation time t of a
% stage, M\(f(t) - C*xd - K*x); same summation order as
% InitialAcceleration (which is the case t = 0)
[s,ds,dds] = SignalHistory(SignalType,ExciteParam,t);
xdd = solveM(F - C*xd - K*x + Fdt.Neumann*s + Fdt.K*s + Fdt.C*ds + ...
             Fdt.M*dds);
end
