function a0 = InitialAcceleration(solveM,F,Fdt,K,C,u0,v0,SignalType,ExciteParam)
%INITIALACCELERATION Acceleration at t = 0 from the equation of motion.
%
%   a0 = InitialAcceleration(solveM, F, Fdt, K, C, u0, v0, SignalType,
%   ExciteParam) returns M\(f(0) - C*v0 - K*u0), f(0) the load at t = 0
%   (see LoadAtTime). solveM applies the inverse of the mass matrix, e.g.
%   @(r) Minv.*r for a diagonal mass matrix or @(r) dM\r with a
%   decomposition dM.

[s,ds,dds] = SignalHistory(SignalType,ExciteParam,0);
a0 = solveM(F - C*v0 - K*u0 + Fdt.Neumann*s + Fdt.K*s + Fdt.C*ds + ...
            Fdt.M*dds);

end
