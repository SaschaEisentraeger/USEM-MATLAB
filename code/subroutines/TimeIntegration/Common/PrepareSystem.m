function [K,M,C,Fdt] = PrepareSystem(K,M,C,F_dt,NumDoFActive)
%PREPARESYSTEM Common preprocessing of the time integration methods.
%
%   [K, M, C, Fdt] = PrepareSystem(K, M, C, F_dt, NumDoFActive)
%   - restores full matrices if K is stored as a triangular matrix (one
%     half of a symmetric matrix; the unsymmetric matrices of this code
%     are always full),
%   - returns the time-dependent loads as a struct with the fields
%     Neumann, K, C, M (multiplied by s, s, ds/dt and d2s/dt2 of the
%     excitation signal s). A plain vector F_dt is taken as the Neumann
%     part; the other parts are then zero.

if (istriu(K) || istril(K)) && ~isscalar(K)
    K = K + K';
    M = M + M';
    C = C + C';
end

if isstruct(F_dt)
    Fdt.Neumann = F_dt.Neumann;
    Fdt.K       = F_dt.K;
    Fdt.C       = F_dt.C;
    Fdt.M       = F_dt.M;
else
    Fdt.Neumann = F_dt;
    Fdt.K       = zeros(NumDoFActive,1);
    Fdt.C       = zeros(NumDoFActive,1);
    Fdt.M       = zeros(NumDoFActive,1);
end

end
