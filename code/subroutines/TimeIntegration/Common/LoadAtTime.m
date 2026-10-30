function Ft = LoadAtTime(F,Fdt,s,ds,dds)
%LOADATTIME Load vector at one instant of time.
%
%   Ft = LoadAtTime(F, Fdt, s, ds, dds) returns
%       F + Fdt.Neumann*s + Fdt.K*s + Fdt.C*ds + Fdt.M*dds
%   for the values s, ds, dds of the excitation signal and its first two
%   derivatives at that instant (SignalHistory). Fdt.K, Fdt.C and Fdt.M
%   are the lifting terms of time-dependent prescribed displacements.

Ft = F + Fdt.Neumann*s + Fdt.K*s + Fdt.C*ds + Fdt.M*dds;

end
