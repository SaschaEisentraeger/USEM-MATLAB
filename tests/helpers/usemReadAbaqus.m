function S = usemReadAbaqus(inpFile, p)
%USEMREADABAQUS Read an ABAQUS input file like Start_Sim.m does.
%
%   S = usemReadAbaqus(inpFile, p) returns the model struct of
%   ReadModel(inpFile, p) without progress output: mesh, materials, step,
%   loads and boundary conditions, lifted to spectral elements of degree
%   p (polynomial degree structure with pxi, peta, ndist). Kept for the
%   studies that call it; new code can call ReadModel directly.

S = ReadModel(inpFile, p, false);

end
