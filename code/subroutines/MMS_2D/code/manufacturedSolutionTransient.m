function [fh,fSym] = manufacturedSolutionTransient(fFun,saveFile)
% Compute body load fFun based on an analytical solution defined by a cell
% of function handles in varFuns.
syms x y t real

if nargin == 1
    saveFile = false;
end

fSym = fFun(x,y,t);
fh = matlabFunction(fSym,'Vars',[x,y,t]);	% function handle

if saveFile
    dir = ['.' filesep 'subroutines' filesep 'BoundaryConditions' filesep 'Neumann' filesep];
    funcName = 'BodyLoadFuncTransient.m';
    matlabFunction(fSym,'File',[dir funcName]);
end

end