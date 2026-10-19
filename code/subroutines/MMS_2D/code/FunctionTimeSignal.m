function [df,ddf] = FunctionTimeSignal(f,saveFile)
% Time-dependent function f.
syms t real

if nargin == 1
    saveFile = false;
end

F   = str2sym(char(f));
dF  = diff(F,t);
ddF = diff(F,t,2);

df  = matlabFunction(dF,'Vars',t);     % function handle
ddf = matlabFunction(ddF,'Vars',t);

if saveFile
    dir = ['.' filesep 'subroutines' filesep 'BoundaryConditions' filesep 'Dirichlet' filesep];
    funcName = 'TimeSignalFunc.m';
    matlabFunction(F,'File',[dir funcName]);
    funcName = 'dTimeSignalFunc.m';
    matlabFunction(dF,'File',[dir funcName]);
    funcName = 'ddTimeSignalFunc.m';
    matlabFunction(ddF,'File',[dir funcName]);
end

end