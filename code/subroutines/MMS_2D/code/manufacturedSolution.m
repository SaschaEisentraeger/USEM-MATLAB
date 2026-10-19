function [fh,fSym] = manufacturedSolution(varFuns,fFun,saveFile,funcName)
% Compute body load fFun based on an analytical solution defined by a cell
% of function handles in varFuns.
syms x y real

if nargin < 4
	funcName = 'BodyLoadFunc.m';
	if nargin == 2
		saveFile = false;
	end
end

varSym{numel(varFuns)} = {};
for i = 1:numel(varFuns)
    varSym{i} = str2sym(char(varFuns{i}));
end
fSym = fFun(varSym{:});                     % symbolic expression

fh = matlabFunction(fSym,'Vars',[x,y]);     % function handle

if saveFile
    vars = symvar(fSym);
    dir = ['.' filesep 'subroutines' filesep 'BoundaryConditions' filesep 'Neumann' filesep];
    if isempty(vars)
        matlabFunction(fSym,'File',[dir funcName], 'Vars',[x,y]);
    else
        matlabFunction(fSym,'File',[dir funcName]);
    end
end

end