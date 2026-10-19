function [uh] = FunctionDisplacement(ux,uy,saveFile)
% Vector-valued function handle for the displacement field.

if nargin == 2
    saveFile = false;
end

if nargin(ux) == 2
	syms x y real
    u = [str2sym(char(ux)), str2sym(char(uy))];
    uh = matlabFunction(u,'Vars',[x,y]);        % function handle
    funcName = 'DisplacementFunc.m';
elseif nargin(ux) == 3
    syms x y t real
    uxSym = ux(x,y,t);
    uySym = uy(x,y,t);
    u = [uxSym, uySym];
    uh = matlabFunction(u,'Vars',[x,y,t]);    % function handle
    funcName = 'DisplacementFuncTransient.m';
end

if saveFile
    dir = ['.' filesep 'subroutines' filesep 'BoundaryConditions' filesep 'Dirichlet' filesep];
    matlabFunction(u,'File',[dir funcName]);
end

end