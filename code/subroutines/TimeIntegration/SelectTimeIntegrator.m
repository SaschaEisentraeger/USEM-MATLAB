function [method,extra] = SelectTimeIntegrator(code,opt)
%SELECTTIMEINTEGRATOR Time integration routine and its trailing arguments for a method code.

% Look the method up in TimeIntegrationMethods
methods = TimeIntegrationMethods();
k = find([methods.code] == code,1);
if isempty(k)
    % Unknown code: Newmark method (as before), but say so
    warning('Solver_Transient:unknownMethod', ['Unknown time integration ' ...
            'method %d; the Newmark method (301) is used.'],code);
    k = find([methods.code] == 301,1);
end
method = methods(k).handle;
extra = methods(k).args(opt);
disp(['Time integration: ' methods(k).name]);
end
