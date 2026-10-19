function [forceHistory] = ForcePadeTimeIntegration(porder,SignalType,...
                          ExciteParam,ns,dt)
%-------------------------------------------------------------------------%
% University of New South Wales (UNSW Sydney)
% Written by Sascha Eisentraeger (12/2020)
%-------------------------------------------------------------------------%
% Input:
% porder:           order of Pade expansion: = 11; 22; 33; 44; 55
% SignalType:       signal type (depending on your own definition)
% ExciteParam:      parameter for the excitation (corresponding to the type
%                   of excitation signal)
% ns:               number of time steps
% dt:               time steps size
%
% Output:
% forceHistory:     external loading; the force is computed at different
%                   points within the time step --> Use GLL-points
%                   [time amp_1 amp_2 ... amp_porder+1];
%               
%-------------------------------------------------------------------------%

%% Compute the force amplitude values at the different times
% Time vector
t_vec = (0:1:ns)'*dt;
% Calculate the force at the end of each time step
amp = ExcitationSignal(SignalType,ExciteParam,t_vec)';

switch porder
    case 11
        % Output
        forceHistory = [t_vec,amp];
    case 22
        % GLL-points in the interval [-1,1]
        xi = lglnodes(2);
        xi = flip(xi);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
        t_vec1 = (s(2):1:ns+1)*dt;
        amp_1 = ExcitationSignal(SignalType,ExciteParam,t_vec1)';
        % Output
        forceHistory = [t_vec,amp,amp_1];
    case 33
        % GLL-points in the interval [-1,1]
        xi = lglnodes(3);
        xi = flip(xi);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
        t_vec1 = (s(2):1:ns+1)*dt;
        t_vec2 = (s(3):1:ns+1)*dt;
        amp_1 = ExcitationSignal(SignalType,ExciteParam,t_vec1)';
        amp_2 = ExcitationSignal(SignalType,ExciteParam,t_vec2)';
        % Output
        forceHistory = [t_vec,amp,amp_1,amp_2];
    case 44
        % GLL-points in the interval [-1,1]
        xi = lglnodes(4);
        xi = flip(xi);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
        t_vec1 = (s(2):1:ns+1)*dt;
        t_vec2 = (s(3):1:ns+1)*dt;
        t_vec3 = (s(4):1:ns+1)*dt;
        amp_1 = ExcitationSignal(SignalType,ExciteParam,t_vec1)';
        amp_2 = ExcitationSignal(SignalType,ExciteParam,t_vec2)';
        amp_3 = ExcitationSignal(SignalType,ExciteParam,t_vec3)';
        % Output
        forceHistory = [t_vec,amp,amp_1,amp_2,amp_3];
    case 55
        % GLL-points in the interval [-1,1]
        xi = lglnodes(5);
        xi = flip(xi);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
        t_vec1 = (s(2):1:ns+1)*dt;
        t_vec2 = (s(3):1:ns+1)*dt;
        t_vec3 = (s(4):1:ns+1)*dt;
        t_vec4 = (s(5):1:ns+1)*dt;
        amp_1 = ExcitationSignal(SignalType,ExciteParam,t_vec1)';
        amp_2 = ExcitationSignal(SignalType,ExciteParam,t_vec2)';
        amp_3 = ExcitationSignal(SignalType,ExciteParam,t_vec3)';
        amp_4 = ExcitationSignal(SignalType,ExciteParam,t_vec4)';
        % Output
        forceHistory = [t_vec,amp,amp_1,amp_2,amp_3,amp_4];
    otherwise
        % Output
        forceHistory = [t_vec,amp];
end

end