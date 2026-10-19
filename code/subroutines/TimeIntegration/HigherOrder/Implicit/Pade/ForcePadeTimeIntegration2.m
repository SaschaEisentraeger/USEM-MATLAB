function [forceHistory] = ForcePadeTimeIntegration2(porder,SignalType,...
                          ExciteParam,ns,dt,type)
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
% type:				type of interpolation nodes
%
% Output:
% forceHistory:     external loading; the force is computed at different
% 					points within the time step;
% 					[time_1 amp_1 time_2 amp_2 ...];
%               
%-------------------------------------------------------------------------%

switch porder
    case 11
		% Interpolation points in the interval [-1,1]
        xi = InterpolationPointsForce(1,type);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
		% Calculate the force values within each time step
		t_vec1 = (s(1):1:ns)'*dt;
        t_vec2 = (s(2):1:ns)'*dt;
        amp_1 = ExcitationSignal(SignalType,ExciteParam,t_vec1)';
		amp_2 = ExcitationSignal(SignalType,ExciteParam,t_vec2)';
        % Output
        forceHistory = [t_vec1(1:ns),amp_1(1:ns),...
                        t_vec2(1:ns),amp_2(1:ns)];
    case 22
        % Interpolation points in the interval [-1,1]
        xi = InterpolationPointsForce(2,type);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
		t_vec1 = (s(1):1:ns)'*dt;
        t_vec2 = (s(2):1:ns)'*dt;
		t_vec3 = (s(3):1:ns)'*dt;
        amp_1 = ExcitationSignal(SignalType,ExciteParam,t_vec1)';
		amp_2 = ExcitationSignal(SignalType,ExciteParam,t_vec2)';
		amp_3 = ExcitationSignal(SignalType,ExciteParam,t_vec3)';
        % Output
        forceHistory = [t_vec1(1:ns),amp_1(1:ns),...
                        t_vec2(1:ns),amp_2(1:ns),...
                        t_vec3(1:ns),amp_3(1:ns)];
    case 33
        % Interpolation points in the interval [-1,1]
        xi = InterpolationPointsForce(3,type);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
		t_vec1 = (s(1):1:ns)'*dt;
        t_vec2 = (s(2):1:ns)'*dt;
		t_vec3 = (s(3):1:ns)'*dt;
		t_vec4 = (s(4):1:ns)'*dt;
        amp_1 = ExcitationSignal(SignalType,ExciteParam,t_vec1)';
		amp_2 = ExcitationSignal(SignalType,ExciteParam,t_vec2)';
		amp_3 = ExcitationSignal(SignalType,ExciteParam,t_vec3)';
		amp_4 = ExcitationSignal(SignalType,ExciteParam,t_vec4)';
        % Output
        forceHistory = [t_vec1(1:ns),amp_1(1:ns),...
                        t_vec2(1:ns),amp_2(1:ns),...
                        t_vec3(1:ns),amp_3(1:ns),...
                        t_vec4(1:ns),amp_4(1:ns)];
    case 44
        % Interpolation points in the interval [-1,1]
        xi = InterpolationPointsForce(4,type);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
		t_vec1 = (s(1):1:ns)'*dt;
        t_vec2 = (s(2):1:ns)'*dt;
		t_vec3 = (s(3):1:ns)'*dt;
		t_vec4 = (s(4):1:ns)'*dt;
		t_vec5 = (s(5):1:ns)'*dt;
        amp_1 = ExcitationSignal(SignalType,ExciteParam,t_vec1)';
		amp_2 = ExcitationSignal(SignalType,ExciteParam,t_vec2)';
		amp_3 = ExcitationSignal(SignalType,ExciteParam,t_vec3)';
		amp_4 = ExcitationSignal(SignalType,ExciteParam,t_vec4)';
		amp_5 = ExcitationSignal(SignalType,ExciteParam,t_vec5)';
        % Output
        forceHistory = [t_vec1(1:ns),amp_1(1:ns),...
                        t_vec2(1:ns),amp_2(1:ns),...
                        t_vec3(1:ns),amp_3(1:ns),...
                        t_vec4(1:ns),amp_4(1:ns),...
                        t_vec5(1:ns),amp_5(1:ns)];
    case 55
        % Interpolation points in the interval [-1,1]
        xi = InterpolationPointsForce(5,type);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
        % Calculate the force values within each time step
		t_vec1 = (s(1):1:ns)'*dt;
        t_vec2 = (s(2):1:ns)'*dt;
		t_vec3 = (s(3):1:ns)'*dt;
		t_vec4 = (s(4):1:ns)'*dt;
		t_vec5 = (s(5):1:ns)'*dt;
		t_vec6 = (s(6):1:ns)'*dt;
        amp_1 = ExcitationSignal(SignalType,ExciteParam,t_vec1)';
		amp_2 = ExcitationSignal(SignalType,ExciteParam,t_vec2)';
		amp_3 = ExcitationSignal(SignalType,ExciteParam,t_vec3)';
		amp_4 = ExcitationSignal(SignalType,ExciteParam,t_vec4)';
		amp_5 = ExcitationSignal(SignalType,ExciteParam,t_vec5)';
		amp_6 = ExcitationSignal(SignalType,ExciteParam,t_vec6)';
        % Output
        forceHistory = [t_vec1(1:ns),amp_1(1:ns),...
                        t_vec2(1:ns),amp_2(1:ns),...
                        t_vec3(1:ns),amp_3(1:ns),...
                        t_vec4(1:ns),amp_4(1:ns),...
                        t_vec5(1:ns),amp_5(1:ns),...
                        t_vec6(1:ns),amp_6(1:ns)];
    otherwise
        % Interpolation points in the interval [-1,1]
        xi = InterpolationPointsForce(1,type);
        % Mapping to the interval [0,1]
        s  = 1/2*(xi + 1);
		% Calculate the force values within each time step
		t_vec1 = (s(1):1:ns)'*dt;
        t_vec2 = (s(2):1:ns)'*dt;
        amp_1 = ExcitationSignal(SignalType,ExciteParam,t_vec1)';
		amp_2 = ExcitationSignal(SignalType,ExciteParam,t_vec2)';
        % Output
        forceHistory = [t_vec1,amp_1,t_vec2,amp_2];
end

end