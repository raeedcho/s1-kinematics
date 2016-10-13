function [time,state] = move_point_to_target(target_pos,target_size,hold_time,init_state)
% MOVE_POINT_TO_TARGET moves point mass to a target like a damped spring
%   Outputs-
%       time - list of times for the state
%       state - state vector in columns, each row is new time point ([x y
%       x' y'])
%   Inputs-
%       target_pos - position of single target ([x y])
%       target_size - size of target, determines when to stop
%       init_state - position and velocity of point initially ([x y x' y'])

if length(target_pos)~=2
    error('Give single target position in [x y]')
end

% dynamics constants
m = 0.1; %kg
beta = 0.5; %N/(cm/s)
k = 2; %N/cm

dynamics_mat = [0 0 1 0; 0 0 0 1; -k/m 0 -beta/m 0; 0 -k/m 0 -beta/m];

% simulate
[full_time,full_state] = ode45(@(t,x)(dynamics_mat*x+k/m*[0;0;target_pos(:)]),[0 5],init_state);

% cut off after target is reached (circle targets)
pos = full_state(:,1:2);
pos_err = repmat(target_pos,length(pos),1)-pos;
enter_idx = find(sum(pos_err.^2,2)<target_size.^2,1,'first');
enter_time = full_time(enter_idx);
movement_idx = full_time<enter_time+hold_time;
time = full_time(movement_idx);
state = full_state(movement_idx,:);