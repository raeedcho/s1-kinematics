function [targets, movements] = randomwalk_sim(boundaries, num_targets, target_size, spacing, hold_time, time_step)

targets = generate_targets(num_targets,boundaries,spacing);

old_state = [0 0 0 0];
old_time = 0;
states = [];
times = [];
for i = 1:length(targets);
    [time,state] = move_point_to_target(targets(i,:),target_size,hold_time,old_state);
    times = [times(1:end-1);time+old_time];
    states = [states(1:end-1,:);state];
    old_time = times(end);
    old_state = states(end,:);
end

% interpolate state vectors
interp_times = (0:time_step:times(end))';
states_interp = interp1(times,states,interp_times);
movement_dir = atan2(states_interp(:,4),states_interp(:,3));

movements = [interp_times states_interp];