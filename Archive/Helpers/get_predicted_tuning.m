function [modeled_tuning,curves] = get_predicted_tuning(bdf,xlim,ylim,weights,frame,gen_model,fit_model)
% GET_MUSCLE_PREDICTED_TUNING Get predicted tuning, assuming neurons draw
% directly from muscle inputs. Assumes GLM for generative model and tuning.
%   Inputs:
%       BDF - postprocessed bdf, including OpenSim kinematics as joint_pos
%       and muscle_pos
%       BOTTOM_LEFT - pair of numbers indicating bottom left corner of
%       workspace. Leave blank to use whole workspace
%       TOP_RIGHT - pair of numbers indicating top right corner of
%       workspace. Leave blank to use whole workspace
%       WEIGHTS - weights on 39 muscles. Includes first column of baseline
%       weights
%       FRAME - Coordinate frame: 'joint', 'muscle', 'endpoint', or
%       'egocentric'
%       GEN_MODEL - Generative model: 'linear' or 'GLM'
%       FIT_MODEL - Fitting model: 'linear' or 'GLM'

%% setup
num_neurons = size(weights,2);

% choose correct frame
if(strcmp(frame,'muscle'))
    frame_pos = bdf.muscle_pos;
elseif(strcmp(frame,'joint'))
    frame_pos = bdf.joint_pos;
elseif strcmp(frame,'endpoint')
    frame_pos = array2table(bdf.pos,'VariableNames',{'time','x','y'});
%elseif strcmp(frame,'egocentric')
else
    error('Invalid frame chosen')
end

% kinematics
t = bdf.pos(:,1);
pos = bdf.pos(:,2:3);
vel = bdf.vel(:,2:3);
spd = sqrt(sum(vel.^2,2));
endpoint_kin = [pos vel spd];

% Get times in workspace
if(~isempty(xlim) && ~isempty(ylim))
    times = extract_workspace_times(bdf,xlim,ylim);
else
    warning('No workspace corners provided. Using full workspace')
    times = [t(1) t(end)];
end

%% get intrinsic velocities
% frame velocities
frame_vel = frame_pos;
for i=2:size(frame_pos,2)
    frame_vel{:,i} = gradient(frame_pos{:,i},frame_pos.time);
end

%% interpolate everything to whatever firing rate bin size is at
bin_times = bdf.units(1).FR(:,1);

% throw out times outside data
bin_times(bin_times<frame_vel.time(1) | bin_times>frame_vel.time(end)) = [];

endpoint_kin_interp = interp1(t,endpoint_kin,bin_times);
temp_arr = [bin_times interp1(frame_vel.time,frame_vel{:,2:end},bin_times)];
frame_vel_interp = array2table(temp_arr,'VariableNames',frame_vel.Properties.VariableNames);

%% frame kinematics for the workspace
frame_vel_sub = table;
for i = 1:size(times,1)
    reach_ind = frame_vel_interp.time>times(i,1) & frame_vel_interp.time<times(i,2);
    frame_vel_sub = [frame_vel_sub; frame_vel_interp(reach_ind,:)];
end

if(strcmp(gen_model,'linear'))
    % Linear model
%     sim_FR = [ones(size(frame_vel_sub,1),1) frame_vel_sub{:,2:end}]*weights;
    sim_FR = frame_vel_sub{:,2:end}*weights;
elseif(strcmp(gen_model,'GLM'))
    % GLM Poisson model
    sim_FR = exp([ones(size(frame_vel_sub,1),1) frame_vel_sub{:,2:end}]*weights);
else
    error('Invalid generative model')
end

% interpolate endpoint kinematics to muscle times
endpoint_kin_sub = interp1(t,endpoint_kin,frame_vel_sub.time);

%% Calculate simulated PDs
if(strcmp(fit_model,'linear'))
    % Linear model
    bootfunc = @(X,y) LinearModel.fit(X,y);
elseif(strcmp(fit_model,'GLM'))
    % GLM Poisson model
    bootfunc = @(X,y) GeneralizedLinearModel.fit(X,y,'Distribution','poisson');
else
    error('Invalid fitting model')
end

modeled_tuning = table;
tic;
for i = 1:num_neurons
    modeled_tuning(i,:) = struct2table(calc_PD_helper(bootfunc,endpoint_kin_sub,sim_FR(:,i)));
    disp(['Processed ' frame ' neuron ' num2str(i) ' (Time: ' num2str(toc) ')'])
end

%% Calculate simulated tuning curves
curves = table;
for i = 1:num_neurons
    curves(i,:) = struct2table(get_single_tuning_curve(endpoint_kin_sub(:,3:4),sim_FR(:,i)));
end
