function [figure_handles, output_data] = plot_PD_predictions(folder,options)
% Plot PD predictions based on different coordinate frames for neurons. Use
% only RW tasks for this function.

figure_handles = [];
output_data = struct;


%% load bdf
opensim_prefix = options.opensim_prefix;

bdf = get_nev_mat_data([folder options.prefix],options.labnum);

% extract separate workspaces
if(isfield(options,'xlim_PM') && isfield(options,'ylim_PM') && isfield(options,'xlim_DL') && isfield(options,'ylim_DL'))
    xlim_PM = options.xlim_PM;
    xlim_DL = options.xlim_DL;
    ylim_PM = options.ylim_PM;
    ylim_DL = options.ylim_DL;
else
    warning('No workspace limits specified; using defaults of [-10 10] and [-55 -35]')
    xlim_PM = [-10 0];
    xlim_DL = [0 10];
    ylim_PM = [-55 -45];
    ylim_DL = [-45 -35];
end

bdf.meta.task = 'RW';
opts.binsize=0.05;
opts.offset=-.015;
opts.do_trial_table=1;
opts.do_firing_rate=1;
bdf=postprocess_bdf(bdf,opts);

%% load joint kinematics
joint_pos_mat = csvread([folder 'Analysis/' opensim_prefix '_Kinematics_q.sto'],11,0);
bdf.joint_pos = array2table(joint_pos_mat,'VariableNames',{'time','shoulder_adduction','shoulder_rotation','shoulder_flexion','elbow_flexion','radial_pronation','wrist_flexion','wrist_abduction'});
clear joint_kin_mat

%% load muscle kinematics
muscle_pos_mat = csvread([folder 'Analysis/' opensim_prefix '_MuscleAnalysis_Length.sto'],12,0);
bdf.muscle_pos = array2table(muscle_pos_mat,'VariableNames',{'time','abd_poll_longus','anconeus','bicep_lh','bicep_sh','brachialis','brachioradialis','coracobrachialis','deltoid_ant','deltoid_med','deltoid_pos','dorsoepitrochlearis','ext_carpi_rad_longus','ext_carp_rad_brevis','ext_carpi_ulnaris','ext_digitorum','ext_digiti','ext_indicis','flex_carpi_radialis','flex_carpi_ulnaris','flex_digit_profundus','flex_digit_superficialis','flex_poll_longus','infraspinatus','lat_dorsi_sup','lat_dorsi_cen','lat_dorsi_inf','palmaris_longus','pectoralis_sup','pectoralis_inf','pronator_quad','pronator_teres','subscapularis','supinator','supraspinatus','teres_major','teres_minor','tricep_lat','tricep_lon','tricep_sho'});
clear muscle_kin_mat

%% filter kinematics (PROBABLY UNNECESSARY, AND UNUSED IN THIS SCRIPT)
% nyq_frq = 1/2*1/mode(diff(joint_kin.time));
% cuttoff = 10/nyq_frq;
% [b,a] = butter(4,cuttoff);
% 
% % filter joint angles
% joint_kin_filt = joint_kin;
% for i=2:size(joint_kin,2)
%     joint_kin{:,i} = filtfilt(b,a,joint_kin{:,i});
% end
% 
% % filter muscle lengths
% muscle_kin_filt = muscle_kin;
% for i=2:size(muscle_kin,2)
%     muscle_kin{:,i} = filtfilt(b,a,muscle_kin{:,i});
% end

%% choose fake neural weights
num_neurons = 100;
joint_weights = randn(size(bdf.joint_pos,2)-1,num_neurons);
muscle_weights = randn(size(bdf.muscle_pos,2)-1,num_neurons);
end_weights = randn(size(bdf.vel,2)-1,num_neurons);
% muscle_weights = eye(num_neurons);

%% run simulations
[joint_tuning_PM,joint_curve_PM] = get_predicted_tuning(bdf,xlim_PM,ylim_PM,joint_weights,'joint','linear','linear');
[joint_tuning_DL,joint_curve_DL] = get_predicted_tuning(bdf,xlim_DL,ylim_DL,joint_weights,'joint','linear','linear');

[muscle_tuning_PM,muscle_curve_PM] = get_predicted_tuning(bdf,xlim_PM,ylim_PM,muscle_weights,'muscle','linear','linear');
[muscle_tuning_DL,muscle_curve_DL] = get_predicted_tuning(bdf,xlim_DL,ylim_DL,muscle_weights,'muscle','linear','linear');

[end_tuning_PM,end_curve_PM] = get_predicted_tuning(bdf,xlim_PM,ylim_PM,end_weights,'endpoint','linear','linear');
[end_tuning_DL,end_curve_DL] = get_predicted_tuning(bdf,xlim_DL,ylim_DL,end_weights,'endpoint','linear','linear');

%% Iris plots
% joint_tuning_PM = struct2table(joint_tuning_PM);
% joint_tuning_DL = struct2table(joint_tuning_DL);
% muscle_tuning_PM = struct2table(muscle_tuning_PM);
% muscle_tuning_DL = struct2table(muscle_tuning_DL);
% end_tuning_PM = struct2table(end_tuning_PM);
% end_tuning_DL = struct2table(end_tuning_DL);

h = figure('name','joint_PD_diff');
figure_handles = [figure_handles;h];
iris_plot(joint_tuning_PM,joint_tuning_DL)
title('Plot of PD changes (joint)')

h = figure('name','muscle_PD_diff');
figure_handles = [figure_handles;h];
iris_plot(muscle_tuning_PM,muscle_tuning_DL)
title('Plot of PD changes (muscle)')

h = figure('name','end_PD_diff');
figure_handles = [figure_handles;h];
iris_plot(end_tuning_PM,end_tuning_DL)
title('Plot of PD changes (end)')

%% Check kinematics
% move_corr = bdf.vel(:,2:3);
% dir = atan2(move_corr(:,2),move_corr(:,1));
% spd = sqrt(sum(move_corr.^2,2));
% 
% % bin directions
% dir_bins = round(dir/(pi/4))*(pi/4);
% dir_bins(dir_bins==-pi) = pi;
% 
% % find baseline move_corr...somehow
% 
% 
% % average firing rates for directions
% bins = -3*pi/4:pi/4:pi;
% bins = bins';
% 
% full_binned_FR = [];
% groups = [];
% for i = 1:length(bins)
%     vel_in_bin = joint_vel_PM{dir_bins==bins(i),2:end};
%     spd_in_bin = spd(dir_bins==bins(i));
%     
%     % Mean binned FR has normal-looking distribution (checked with
%     % bootstrapping)
%     binned_FR(i,:) = mean(vel_in_bin); % mean firing rate
%     binned_spd(i,:) = mean(spd_in_bin); % mean speed
%     binned_stderr(i,:) = std(vel_in_bin)/sqrt(length(vel_in_bin)); % standard error
%     binned_spd_err(i,:) = std(spd_in_bin)/sqrt(length(spd_in_bin)); % standard error of speed
%     tscore = tinv(0.975,length(vel_in_bin)-1); % t-score for 95% CI
%     binned_CI_high(i,:) = binned_FR(i,:)+tscore*binned_stderr(i,:); %high CI
%     binned_CI_low(i,:) = binned_FR(i,:)-tscore*binned_stderr(i,:); %low CI
%     binned_CI_high_spd(i,:) = binned_spd(i,:)+tscore*binned_spd_err(i,:); %high CI
%     binned_CI_low_spd(i,:) = binned_spd(i,:)-tscore*binned_spd_err(i,:); %low CI
% end
% 
% % plot tuning curves
% % plot speed curves
% h = figure('name','Binned Speed');
% figure_handles = [figure_handles;h];
% polar(repmat(bins,2,1),repmat(binned_spd,2,1))
% 
% % plot confidence intervals 
% th_fill = [flipud(bins); bins(end); bins(end); bins];
% r_fill = [flipud(binned_CI_high_spd); binned_CI_high_spd(end); binned_CI_low_spd(end); binned_CI_low_spd];
% [x_fill,y_fill] = pol2cart(th_fill,r_fill);
% patch(x_fill,y_fill,[0 0 1],'facealpha',0.3,'edgealpha',0);

output_data.joint_tuning_PM = joint_tuning_PM;
output_data.joint_tuning_DL = joint_tuning_DL;
output_data.muscle_tuning_PM = muscle_tuning_PM;
output_data.muscle_tuning_DL = muscle_tuning_DL;
output_data.end_tuning_PM = end_tuning_PM;
output_data.end_tuning_DL = end_tuning_DL;

output_data.joint_curve_PM = joint_curve_PM;
output_data.joint_curve_DL = joint_curve_DL;
output_data.muscle_curve_PM = muscle_curve_PM;
output_data.muscle_curve_DL = muscle_curve_DL;
output_data.end_curve_PM = end_curve_PM;
output_data.end_curve_DL = end_curve_DL;
end