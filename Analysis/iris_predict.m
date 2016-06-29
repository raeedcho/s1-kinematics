function [figure_handles, output_data] = iris_predict(folder,options)
% Plot PD predictions based on different coordinate frames for neurons. Use
% only RW tasks for this function.

figure_handles = [];
output_data = struct;


%% load bdf
bdf_PM = get_nev_mat_data([folder options.prefix '_PM'],options.labnum);
bdf_DL = get_nev_mat_data([folder options.prefix '_DL'],options.labnum);

bdf_PM.meta.task = 'RW';
opts.binsize=0.05;
opts.offset=-.015;
opts.do_trial_table=1;
opts.do_firing_rate=1;
bdf_PM=postprocess_bdf(bdf_PM,opts);

bdf_DL.meta.task = 'RW';
opts.binsize=0.05;
opts.offset=-.015;
opts.do_trial_table=1;
opts.do_firing_rate=1;
bdf_DL=postprocess_bdf(bdf_DL,opts);

%% load joint kinematics
joint_pos_mat = dlmread([folder 'Analysis/' options.prefix '_PM_Kinematics_q.sto'],'\t',11,0);
bdf_PM.joint_pos = array2table(joint_pos_mat,'VariableNames',{'time','shoulder_adduction','shoulder_rotation','shoulder_flexion','elbow_flexion','radial_pronation','wrist_flexion','wrist_abduction'});
joint_pos_mat = dlmread([folder 'Analysis/' options.prefix '_DL_Kinematics_q.sto'],'\t',11,0);
bdf_DL.joint_pos = array2table(joint_pos_mat,'VariableNames',{'time','shoulder_adduction','shoulder_rotation','shoulder_flexion','elbow_flexion','radial_pronation','wrist_flexion','wrist_abduction'});
clear joint_kin_mat

%% load muscle kinematics
muscle_pos_mat = dlmread([folder 'Analysis/' options.prefix '_PM_MuscleAnalysis_Length.sto'],'\t',12,0);
bdf_PM.muscle_pos = array2table(muscle_pos_mat,'VariableNames',{'time','abd_poll_longus','anconeus','bicep_lh','bicep_sh','brachialis','brachioradialis','coracobrachialis','deltoid_ant','deltoid_med','deltoid_pos','dorsoepitrochlearis','ext_carpi_rad_longus','ext_carp_rad_brevis','ext_carpi_ulnaris','ext_digitorum','ext_digiti','ext_indicis','flex_carpi_radialis','flex_carpi_ulnaris','flex_digit_profundus','flex_digit_superficialis','flex_poll_longus','infraspinatus','lat_dorsi_sup','lat_dorsi_cen','lat_dorsi_inf','palmaris_longus','pectoralis_sup','pectoralis_inf','pronator_quad','pronator_teres','subscapularis','supinator','supraspinatus','teres_major','teres_minor','tricep_lat','tricep_lon','tricep_sho'});
muscle_pos_mat = dlmread([folder 'Analysis/' options.prefix '_DL_MuscleAnalysis_Length.sto'],'\t',12,0);
bdf_DL.muscle_pos = array2table(muscle_pos_mat,'VariableNames',{'time','abd_poll_longus','anconeus','bicep_lh','bicep_sh','brachialis','brachioradialis','coracobrachialis','deltoid_ant','deltoid_med','deltoid_pos','dorsoepitrochlearis','ext_carpi_rad_longus','ext_carp_rad_brevis','ext_carpi_ulnaris','ext_digitorum','ext_digiti','ext_indicis','flex_carpi_radialis','flex_carpi_ulnaris','flex_digit_profundus','flex_digit_superficialis','flex_poll_longus','infraspinatus','lat_dorsi_sup','lat_dorsi_cen','lat_dorsi_inf','palmaris_longus','pectoralis_sup','pectoralis_inf','pronator_quad','pronator_teres','subscapularis','supinator','supraspinatus','teres_major','teres_minor','tricep_lat','tricep_lon','tricep_sho'});
clear muscle_kin_mat

%% choose fake neural weights
if(isfield(options,'joint_weights'))
    joint_weights = options.joint_weights;
else
    joint_weights = randn(size(bdf_PM.joint_pos,2)-1,100);
end
if(isfield(options,'muscle_weights'))
    muscle_weights = options.muscle_weights;
else
    muscle_weights = randn(size(bdf_PM.muscle_pos,2)-1,100);
end
% joint_weights = eye(size(bdf_PM.joint_pos,2));
% muscle_weights = eye(num_neurons);

%% predict tunings
[joint_tuning_PM,joint_curves_PM] = get_predicted_tuning(bdf_PM,[],[],joint_weights,'joint','linear','linear');
[joint_tuning_DL,joint_curves_DL] = get_predicted_tuning(bdf_DL,[],[],joint_weights,'joint','linear','linear');

[muscle_tuning_PM,muscle_curves_PM] = get_predicted_tuning(bdf_PM,[],[],muscle_weights,'muscle','linear','linear');
[muscle_tuning_DL,muscle_curves_DL] = get_predicted_tuning(bdf_DL,[],[],muscle_weights,'muscle','linear','linear');

%% iris plots
h = figure('name','joint_PD_diff');
figure_handles = [figure_handles;h];
iris_plot(joint_tuning_PM,joint_tuning_DL)
title('Plot of PD changes (joint)')

h = figure('name','muscle_PD_diff');
figure_handles = [figure_handles;h];
iris_plot(muscle_tuning_PM,muscle_tuning_DL)
title('Plot of PD changes (muscle)')

%% tuning curve plots
% for i = 1:num_neurons
%     figure
%     maxFR = max([joint_curves_PM.FR{i};joint_curves_DL.FR{i}]);
%     plot_tuning(joint_tuning_PM(i,:),joint_curves_PM(i,:),maxFR,[0.6 0.5 0.7])
%     plot_tuning(joint_tuning_DL(i,:),joint_curves_PM(i,:),maxFR,[1 0 0])
% end

%% Real iris plot 


% h = figure('name','real_PD_diff');
% figure_handles = [figure_handles;h];
% iris_plot(real_tuning_PM,real_tuning_DL)
% title('Plot of PD changes (real)')

%% Check kinematics
% move_corr = endpoint_kin_sim_PM(:,3:4);
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

%% Save data
output_data.bdf_PM = bdf_PM;
output_data.bdf_DL = bdf_DL;
output_data.joint_weights = joint_weights;
output_data.muscle_weights = muscle_weights;

output_data.joint_tuning_PM = joint_tuning_PM;
output_data.joint_tuning_DL = joint_tuning_DL;
output_data.muscle_tuning_PM = muscle_tuning_PM;
output_data.muscle_tuning_DL = muscle_tuning_DL;

output_data.joint_curves_PM = joint_curves_PM;
output_data.joint_curves_DL = joint_curves_DL;
output_data.muscle_curves_PM = muscle_curves_PM;
output_data.muscle_curves_DL = muscle_curves_DL;

end