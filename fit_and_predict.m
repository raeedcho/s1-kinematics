% fit to each coordinate in one workspace and predict the PD in the next

%% read files in
folder = '/home/raeed/Projects/limblab/FSMRes/limblab/User_folders/Raeed/Arm Model/Data/Chips/experiment_20151203_RW_002/';
% folder = 'Z:\limblab\User_folders\Raeed\Arm Model\Data\Chips\experiment_20151203_RW_002\';
prefix =  'Chips_20151203_RW_002';
labnum =  6;
opensim_prefix = 'Chips_20151203_0-320';
time_window = [1 319];

bdf = get_nev_mat_data([folder prefix],labnum);

% extract separate workspaces
times_PM = extract_workspace_times(bdf,[-10 -55],[0 -45]);
times_DL = extract_workspace_times(bdf,[0 -45],[10 -35]);

bdf.meta.task = 'RW';
opts.binsize=0.05;
opts.offset=-.015;
opts.do_trial_table=1;
opts.do_firing_rate=1;
bdf=postprocess_bdf(bdf,opts);

joint_pos_mat = csvread([folder 'Analysis/' opensim_prefix '_Kinematics_q.sto'],11,0);
bdf.joint_pos = array2table(joint_pos_mat,'VariableNames',{'time','shoulder_adduction','shoulder_rotation','shoulder_flexion','elbow_flexion','radial_pronation','wrist_flexion','wrist_abduction'});
clear joint_pos_mat

muscle_pos_mat = csvread([folder 'Analysis/' opensim_prefix '_MuscleAnalysis_Length.sto'],12,0);
bdf.muscle_pos = array2table(muscle_pos_mat,'VariableNames',{'time','abd_poll_longus','anconeus','bicep_lh','bicep_sh','brachialis','brachioradialis','coracobrachialis','deltoid_ant','deltoid_med','deltoid_pos','dorsoepitrochlearis','ext_carpi_rad_longus','ext_carp_rad_brevis','ext_carpi_ulnaris','ext_digitorum','ext_digiti','ext_indicis','flex_carpi_radialis','flex_carpi_ulnaris','flex_digit_profundus','flex_digit_superficialis','flex_poll_longus','infraspinatus','lat_dorsi_sup','lat_dorsi_cen','lat_dorsi_inf','palmaris_longus','pectoralis_sup','pectoralis_inf','pronator_quad','pronator_teres','subscapularis','supinator','supraspinatus','teres_major','teres_minor','tricep_lat','tricep_lon','tricep_sho'});
clear muscle_pos_mat

%% get muscle PCA
[coeff,score,latent] = pca(bdf.muscle_pos{:,2:end});

% use first 5 PCs
muscle_scores = score(:,1:5);

%% Get kinematics
t = bdf.pos(:,1);
pos = bdf.pos(:,2:3);
vel = bdf.vel(:,2:3);
spd = sqrt(sum(vel.^2,2));
endpoint_kin = [pos vel spd];

% interpolate to FR and toss everything outside of time window
bin_times = bdf.units(1).FR(:,1);
bin_times( bin_times<time_window(1) | bin_times>time_window(2), :) = [];
endpoint_kin_full = interp1(t,endpoint_kin,bin_times);

% find time bins to look in for PM and DL
% PM first
is_PM_time = zeros(size(bin_times));
for i = 1:length(times_PM)
    is_PM_time = is_PM_time | (bin_times>times_PM(i,1) & bin_times<times_PM(i,2));
end
% then DL
is_DL_time = zeros(size(bin_times));
for i = 1:length(times_DL)
    is_DL_time = is_DL_time | (bin_times>times_DL(i,1) & bin_times<times_DL(i,2));
end

% clear unneeded things
clear t
clear pos
clear vel
clear spd

%% Fit neurons to endpoint
bootfunc = @(X,y) GeneralizedLinearModel.fit(X,y,'Distribution','poisson');

tic
for unit_ctr = 1:length(bdf.units)
    % skip unsorted or invalid units
    if bdf.units(unit_ctr).id(2)==0 || bdf.units(unit_ctr).id(2)==255
        continue;
    end
    
    % get unit id
    uid = bdf.units(unit_ctr).id;
    
    % get firing rate
    FR = bdf.units(unit_ctr).FR(:,2)/0.05;
    bin_times = bdf.units(unit_ctr).FR(:,1);
    FR(bin_times<time_window(1) | bin_times>time_window(2)) = [];
    
    endpoint_tuning = calc_PD_helper(bootfunc,endpoint_kin_full,FR,['Calculating endpoint chan ' num2str(uid(1)) ', unit ' num2str(uid(2)) ' (Time: ' num2str(toc) ')']);
    endpoint_curve = get_single_tuning_curve(endpoint_kin_full(:,3:4),FR);
    
    endpoint_tuning_PM = calc_PD_helper(bootfunc,endpoint_kin_full(is_PM_time,:),FR(is_PM_time),'PM tuning');
    endpoint_curve_PM = get_single_tuning_curve(endpoint_kin_full(is_PM_time,3:4),FR(is_PM_time));
    
    endpoint_tuning_DL = calc_PD_helper(bootfunc,endpoint_kin_full(is_DL_time,:),FR(is_DL_time),'DL tuning');
    endpoint_curve_DL = get_single_tuning_curve(endpoint_kin_full(is_DL_time,3:4),FR(is_DL_time));
    
    % set up figure
    handle = figure('name',['chan ' num2str(uid(1)) ', unit ' num2str(uid(2)) 'tuning curve']);
    
    % find maximum curve for plot
    max_curve = max([endpoint_curve.CI_high;endpoint_curve_PM.CI_high;endpoint_curve_DL.CI_high]);
    
    % plot all three tunings
    subplot(2,2,1)
    plot_tuning(endpoint_tuning,endpoint_curve,max_curve,[0 0 0])
    subplot(2,2,2)
    plot_tuning(endpoint_tuning_PM,endpoint_curve_PM,max_curve,[0.6 0.5 0.7]);
    subplot(2,2,3)
    plot_tuning(endpoint_tuning_DL,endpoint_curve_DL,max_curve,[1 0 0]);
    subplot(2,2,4)
    plot_tuning(endpoint_tuning,endpoint_curve,max_curve,[0 0 0])
    plot_tuning(endpoint_tuning_PM,endpoint_curve_PM,max_curve,[0.6 0.5 0.7]);
    plot_tuning(endpoint_tuning_DL,endpoint_curve_DL,max_curve,[1 0 0]);
end

%% Fit neuron to muscle in one workspace (bootstrap)

%% Predict tuning in other workspace

%% Plot tuning in both workspaces