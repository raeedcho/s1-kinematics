%% play muscle tuning

%% Read files in
folder = '/home/raeed/Projects/limblab/FSMRes/limblab/User_folders/Raeed/Arm Model/Data/Chips/experiment_20151203_RW_002/';
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
joint_pos = array2table(joint_pos_mat,'VariableNames',{'time','shoulder_adduction','shoulder_rotation','shoulder_flexion','elbow_flexion','radial_pronation','wrist_flexion','wrist_abduction'});
clear joint_pos_mat

muscle_pos_mat = csvread([folder 'Analysis/' opensim_prefix '_MuscleAnalysis_Length.sto'],12,0);
muscle_pos = array2table(muscle_pos_mat,'VariableNames',{'time','abd_poll_longus','anconeus','bicep_lh','bicep_sh','brachialis','brachioradialis','coracobrachialis','deltoid_ant','deltoid_med','deltoid_pos','dorsoepitrochlearis','ext_carpi_rad_longus','ext_carp_rad_brevis','ext_carpi_ulnaris','ext_digitorum','ext_digiti','ext_indicis','flex_carpi_radialis','flex_carpi_ulnaris','flex_digit_profundus','flex_digit_superficialis','flex_poll_longus','infraspinatus','lat_dorsi_sup','lat_dorsi_cen','lat_dorsi_inf','palmaris_longus','pectoralis_sup','pectoralis_inf','pronator_quad','pronator_teres','subscapularis','supinator','supraspinatus','teres_major','teres_minor','tricep_lat','tricep_lon','tricep_sho'});
clear muscle_pos_mat

%% get muscle PCA
[coeff,score,latent] = pca(muscle_pos{:,2:end});

% use first 5 PCs
muscle_scores = score(:,1:5);

%% get velocities
% joint velocities
joint_vel = joint_pos;
for i=2:size(joint_pos,2)
    joint_vel{:,i} = gradient(joint_pos{:,i},joint_pos.time);
end

% muscle velocities
muscle_vel = muscle_pos;
for i=2:size(muscle_pos,2)
    muscle_vel{:,i} = gradient(muscle_pos{:,i},muscle_pos.time);
end

% muscle PC velocities
muscle_PC_vel_mat = muscle_scores;
for i=1:size(muscle_scores,2)
    muscle_PC_vel_mat(:,i) = gradient(muscle_scores(:,i),muscle_pos.time);
end

joint_vel_mat = joint_vel{:,2:end};
muscle_vel_mat = muscle_vel{:,2:end};
endpoint_vel_mat = bdf.vel(:,2:3);

clear i

%% interpolate intrinsic velocities to 50 ms time points
% get time vector
time_vect = bdf.units(1).FR(:,1);

% cut out times outside relevant window
if( exist('time_window','var') && length(time_window(:))==2 )
    time_vect(time_vect<time_window(1) | time_vect>time_window(2)) = [];
else
    warning('No time window defined')
end

joint_vel_interp = interp1(joint_vel.time,joint_vel_mat,time_vect);
muscle_vel_interp = interp1(muscle_vel.time,muscle_vel_mat,time_vect);
muscle_PC_vel_interp = interp1(muscle_vel.time,muscle_PC_vel_mat,time_vect);
endpoint_vel_interp = interp1(bdf.vel(:,1),endpoint_vel_mat,time_vect);

%% plot each muscle PC against endpoint
figure
for i = 1:5
    subplot(2,5,i)
    scatter(muscle_PC_vel_interp(:,i)/var(muscle_PC_vel_interp(:,i)),endpoint_vel_interp(:,1)/var(endpoint_vel_interp(:,1)))
    subplot(2,5,i+5)
    scatter(muscle_PC_vel_interp(:,i)/var(muscle_PC_vel_interp(:,i)),endpoint_vel_interp(:,2)/var(endpoint_vel_interp(:,2)))
end

%% find correlations
R = corrcoef([muscle_PC_vel_interp endpoint_vel_interp]);
imagesc(R)
colormap jet
colorbar