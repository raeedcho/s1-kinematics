%%% Make simulation of two versus one workspace (partitioned later in analysis)
clear

%% Define workspaces (assume right handed)
full_limits = [-10 10 -10 10]; %[xmin xmax ymin ymax]
DL_limits = [0 10 0 10];
PM_limits = [-10 0 -10 0];
spacing = [0 inf];
target_size = 0.5;
hold_time = 0.05;
num_targets = 10000;

%% Generate targets for each workspace
[full_targets,full_movements] = randomwalk_sim(full_limits,num_targets,target_size,spacing,hold_time,0.01);
[DL_targets, DL_movements] = randomwalk_sim(DL_limits,uint16(num_targets/10),target_size,spacing,hold_time,0.01);
[PM_targets, PM_movements] = randomwalk_sim(PM_limits,uint16(num_targets/10),target_size,spacing,hold_time,0.01);

%% Find partitioned targets
DL_part_target_idx = target_in_workspace(full_targets,DL_limits);
PM_part_target_idx = target_in_workspace(full_targets,PM_limits);

%% Analyze target directions
full_vectors = diff([0 0;full_targets]);
full_directions = atan2(full_vectors(:,2),full_vectors(:,1));

DL_vectors = diff([0 0;DL_targets]);
DL_directions = atan2(DL_vectors(:,2),DL_vectors(:,1));

PM_vectors = diff([0 0;PM_targets]);
PM_directions = atan2(PM_vectors(:,2),PM_vectors(:,1));

DL_part_directions = full_directions(DL_part_target_idx);
PM_part_directions = full_directions(PM_part_target_idx);

%% Find Movement directions
full_movement_dir = atan2(full_movements(:,5),full_movements(:,4));
DL_movement_dir = atan2(DL_movements(:,5),DL_movements(:,4));
PM_movement_dir = atan2(PM_movements(:,5),PM_movements(:,4));

DL_movement_idx = target_in_workspace(full_movements(:,2:3),DL_limits);
PM_movement_idx = target_in_workspace(full_movements(:,2:3),PM_limits);
DL_part_movement_dir = atan2(full_movements(DL_movement_idx,5),full_movements(DL_movement_idx,4));
PM_part_movement_dir = atan2(full_movements(PM_movement_idx,5),full_movements(PM_movement_idx,4));

%% Plot target locations
figure
plot(full_targets(:,1),full_targets(:,2),'ro','Markersize',15)
hold on
plot(full_movements(:,2),full_movements(:,3))
axis(full_limits)

%% Plot direction distribution
figure
subplot(221)
rose(full_directions)
title 'Full Target Directions'
subplot(222)
rose(full_movement_dir)
title 'Full Movement Directions'

figure
subplot(221)
rose(DL_directions)
title 'DL Target Directions'
subplot(223)
rose(DL_part_directions)
title 'Partitioned DL Target Directions'
subplot(222)
rose(DL_movement_dir)
title 'DL Movement Directions'
subplot(224)
rose(DL_part_movement_dir)
title 'Partitioned DL Movement Directions'

figure
subplot(221)
rose(PM_directions)
title 'PM Target Directions'
subplot(223)
rose(PM_part_directions)
title 'Partitioned PM Target Directions'
subplot(222)
rose(PM_movement_dir)
title 'PM Movement Directions'
subplot(224)
rose(PM_part_movement_dir)
title 'Partitioned PM Movement Directions'