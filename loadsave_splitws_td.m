%% Create TD
params.array_alias = {'LeftS1Area2','S1'};
params.exclude_units = [0,255];
params.trial_results = {'R','A','F','I'};
meta = struct('task','RW','epoch','DL');
params.meta = meta;

trial_data_DL = parseFileByTrial(cds_DL,params);

clear params
params.array_alias = {'LeftS1Area2','S1'};
params.exclude_units = [0,255];
params.trial_results = {'R','A','F','I'};
meta = struct('task','RW','epoch','PM');
params.meta = meta;

trial_data_PM = parseFileByTrial(cds_PM,params);

trial_data = [trial_data_DL trial_data_PM];

%% Save TD

save('C:\Users\rhc307\Projects\limblab\data-preproc\MultiWorkspace\SplitWS\Han\20160322\TD\Han_20160322_SplitWS_area2_TD.mat','trial_data','-v7.3')
