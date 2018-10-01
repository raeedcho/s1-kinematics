datadir = '/home/raeed/Wiki/Projects/limblab/multiworkspace/data';
fileprefix = {'Han_20171101_TRT','Chips_20170915_TRT','Lando_20170802_RWTW'};
savesuffix = '_encodingResults_9MuscPCs_run20180924.mat';

model_aliases = {'ext','ego','joint','musc','handelbow','ego_handelbow'};
arrayname = 'S1';
num_musc_pcs = 9;

for filenum = 1:length(fileprefix)
    clear encoderResults

    %% Load data
    load(fullfile(datadir,[fileprefix{filenum} '_TD.mat']))

    % prep trial data by getting only rewards and trimming to only movements
    % first process marker data
    td = trial_data;
    [~,td] = getTDidx(td,'result','R');
    td = smoothSignals(td,struct('signals','markers'));
    td = getDifferential(td,struct('signals','markers','alias','marker_vel'));
    % add firing rates rather than spike counts
    td = addFiringRates(td,struct('array',arrayname));

    % for active movements
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    % for bumps
    % td = td(~isnan(cat(1,td.idx_bumpTime)));
    % td = trimTD(td,{'idx_bumpTime',0},{'idx_bumpTime',15});

    % bin data at 50ms
    td = binTD(td,5);

    %% Get encoding models
    encoderResults = mwEncoders(td,struct('model_aliases',{model_aliases},'arrayname',arrayname,'num_tuning_bins',16,'num_repeats',20,'num_folds',5));

    %% save
    save(fullfile(datadir,'Results','Encoding',[fileprefix{filenum} savesuffix]),'encoderResults')
end
