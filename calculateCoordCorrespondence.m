datadir = '/home/raeed/codebase/limblab/data-td/MultiWorkspace';
% fileprefix = {'Han_20160325_RWhold','Chips_20151211_RW'};
fileprefix = {'Han_20171101_TRT','Chips_20170915_TRT','Lando_20170802_RWTW'};
savesuffix = '_coordResults_run20180813.mat';

for filenum = 1:length(fileprefix)
    clear coordResults

    %% Load data
    load(fullfile(datadir,[fileprefix{filenum} '_TD.mat']))

    % prep trial data by getting only rewards and trimming to only movements
    % first process marker data
    td = trial_data;
    [~,td] = getTDidx(td,'result','R');
    % remove trials with no go cue
    if any(isnan([td.idx_targetStartTime]))
        warning('Some trials have no go cue time, deleting trials...')
        td(isnan([td.idx_targetStartTime])) = [];
    end

    td = smoothSignals(td,struct('signals','markers'));
    td = getDifferential(td,struct('signals','markers','alias','marker_vel'));

    % bin data at 50ms
    td = binTD(td,5);

    % Trim td
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    %% Get VAF for model of muscle length from ext or hand/elbow
    coordResults = coordCorr(td,struct('model_aliases',{{'ext','handelbow','joint','musc'}},'num_musc_pcs',5,'num_repeats',10,'num_folds',10));

    %% save
    save(fullfile(datadir,'Results','CoordCorr',[fileprefix{filenum} savesuffix]),'coordResults')
end

