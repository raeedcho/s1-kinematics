datadir = '/home/raeed/Projects/limblab/data-td/FullWS';
fileprefix = {'Chips_20151211_RW'};
savesuffix = '_decodingResults_run20180813.mat';

for filenum = 1:length(fileprefix)
    clear decoderResults

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
    td = getDifferential(td,struct('signal','markers','alias','marker_vel'));

    % bin data at 50ms
    td = binTD(td,5);

    % duplicate and shift spike data before chopping data up
    % Get future 150 ms of S1 activity to predict current kinematics
    % non-overlapping bins...
    neur_name = 'S1_spikes';
    td = dupeAndShift(td,neur_name,-3);

    % Trim td
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    %% Get encoding models
    decoderResults = mwDecoders(td,struct('neur_name','S1_spikes','add_shift',true,'num_repeats',10,'num_folds',10));

    %% save
    save(fullfile(datadir,'Results','Decoding',[fileprefix{filenum} savesuffix]),'decoderResults')
end

