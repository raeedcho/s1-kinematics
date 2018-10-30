% datadir = '/home/raeed/data/limblab/data-td/FullWS';
% fileprefix = {'Chips_20151211_RW'};
datadir = '/home/raeed/data/limblab/data-td/ActPas';
fileprefix = {'Han_20170203_COactpas','Chips_20170913_COactpas'};
savesuffix = '_decodingResults_run20181029.mat';

for filenum = 1:length(fileprefix)
    clear decoderResults

    %% Load data
    load(fullfile(datadir,[fileprefix{filenum} '_TD.mat']))

    % prep trial data by getting only rewards and trimming to only movements
    % first process marker data
    td = trial_data;
    [~,td] = getTDidx(td,'result','R');
    % Remove unsorted channels
    keepers = (td(1).S1_unit_guide(:,2)~=0);
    for trial = 1:length(td)
        td(trial).S1_unit_guide = td(trial).S1_unit_guide(keepers,:);
        td(trial).S1_spikes = td(trial).S1_spikes(:,keepers);
    end

    % remove low firing neurons
    td = removeBadNeurons(td,struct('min_fr',0.1));
    
    td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
    % remove trials with no go cue
    % if any(isnan([td.idx_targetStartTime]))
    %     warning('Some trials have no go cue time, deleting trials...')
    %     td(isnan([td.idx_targetStartTime])) = [];
    % end

    td = smoothSignals(td,struct('signals','markers'));
    td = getDifferential(td,struct('signals','markers','alias','marker_vel'));

    % bin data at 50ms
    td = binTD(td,5);

    % duplicate and shift spike data before chopping data up
    % Get future 150 ms of S1 activity to predict current kinematics
    % non-overlapping bins...
    neur_name = 'S1_spikes';
    td = dupeAndShift(td,neur_name,-3);

    % Trim td
    % td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});
    % for actpas
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    % clean nans out...?
    nanners = isnan(cat(1,td_act.target_direction));
    td_act = td_act(~nanners);
    td_act = trimTD(td_act,{'idx_movement_on',-3},{'idx_movement_on',5});

    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    td_pas = trimTD(td_pas,{'idx_bumpTime',-3},{'idx_bumpTime',5});
    td = cat(2,td_act,td_pas);

    %% Get encoding models
    decoderResults = mwDecoders(td,struct('neur_name','S1_spikes','add_shift',true,'num_repeats',10,'num_folds',10));

    %% save
    save(fullfile(datadir,'Results','Decoding',[fileprefix{filenum} savesuffix]),'decoderResults')
end

