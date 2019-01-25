%% Set up meta info
if ispc
    homefolder = 'C:\Users\rhc307';
else
    homefolder = '/home/raeed';
end

datadir = fullfile(homefolder,'data','cds-library');
fileprefix = {'Chips_20170907'};
filename = strcat(fileprefix,'_CDS_TRT_001.mat');
savedir = fullfile(homefolder,'data','project-data','limblab','s1-kinematics','Results','Encoding');
savesuffix = '_encodingResults_allModels_run20190125.mat';

model_aliases = {'ext','ego','joint','musc','handelbow','ego_handelbow'};
num_musc_pcs = 5;

% for TRT files...
load_params = struct(...
    'array_name','S1',...
    'cont_signal_names',{{...
        'pos',...
	    'vel',...
	    'markers',...
	    'joint_ang',...
	    'joint_vel',...
	    'muscle_len',...
	    'muscle_vel',...
        }},...
    'event_names',{{...
        'startTime',...
        'endTime',...
        'goCueTime',...
        'bumpTime',...
        'targetStartTime',...
        'ctHoldTime',...
        'otHoldTime',...
        }},...
    'trial_meta',{{...
        'bumpDir',...
        'spaceNum',...
        }});

%% Loop through files
for filenum = 1:length(fileprefix)
    clear encoderResults

    %% Load data
    td = loadTDfromCDS(fullfile(datadir,filename{filenum}),load_params);

    % to save memory, get rid of some continuous signals we don't need
    cont_fields = getTDfields(td,'cont');
    unused_fields = cont_fields(~ismember(cont_fields,required_signals));
    td = rmfield(td,unused_fields);

    % first process marker data
    % find times when markers are NaN and replace with zeros temporarily
    for trialnum = 1:length(td)
        markernans = isnan(td(trialnum).markers);
        td(trialnum).markers(markernans) = 0;
        td(trialnum) = smoothSignals(td(trialnum),struct('signals','markers'));
        td(trialnum).markers(markernans) = NaN;
        clear markernans
    end

    % get marker velocity
    td = getDifferential(td,struct('signals','markers','alias','marker_vel'));
    
    % remove unsorted neurons
    unit_ids = td(1).S1_unit_guide;
    unsorted_units = (unit_ids(:,2)==0);
    new_unit_guide = unit_ids(~unsorted_units,:);
    
    for trialnum = 1:length(td)
        td(trialnum).(sprintf('%s_unit_guide',arrayname)) = new_unit_guide;
        
        spikes = td(trialnum).(sprintf('%s_spikes',arrayname));
        spikes(:,unsorted_units) = [];
        td(trialnum).(sprintf('%s_spikes',arrayname)) = spikes;
    end
    
    % add firing rates in addition to spike counts
    td = addFiringRates(td,struct('array',arrayname));

    % prep trial data by getting only rewards and trimming to only movements
    % split into trials
    td = splitTD(...
        td,...
        struct(...
            'split_idx_name','idx_startTime',...
            'linked_fields',{{...
                'trialID',...
                'result',...
                'spaceNum',...
                'bumpDir',...
                }},...
            'start_name','idx_startTime',...
            'end_name','idx_endTime'));
    [~,td] = getTDidx(td,'result','R');
    td = reorderTDfields(td);

    % for active movements
    % remove trials without a target start (for whatever reason)
    td(isnan(cat(1,td.idx_targetStartTime))) = [];
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    % remove trials where markers aren't present
    bad_trial = false(length(td),1);
    for trialnum = 1:length(td)
        if any(any(isnan(td(trialnum).markers)))
            bad_trial(trialnum) = true;
        end
    end
    td(bad_trial) = [];
    fprintf('Removed %d trials because of missing markers\n',sum(bad_trial))

    % for bumps
    % td = td(~isnan(cat(1,td.idx_bumpTime)));
    % td = trimTD(td,{'idx_bumpTime',0},{'idx_bumpTime',15});

    % bin data at 50ms
    td = binTD(td,50);

    %% Get encoding models
    encoderResults = mwEncoders(td,struct('model_aliases',{model_aliases},'arrayname',arrayname,'num_tuning_bins',16,'num_repeats',20,'num_folds',5));

    %% save
    save(fullfile(savedir,[fileprefix{filenum} savesuffix]),'encoderResults')
end
