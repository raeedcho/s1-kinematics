%% Set up meta info
    if ispc
        dataroot = 'G:\raeed\';
    else
        dataroot = '/data/raeed/';
    end
    
    datadir = fullfile(dataroot,'project-data','limblab','s1-kinematics','td-library');
    file_info = dir(fullfile(datadir,'*TRT*'));
    filenames = horzcat({file_info.name})';
    savedir = fullfile(dataroot,'project-data','limblab','s1-kinematics','td-library','lib-s1-kin-paper');
    arrayname = 'S1';
    
%% Loop through files
    trial_data_cell = cell(1,length(filenames));
    for filenum = 1:length(filenames)
        %% Load data
        td = load(fullfile(datadir,[filenames{filenum}]));
    
        % rename trial_data for ease
        td = td.trial_data;

        % anonymize trial data
        for trialnum = 1:length(td)
            td(trialnum).monkey = td(trialnum).monkey(1);
        end
    
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
        
        % remove trials where muscles aren't present
        bad_trial = false(length(td),1);
        for trialnum = 1:length(td)
            if any(any(isnan(td(trialnum).muscle_len) | isnan(td(trialnum).muscle_vel)))
                bad_trial(trialnum) = true;
            end
        end
        td(bad_trial) = [];
        fprintf('Removed %d trials because of missing muscles\n',sum(bad_trial))

        trial_data_cell{filenum} = td;
    end

%% save
save(fullfile(savedir,'twoworkspace_trim_TD.mat'),'trial_data_cell','-v7.3')
