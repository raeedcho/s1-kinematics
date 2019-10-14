% This script will compare the target decoding capabilities of neural activity and kinematics over time
% for active and passive trials

%% Set up
    if ispc
        homefolder = 'C:\Users\rhc307';
        dataroot = '';
    else
        homefolder = '/home/raeed';
        dataroot = '/data/raeed';
    end
    
    datadir = fullfile(dataroot,'project-data','limblab','s1-kinematics','td-library');
    file_info = dir(fullfile(datadir,'*COactpas*.mat'));
    filenames = horzcat({file_info.name})';
    arrayname = 'S1';
    crossval_lookup = [];

%% Loop through files
    for filenum = 1:4%length(filenames)
        %% load and preprocess data
        td = load(fullfile(datadir,[filenames{filenum}]));
    
        % rename trial_data for ease
        td = td.trial_data;
    
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
        
        % get speed and ds
        td = getNorm(td,struct('signals','vel','field_extra','_norm'));
        td = getDifferential(td,struct('signals','vel_norm','alias','dvel_norm'));
        
        % get norm planar force and derivative
        td = getNorm(td,struct('signals',{{'force',1:2}},'field_extra','_plane_norm'));
        td = getDifferential(td,struct('signals','force_plane_norm','alias','dforce_plane_norm'));
        
        % remove unsorted neurons
        unit_ids = td(1).([arrayname '_unit_guide']);
        unsorted_units = (unit_ids(:,2)==0);
        new_unit_guide = unit_ids(~unsorted_units,:);
        
        for trialnum = 1:length(td)
            td(trialnum).(sprintf('%s_unit_guide',arrayname)) = new_unit_guide;
            
            spikes = td(trialnum).(sprintf('%s_spikes',arrayname));
            spikes(:,unsorted_units) = [];
            td(trialnum).(sprintf('%s_spikes',arrayname)) = spikes;
        end
    
        % prep trial data by getting only rewards and trimming to only movements
        % split into trials
        td = splitTD(...
            td,...
            struct(...
                'split_idx_name','idx_startTime',...
                'linked_fields',{{...
                    'trialID',...
                    'result',...
                    'bumpDir',...
                    'tgtDir',...
                    'ctrHoldBump',...
                    'ctrHold',...
                    }},...
                'start_name','idx_startTime',...
                'end_name','idx_endTime'));
        [~,td] = getTDidx(td,'result','R');
        td = reorderTDfields(td);
        
        % clean nans out...?
        nanners = isnan(cat(1,td.tgtDir));
        td = td(~nanners);
        fprintf('Removed %d trials because of missing target direction\n',sum(nanners))
        biggers = cat(1,td.ctrHoldBump) & abs(cat(1,td.bumpDir))>360;
        td = td(~biggers);
        fprintf('Removed %d trials because bump direction makes no sense\n',sum(biggers))
    
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
        
        % for Chips_20170912, file has bumps on 100% of trials
        if strcmpi(td(1).monkey,'Chips') && contains(td(1).date_time,'2017/9/12')
            td_copy = td;
            [td_copy.ctrHoldBump] = deal(false);
            td = cat(2,td,td_copy);
            clear td_copy
        end
        
        % split into active and passive
        [~,td_act] = getTDidx(td,'ctrHoldBump',false);
        [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    
        % find the relevant movmement onsets
        td_act = getMoveOnsetAndPeak(td_act,struct(...
            'start_idx','idx_goCueTime',...
            'start_idx_offset',20,...
            'peak_idx_offset',20,...
            'end_idx','idx_endTime',...
            'method','peak',...
            'peak_divisor',10,...
            'min_ds',1));
        td_pas = getMoveOnsetAndPeak(td_pas,struct(...
            'start_idx','idx_bumpTime',...
            'start_idx_offset',-5,... % give it some wiggle room
            'peak_idx_offset',-5,... % give it some wiggle room
            'end_idx','idx_goCueTime',...
            'method','peak',...
            'peak_divisor',10,...
            'min_ds',1));
        % throw out all trials where bumpTime and movement_on are more than 3 bins apart
        bad_trial = isnan(cat(1,td_pas.idx_movement_on)) | abs(cat(1,td_pas.idx_movement_on)-cat(1,td_pas.idx_bumpTime))>3;
        td_pas = td_pas(~bad_trial);
        fprintf('Removed %d trials because of bad movement onset\n',sum(bad_trial))
    
        % Get 50 ms bins of movement
        td_act = trimTD(td_act,{'idx_movement_on',-15},{'idx_movement_on',20});
        td_act = binTD(td_act,5);
        td_pas = trimTD(td_pas,{'idx_movement_on',-15},{'idx_movement_on',20});
        td_pas = binTD(td_pas,5);

        % even out sizes and put back together
        minsize = min(length(td_act),length(td_pas));
        td_act = td_act(1:minsize);
        td_pas = td_pas(1:minsize);
        td = cat(2,td_act,td_pas);
    
        % remove low firing neurons
        td = removeBadNeurons(td,struct(...
            'min_fr',1,...
            'calc_fr',true));
        
        % add firing rates in addition to spike counts
        td = addFiringRates(td,struct('array',arrayname));

        % train models on active movements, test on passive movements
        act_idx = getTDidx(td,'ctrHoldBump',false);
        pas_idx = getTDidx(td,'ctrHoldBump',true);

        % crossvalidate decoders for active and passive
        repeat_tic = tic;
        fprintf('Starting %dx%d cross-validation at time: %f\n',num_repeats,num_folds,toc(repeat_tic))
        meta_table = cell2table({td_bin(1).monkey,td_bin(1).date_time,td_bin(1).task},...
            'VariableNames',{'monkey','date','task'});
        meta_table.Properties.VariableDescriptions = repmat({'meta'},1,3);
        [trial_table_cell,lda_table_cell,foldEval] = deal(cell(num_repeats,num_folds));
        for repeatnum = 1:num_repeats
            foldidx = crossvalind('kfold',length(td_act),num_folds);
            fold_tic = tic;
            for foldnum = 1:num_folds
                % get crossval table entry
                crossval_table = table(uint16([repeatnum foldnum]),'VariableNames',{'crossvalID'});
                crossval_table.Properties.VariableDescriptions = {'meta'};

                % split into training and testing
                if isempty(crossval_lookup)
                    train_idx = (foldidx~=foldnum);
                    td_train = cat(2,td_act(train_idx),td_pas(train_idx));
                    test_idx = (foldidx==foldnum);
                    td_test = cat(2,td_act(test_idx),td_pas(test_idx));
                else
                    % read the crossval_lookup
                    [~,current_crossval] = getNTidx(crossval_lookup,'crossvalID',[repeatnum foldnum]);
                    td_whole = cat(2,td_act,td_pas);
                    whole_ids = cat(1,td_whole.trialID);
                    whole_isPassive = cat(1,false(length(td_act),1),true(length(td_pas),1));
                    test_idx_whole = ismember([whole_ids whole_isPassive],current_crossval{:,{'trialID' 'isPassive'}},'rows');
                    train_idx_whole = ~test_idx_whole;
                    
                    td_test = td_whole(test_idx_whole);
                    td_train = td_whole(train_idx_whole);
                end

                % train and test models
                %...

                fprintf('\tEvaluated fold %d at time: %f\n',foldnum,toc(fold_tic))
            end
            fprintf('Evaluated repeat %d at time: %f\n',repeatnum,toc(repeat_tic))
        end
    end

