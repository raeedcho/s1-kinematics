function model_diff_table = findActpasModelDiff
% This function looks for a difference in parameters between the active and passive handelbow models
%% Set up meta info
% model_aliases = {'ext','extforce','handelbow','ext_actpasbaseline'};
% model_aliases = {'ext','extforce','handelbow'};
model_type = 'glm';
arrayname = 'S1';
num_boots = 1000;
savefile = true;

if ispc
    dataroot = '';
else
    dataroot = '/data/raeed';
end

datadir = fullfile(dataroot,'project-data','limblab','s1-kinematics','td-library');
file_info = dir(fullfile(datadir,'*COactpas*.mat'));
filenames = horzcat({file_info.name})';
if savefile
    savedir = fullfile(dataroot,'project-data','limblab','s1-kinematics','Results','Separability');
    run_date = char(datetime('today','format','yyyyMMdd'));
    savename = sprintf('handelbow_actpasModelDiff_run%s.mat',run_date);
end

%% Loop through files
model_diff = cell(4,1);
filetic = tic;
fprintf('Starting analysis...\n')
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
    
    % for Chips_20170912, split bumps out of trials (different trial structure)
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

    % even out sizes and put back together
    minsize = min(length(td_act),length(td_pas));
    td_act = td_act(1:minsize);
    td_pas = td_pas(1:minsize);
    td_bin = cat(2,td_act,td_pas);

    % remove low firing neurons
    td_bin = removeBadNeurons(td_bin,struct(...
        'min_fr',1,...
        'fr_window',{{'idx_movement_on',0;'idx_movement_on',11}},...
        'calc_fr',true));
    
    % add firing rates in addition to spike counts
    td_bin = addFiringRates(td_bin,struct('array',arrayname));

    % shift S1 backwards by 50 ms to account for lag between kinematics and S1
    td_bin = dupeAndShift(td_bin,[arrayname '_FR'],5);

    % trim to just movements
    td_bin = trimTD(td_bin,{'idx_movement_on',0},{'idx_movement_on',11});

    % find average over the movement
    td_bin = binTD(td_bin,'average');
    
    % split into active and passive again
    [~,td_act] = getTDidx(td_bin,'ctrHoldBump',false);
    [~,td_pas] = getTDidx(td_bin,'ctrHoldBump',true);

    %% find separabilities
    % suppress getTDfields warning...
    getTDfields(td_bin,'time');
    onetime_warn = warning('query','last'); 
    warning('off',onetime_warn.identifier)
    
    %% fit models and find parameter diffs
    % set up glm_parameters
    % indices for cartesian hand coordinates
    markername = 'Marker_1';
    [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_bin(1).marker_names);
    assert(all(point_exists),'Hand marker does not exist?')

    markername = 'Pronation_Pt1';
    [point_exists,marker_elbow_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_bin(1).marker_names);
    assert(all(point_exists),'Elbow marker does not exist?')

    glm_params = struct(...
        'model_type',model_type,...
        'model_name','handelbow_model',...
        'in_signals',{{'markers',[marker_hand_idx marker_elbow_idx];'marker_vel',[marker_hand_idx marker_elbow_idx]}},...
        'out_signals',[arrayname '_FR']);

    % bootstrap
    if num_boots>1
        model_diff{filenum} = cell(1,num_boots);
        boottic = tic;
        for bootnum = 1:num_boots
            % get bootstrap sample indices (rand sample with replacement)
            bootsamp = randsample(length(td_act),length(td_act),true);
            [~,glm_info_act] = getModel(td_act(bootsamp),glm_params);
            [~,glm_info_pas] = getModel(td_pas(bootsamp),glm_params);

            % parameters are n_params x n_neurons
            param_diff = glm_info_pas.b - glm_info_act.b;

            % make table of parameter diffs
            model_diff{filenum}{bootnum} = horzcat(...
                makeNeuronTableStarter(td_act,struct(...
                    'out_signal_names',td_act(1).([arrayname '_unit_guide']),...
                    'meta',struct('bootID',bootnum))),...
                table(param_diff','VariableNames',{'handelbow_param_diff'}));

            if mod(bootnum,floor(num_boots/10))==0
                fprintf('\tCompleted bootstrap %d at time %f\n\r',bootnum,toc(boottic));
            end
        end
        model_diff{filenum} = vertcat(model_diff{filenum}{:});
    else
        [~,glm_info_act] = getModel(td_act,glm_params);
        [~,glm_info_pas] = getModel(td_pas,glm_params);

        % parameters are n_params x n_neurons
        param_diff = glm_info_pas.b - glm_info_act.b;

        % make table of parameter diffs
        model_diff{filenum} = horzcat(...
            makeNeuronTableStarter(td_act,struct('out_signal_names',td_act(1).([arrayname '_unit_guide']))),...
            table(param_diff','VariableNames',{'handelbow_param_diff'}));
    end

    % turn warning back on
    warning('on',onetime_warn.identifier)

    % diagnostic
    fprintf('Completed file %d at time %f\n',filenum,toc(filetic))
end

% save
model_diff_table = vertcat(model_diff{:});
if savefile
    save(fullfile(savedir,savename),'model_diff_table')
end
end
