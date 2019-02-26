function checkActPasKin
    %% Set up meta info
    if ispc
        homefolder = 'C:\Users\rhc307';
    else
        homefolder = '/home/raeed';
    end

    datadir = fullfile(homefolder,'data','project-data','limblab','s1-kinematics','td-library');
    file_info = dir(fullfile(datadir,'*COactpas*'));
    filenames = horzcat({file_info.name})';
    
    for filenum = 1:length(filenames)
        makeKinFigs(fullfile(datadir,filenames{filenum}))
    end
    
end

function makeKinFigs(filename)
    %% load and preprocess data
    td = load(filename);

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
    td = getNorm(td,struct('signals','vel','norm_name','speed'));
    td = getDifferential(td,struct('signals','speed','alias','dspeed'));

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
    bad_trials = isnan(cat(1,td_pas.idx_movement_on)) | abs(cat(1,td_pas.idx_movement_on)-cat(1,td_pas.idx_bumpTime))>3;
    td_pas = td_pas(~bad_trials);

    % even out sizes and put back together
    minsize = min(length(td_act),length(td_pas));
    td_act = td_act(1:minsize);
    td_pas = td_pas(1:minsize);
    td_bin = cat(2,td_act,td_pas);

    % trim to just movements
    td_bin = trimTD(td_bin,{'idx_movement_on',0},{'idx_movement_on',14});

    % find average over the movement
    td_bin = binTD(td_bin,'average');
    
    % make figure
    [~,td_act_bin] = getTDidx(td_bin,'ctrHoldBump',false);
    [~,td_pas_bin] = getTDidx(td_bin,'ctrHoldBump',true);

    [~,~,dir_idx_act] = unique(cat(1,td_act_bin.tgtDir));
    [~,~,dir_idx_pas] = unique(cat(1,td_pas_bin.bumpDir));
    dir_colors = linspecer(max(dir_idx_act));
    forcename = {'force',1:3};
    kinname = {'marker_vel',7:9};

    figure
    subplot(1,2,1)
    sig_act = getSig(td_act_bin,kinname);
    sig_pas = getSig(td_pas_bin,kinname);
    scatter3(sig_act(:,1),sig_act(:,2),sig_act(:,3),[],dir_colors(dir_idx_act,:),'filled')
    hold on
    scatter3(sig_pas(:,1),sig_pas(:,2),sig_pas(:,3),[],dir_colors(dir_idx_pas,:))
    axis equal
    subplot(1,2,2)
    sig_act = getSig(td_act_bin,forcename);
    sig_pas = getSig(td_pas_bin,forcename);
    scatter3(sig_act(:,1),sig_act(:,2),sig_act(:,3),[],dir_colors(dir_idx_act,:),'filled')
    hold on
    scatter3(sig_pas(:,1),sig_pas(:,2),sig_pas(:,3),[],dir_colors(dir_idx_pas,:))
    axis equal
    [~,figname] = fileparts(filename);
    suptitle(strrep(figname,'_','-'))
end