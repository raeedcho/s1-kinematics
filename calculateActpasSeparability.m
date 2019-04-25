%% Set up meta info
if ispc
    homefolder = 'C:\Users\rhc307';
else
    homefolder = '/home/raeed';
end

datadir = fullfile(homefolder,'data','project-data','limblab','s1-kinematics','td-library');
file_info = dir(fullfile(datadir,'*COactpas*'));
filenames = horzcat({file_info.name})';
savedir = fullfile(homefolder,'data','project-data','limblab','s1-kinematics','Results','Separability');
savesuffix = '_separationResults_run20190425.mat';

model_aliases = {'ext','extforce','joint','musc','handelbow'};
model_type = 'glm';
arrayname = 'S1';
num_musc_pcs = 5;
num_pcs = 3;
num_repeats = 20;
num_folds = 5;

%% Loop through files
for filenum = 2%1:4%length(filenames)
    clear sepResults

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
    td = getNorm(td,struct('signals','vel','norm_name','speed'));
    td = getDifferential(td,struct('signals','speed','alias','dspeed'));
    
    % get norm planar force and derivative
    td = getNorm(td,struct('signals',{{'force',1:2}},'norm_name','plane_force_norm'));
    td = getDifferential(td,struct('signals','plane_force_norm','alias','dplane_force_norm'));
    
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

    % even out sizes and put back together
    minsize = min(length(td_act),length(td_pas));
    td_act = td_act(1:minsize);
    td_pas = td_pas(1:minsize);
    td_bin = cat(2,td_act,td_pas);

    % trim to just movements
    td_bin = trimTD(td_bin,{'idx_movement_on',0},{'idx_movement_on',11});

    % remove low firing neurons
    td_bin = removeBadNeurons(td_bin,struct('min_fr',1));
    
    % add firing rates in addition to spike counts
    td_bin = addFiringRates(td_bin,struct('array',arrayname));

    % find average over the movement
    td_bin = binTD(td_bin,'average');

    % temporary hack to see what happens if only tuned neurons are used...
    [~,session_table] = getNTidx(pdTable,...
        'monkey','Chips',...
        'date','2017/9/13',...
        'act_movevecTuned',true,...
        'pas_movevecTuned',true);
    % which_units = find(ismember(td_bin(1).([arrayname '_unit_guide']),session_table.signalID,'rows'));
    which_units = 'all';
    
    %% find separabilities
    % suppress getTDfields warning...
    getTDfields(td_bin,'time');
    onetime_warn = warning('query','last'); 
    warning('off',onetime_warn.identifier)
    
    sepResults = actpasSep(td_bin,struct(...
        'num_repeats',num_repeats,...
        'num_folds',num_folds,...
        'model_aliases',{model_aliases},...
        'model_type',model_type,...
        'which_units',which_units,...
        'num_musc_pcs',num_musc_pcs,...
        'num_pcs',num_pcs));

    % turn warning back on
    warning('on',onetime_warn.identifier)

    %% save
    save(fullfile(savedir,strrep(filenames{filenum},'_TD.mat',savesuffix)),'sepResults')
end

%% bar plot for separabilities of models
%     % bar colors
%     bar_colors(1,:) = [0.5 0.5 0.5];
%     bar_colors(2,:) = [247, 148, 30]/255;
%     bar_colors(3,:) = [247, 192, 30]/255;
%     bar_colors(4,:) = [193, 25, 47]/255;
% 
%     % make plot
%     datadir = '~/Projects/limblab/data-td/ForceKin/Results/separability';
%     filenames = {'Han_20170203_actpasSeps_run20180818.mat','Chips_20170913_actpasSeps_run20180818.mat'};
%     num_monks = length(filenames);
%     monk_x = (2:3:((num_monks-1)*3+2))/10;
%     template_x = linspace(-0.5,0.5,4)/10;
%     model_spacing = mode(diff(template_x));
%     figure('defaultaxesfontsize',18)
%     for monkeynum = 1:num_monks
%         load(fullfile(datadir,filenames{monkeynum}))
% 
%         % calculate stats
%         mean_seps = mean(seps{:,:});
%         var_seps = var(seps{:,:});
%         correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
%         std_err_seps = sqrt(var_seps*correction);
% 
%         for modelnum = 1:length(mean_seps)
%             xval = monk_x(monkeynum) + template_x(modelnum);
%             bar(xval,mean_seps(modelnum),model_spacing,'facecolor',bar_colors(modelnum,:),'edgecolor','none')
%             hold on
%             plot([xval xval],[mean_seps(modelnum)-std_err_seps(modelnum) mean_seps(modelnum)+std_err_seps(modelnum)],'k','linewidth',3)
%         end
%     end
%     plot([0 monk_x(end)+0.2],[0.5 0.5],'--k','linewidth',2)
%     plot([0 monk_x(end)+0.2],[1 1],'--k','linewidth',2)
%     set(gca,'box','off','tickdir','out',...
%         'xtick',monk_x,'xticklabel',filenames,'xlim',[0 monk_x(end)+0.2],...
%         'ytick',[0 0.5 1],'yticklabel',{'','50%','100%'},...
%         'ticklabelinterpreter','none')
% 
