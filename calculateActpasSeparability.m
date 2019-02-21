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
savesuffix = '_separationResults_allModels_run20190127.mat';

model_aliases = {'ext','extforce','handle_ext','joint','musc','handelbow'};
arrayname = 'S1';
num_musc_pcs = 5;

%% Loop through files
for filenum = 1:length(filenames)
    clear sepResults

    %% load data
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

    % remove low firing neurons
    td = removeBadNeurons(td,struct('min_fr',0.1));
    
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
                'bumpDir',...
                'tgtDir',...
                'ctrHoldBump',...
                'ctrHold',...
                }},...
            'start_name','idx_startTime',...
            'end_name','idx_endTime'));
    [~,td] = getTDidx(td,'result','R');
    td = reorderTDfields(td);

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
    
    % suppress getTDfields warning...
    getTDfields(td,'time');
    onetime_warn = warning('query','last'); 
    warning('off',onetime_warn.identifier)
    
    sepResults = actpasSep(td,struct(...
        'num_repeats',20,...
        'num_folds',5,...
        'model_aliases',{model_aliases},...
        'num_musc_pcs',num_musc_pcs));

    % turn warning back on
    warning('on',onetime_warn.identifier)

    %% save
    save(fullfile(savedir,strrep(filenames{filenum},'_TD.mat',savesuffix)),'sepResults')
end

%% bar plot for separabilities of models
    % bar colors
    bar_colors(1,:) = [0.5 0.5 0.5];
    bar_colors(2,:) = [247, 148, 30]/255;
    bar_colors(3,:) = [247, 192, 30]/255;
    bar_colors(4,:) = [193, 25, 47]/255;

    % make plot
    datadir = '~/Projects/limblab/data-td/ForceKin/Results/separability';
    filenames = {'Han_20170203_actpasSeps_run20180818.mat','Chips_20170913_actpasSeps_run20180818.mat'};
    num_monks = length(filenames);
    monk_x = (2:3:((num_monks-1)*3+2))/10;
    template_x = linspace(-0.5,0.5,4)/10;
    model_spacing = mode(diff(template_x));
    figure('defaultaxesfontsize',18)
    for monkeynum = 1:num_monks
        load(fullfile(datadir,filenames{monkeynum}))

        % calculate stats
        mean_seps = mean(seps{:,:});
        var_seps = var(seps{:,:});
        correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
        std_err_seps = sqrt(var_seps*correction);

        for modelnum = 1:length(mean_seps)
            xval = monk_x(monkeynum) + template_x(modelnum);
            bar(xval,mean_seps(modelnum),model_spacing,'facecolor',bar_colors(modelnum,:),'edgecolor','none')
            hold on
            plot([xval xval],[mean_seps(modelnum)-std_err_seps(modelnum) mean_seps(modelnum)+std_err_seps(modelnum)],'k','linewidth',3)
        end
    end
    plot([0 monk_x(end)+0.2],[0.5 0.5],'--k','linewidth',2)
    plot([0 monk_x(end)+0.2],[1 1],'--k','linewidth',2)
    set(gca,'box','off','tickdir','out',...
        'xtick',monk_x,'xticklabel',filenames,'xlim',[0 monk_x(end)+0.2],...
        'ytick',[0 0.5 1],'yticklabel',{'','50%','100%'},...
        'ticklabelinterpreter','none')

