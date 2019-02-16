%% Set up meta info
if ispc
    homefolder = 'C:\Users\rhc307';
else
    homefolder = '/home/raeed';
end

datadir = fullfile(homefolder,'data','project-data','limblab','s1-kinematics','td-library');
file_info = dir(fullfile(datadir,'*TRT*'));
filenames = horzcat({file_info.name})';

arrayname = 'S1';
num_trials = 2;
% colors for pm, dl conditions
cond_colors = [...
    231,138,195;...
    166,216,84]/255;

%% Loop through files
for filenum = 4%1:length(filenames)
    %% Load data
    td = load(fullfile(datadir,[filenames{filenum}]));

    % rename trial_data for ease
    td = td.trial_data;

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

    %% choose a random few trials and plot
    figure('defaultaxesfontsize',18)
    max_x = 0;
    for spacenum = 1:2
        [~,td_temp] = getTDidx(td,'spaceNum',spacenum,'rand',num_trials);
        spikes = getSig(td_temp,'S1_spikes')';
        timevec = (1:size(spikes,2))*td_temp(1).bin_size;
        subplot(1,2,spacenum)
        for neuronnum = 1:size(spikes,1)
            spike_times = timevec(spikes(neuronnum,:)>0);
            scatter(spike_times,repmat(neuronnum,size(spike_times)),5,'k','filled')
            hold on
        end
        trial_end = 0;
        for trialnum = 1:num_trials
            plot(td_temp(trialnum).bin_size*repmat(td_temp(trialnum).idx_otHoldTime,2,1)+trial_end,...
                repmat([0;size(spikes,1)],1,length(td_temp(trialnum).idx_otHoldTime)),...
                '--','color',cond_colors(spacenum,:))
            plot(td_temp(trialnum).bin_size*repmat(td_temp(trialnum).idx_targetStartTime,2,1)+trial_end,...
                [0;size(spikes,1)],...
                '--k')
            plot(td_temp(trialnum).bin_size*repmat(td_temp(trialnum).idx_endTime,2,1)+trial_end,...
                [0;size(spikes,1)],...
                '--k')
            trial_end = trial_end + (td_temp(trialnum).idx_endTime-1)*td_temp(trialnum).bin_size;
        end
        max_x = max(max_x,trial_end);
        xlabel 'Time (s)'
        set(gca,'box','off','tickdir','out')
    end
    for spacenum = 1:2
        subplot(1,2,spacenum)
        xlim([0 max_x+0.5]);
    end
end

