%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% actpas_analysis - 
%%      script to run active-passive task analyses in manuscript
%%      and plot out relevant details
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up meta info
    if ispc
        dataroot = 'G:\raeed\project-data\limblab\s1-kinematics';
    else
        dataroot = '/data/raeed/project-data/limblab/s1-kinematics';
    end
    
    file_info = dir(fullfile(dataroot,'reaching_experiments','*COactpas*.mat'));
    filenames = horzcat({file_info.name})';

    % plotting variables
    monkey_names = {'C','H'};
    arrayname = 'S1';
    neural_signals = [arrayname '_FR'];
    model_aliases = {'ext','handelbow'};
    models_to_plot = {neural_signals,'ext','handelbow'};
    model_titles = {'Actual Firing','Hand-only','Whole-arm'};
    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;

%% Loop through trial data files to clean them up
    trial_data_cell = cell(1,length(filenames));
    for filenum = 1:length(filenames)
        %% load and preprocess data
        td = load(fullfile(dataroot,'reaching_experiments',[filenames{filenum}]));
    
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
        
        % for C_20170912, trial structure is such that active and passive are part of the same trial--split it up
        if strcmpi(td(1).monkey,'C') && contains(td(1).date_time,'2017/9/12')
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
        td_trim = cat(2,td_act,td_pas);
    
        % remove low firing neurons
        td_trim = removeBadNeurons(td_trim,struct(...
            'min_fr',1,...
            'fr_window',{{'idx_movement_on',0;'idx_movement_on',11}},...
            'calc_fr',true));
    
        trial_data_cell{filenum} = td_trim;
    end

%% Plot trial info (hand speed and example rasters)
    for filenum = 1:length(trial_data_cell)
        %% load and preprocess data
        td = trial_data_cell{filenum};
    
        % trim to just movements
        td = trimTD(td,{'idx_movement_on',-50},{'idx_movement_on',60});
    
        %% Plot out hand speed
        figure('defaultaxesfontsize',18)
        for trial = 1:length(td)
            timevec = ((1:length(td(trial).vel_norm))-td(trial).idx_movement_on)*td(trial).bin_size;
            if td(trial).ctrHoldBump
                plot(timevec,td(trial).vel_norm,'r')
            else
                plot(timevec,td(trial).vel_norm,'k')
            end
            hold on
        end
        plot(zeros(2,1),ylim,'--k','linewidth',2)
        hold on
        plot(repmat(0.12,2,1),ylim,'--k','linewidth',2)
        xlabel('Time from movement onset (s)')
        ylabel('Hand speed (cm/s)')
        set(gca,'box','off','tickdir','out','xtick',[-0.5 0 0.12 0.5])
        set(gcf,'renderer','Painters')
        suptitle(sprintf('Monkey %s %s',td(1).monkey, td(1).date_time))
    
        %% Plot out example rasters for each direction
        dirs = unique(cat(1,td.tgtDir));
        figure('defaultaxesfontsize',18)
        for dirnum = 1:length(dirs)
            % pick a random active and random passive trial with this direction
            act_idx = getTDidx(td,'tgtDir',dirs(dirnum),'ctrHoldBump',false,'rand',1);
            pas_idx = getTDidx(td,'bumpDir',dirs(dirnum),'ctrHoldBump',true,'rand',1);
            td_temp = td([act_idx pas_idx]);
    
            for trialnum = 1:length(td_temp)
                spikes = getSig(td_temp(trialnum),'S1_spikes')';
                timevec = ((1:size(spikes,2))-td_temp(trialnum).idx_movement_on)*td_temp(trialnum).bin_size;
                % active on left, passive on right
                subplot(length(dirs),length(td_temp),(dirnum-1)*length(td_temp)+trialnum)
                % neurons
                for neuronnum = 1:size(spikes,1)
                    spike_times = timevec(spikes(neuronnum,:)>0);
                    scatter(spike_times,repmat(neuronnum,size(spike_times)),5,'k','filled')
                    hold on
                end
                plot(zeros(1,2),[0 size(spikes,1)+1],'--k')
                plot(ones(1,2)*0.12,[0 size(spikes,1)+1],'--k')
                xlabel('Time from movement onset (s)')
                set(gca,'box','off','tickdir','out','xtick',[-0.5 0 0.12 0.5],'ytick',[])
            end
            subplot(length(dirs),length(td_temp),(dirnum-1)*length(td_temp)+1)
            ylabel(sprintf('Direction %f',dirs(dirnum)))
        end
        suptitle(sprintf('Monkey %s %s',td(1).monkey, td(1).date_time))
    end

%% Loop through results to pull out relevant info
    neuron_eval_cell =cell(length(trial_data_cell),1);
    fileclock = tic;
    fprintf('Starting analysis for %d files. This will take a while...',length(trial_data_cell))
    for filenum = 1:length(trial_data_cell)
        td = trial_data_cell{filenum};
    
        % trim to just movements
        td = trimTD(td,{'idx_movement_on',0},{'idx_movement_on',11});
    
        % check to make sure all neurons fire at least once in each condition (pretty rare that one doesn't)
        [~,td_act] = getTDidx(td,'ctrHoldBump',false);
        [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
        firing_units = mean(getSig(td_act,'S1_spikes'))~=0 & mean(getSig(td_pas,'S1_spikes'))~=0;
        if any(~firing_units)
            unit_ids = td(1).([arrayname '_unit_guide']);
            new_unit_guide = unit_ids(firing_units,:);
            
            for trialnum = 1:length(td)
                td(trialnum).(sprintf('%s_unit_guide',arrayname)) = new_unit_guide;
                
                spikes = td(trialnum).(sprintf('%s_spikes',arrayname));
                spikes(:,~firing_units) = [];
                td(trialnum).(sprintf('%s_spikes',arrayname)) = spikes;
            end
            fprintf('Removed %d neurons for not firing in one condition\n',sum(~firing_units))
        end
        
        % add firing rates in addition to spike counts
        td = addFiringRates(td,struct('array',arrayname));
    
        % find average over the movement
        td = binTD(td,'average');
    
        %% find separabilities
        % suppress getTDfields warning...
        getTDfields(td,'time');
        onetime_warn = warning('query','last'); 
        warning('off',onetime_warn.identifier)
        
        sepResults = actpasSep(td,struct(...
            'neural_signals',[arrayname '_FR'],...
            'model_aliases',{model_aliases}));
    
        % turn warning back on
        warning('on',onetime_warn.identifier)

        % extract neuron_eval_table and trial_table
        % replace infs with nans
        numeric_cols = strcmpi(sepResults.neuron_eval_table.Properties.VariableDescriptions,'linear');
        numeric_vals = sepResults.neuron_eval_table(:,numeric_cols).Variables;
        infidx = isinf(numeric_vals);
        numeric_vals(infidx) = NaN;
        sepResults.neuron_eval_table(:,numeric_cols).Variables = numeric_vals;

        % compile neuron eval table together
        neuron_eval_cell{filenum} = sepResults.neuron_eval_table;

        % extract only the columns we want to keep
        neuron_eval_cell{filenum}.Properties.VariableNames = strrep(neuron_eval_cell{filenum}.Properties.VariableNames,'glm_','');
        neuron_eval_cell{filenum}.Properties.VariableNames = strrep(neuron_eval_cell{filenum}.Properties.VariableNames,'model_','');
        cols_to_keep = [...
            {'monkey','date','task','signalID','crossvalID'},...
            strcat(models_to_plot(2:end),'_eval'),...
            strcat(models_to_plot(2:end),'_act_eval'),...
            strcat(models_to_plot(2:end),'_pas_eval'),...
            strcat(models_to_plot(2:end),'_train_act_eval'),...
            strcat(models_to_plot(2:end),'_train_pas_eval'),...
            strcat(models_to_plot(2:end),'_half_full_train_act_eval'),...
            strcat(models_to_plot(2:end),'_half_full_train_pas_eval'),...
            strcat(models_to_plot(1),'_indiv_sep')];
            
        neuron_eval_cell{filenum} = neuron_eval_cell{filenum}(:,cols_to_keep);

        % make a histogram plot of neural firing rates for active and passive trials
        session_trials = neuronAverage(sepResults.trial_table,struct(...
            'keycols',{{'monkey','date_time','task','trialID','isPassive'}},...
            'do_ci',false));
        figure('defaultaxesfontsize',18)
        ax = zeros(size(session_trials.S1_FR,2),1);
        % here we know that there are only a finite number of possibilities
        possible_FR = unique(session_trials.S1_FR);
        % split into active and passive
        [~,act_trials] = getNTidx(session_trials,'isPassive',false);
        [~,pas_trials] = getNTidx(session_trials,'isPassive',true);
        for neuronnum = 1:size(session_trials.S1_FR,2)
            ax(neuronnum) = subplot(1,size(session_trials.S1_FR,2),neuronnum);
  
            % get counts of fr in the unique bins
            act_counts = histcounts(act_trials.S1_FR(:,neuronnum),[possible_FR;Inf]);
            pas_counts = histcounts(pas_trials.S1_FR(:,neuronnum),[possible_FR;Inf]);
  
            % plot bars for each
            % plot([0 0],[possible_FR(1) possible_FR(end)])
            barh(possible_FR',act_counts,1,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5)
            hold on
            barh(possible_FR',pas_counts,1,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5)
  
            % set(gca,'box','off','tickdir','out')
            axis off
        end
        subplot(1,size(session_trials.S1_FR,2),1)
        axis on
        set(gca,'box','off','tickdir','out','xtick',[])
        ylabel('Firing rate (Hz)')
        suptitle(sprintf('%s %s',session_trials.monkey{1},session_trials.date_time{1}))
        linkaxes(ax,'y')

        % output a counter
        fprintf('Processed file %d of %d at time %f\n',filenum,length(trial_data_cell),toc(fileclock))
    end

    % compile and average
    neuron_eval = vertcat(neuron_eval_cell{:});
    avg_neuron_eval = neuronAverage(neuron_eval,struct(...
        'keycols',{{'monkey','date','task','signalID'}},...
        'do_ci',false,...
        'do_nanmean',true));

%% make plots
    % plot separability of each neuron and save CIs into avg_neuron_eval
    for monkeynum = 1:length(monkey_names)
        [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
        session_dates = unique(monkey_evals.date);
        for sessionnum = 1:length(session_dates)
            [~,session_evals] = getNTidx(monkey_evals,'date',session_dates{sessionnum});
            [~,session_avg] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});

            % create place to save CIs
            signalID = session_avg.signalID;
            S1_FR_indiv_sep_CI_lo = zeros(size(signalID,1),1);
            S1_FR_indiv_sep_CI_hi = zeros(size(signalID,1),1);

            figure('defaultaxesfontsize',18)
            plot([0 size(signalID,1)+1],[0.5 0.5],'--k','linewidth',2)
            hold on
            for neuronnum = 1:size(signalID,1)
                [~,single_neuron_eval] = getNTidx(session_evals,'signalID',signalID(neuronnum,:));

                % figure out error bars and save
                [CI_lo,CI_hi] = crossval_errorbars(single_neuron_eval.S1_FR_indiv_sep,struct(...
                    'num_repeats',double(max(single_neuron_eval.crossvalID(:,1))),...
                    'num_folds',double(max(single_neuron_eval.crossvalID(:,2)))));
                S1_FR_indiv_sep_CI_lo(neuronnum,:) = CI_lo;
                S1_FR_indiv_sep_CI_hi(neuronnum,:) = CI_hi;
                
                % scatter(repmat(neuronnum,1,height(single_neuron_eval)),single_neuron_eval.S1_FR_indiv_sep,25,'k','filled','markerfacealpha',0.2)
                plot(repmat(neuronnum,1,2),[CI_lo CI_hi],'-k','linewidth',2)
                scatter(neuronnum,session_avg.S1_FR_indiv_sep(neuronnum,:),100,'k','filled')
            end
            set(gca,'box','off','tickdir','out','ylim',[0 1],'xlim',[0 size(signalID,1)+1])

            % save error bars into avg_neuron_eval table
            CI_cell{monkeynum,sessionnum} = table(...
                session_avg.monkey,...
                session_avg.date,...
                signalID,...
                S1_FR_indiv_sep_CI_lo,...
                S1_FR_indiv_sep_CI_hi,...
                'VariableNames',{'monkey','date','signalID','S1_FR_indiv_sep_CI_lo','S1_FR_indiv_sep_CI_hi'});
        end
    end
    CI_table = vertcat(CI_cell{:});
    avg_neuron_eval = join(avg_neuron_eval,CI_table);

    % compare pR2 of handelbow vs ext
    figure('defaultaxesfontsize',18)
    model_pairs = {'ext','handelbow'};
    for pairnum = 1:size(model_pairs,1)
        for monkeynum = 1:length(monkey_names)
            % set subplot...
            subplot(size(model_pairs,1),length(monkey_names),...
                (pairnum-1)*length(monkey_names)+monkeynum)
            plot([-0.4 0.6],[-0.4 0.6],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-0.4 0.6],'k-','linewidth',0.5)
            plot([-0.4 0.6],[0 0],'k-','linewidth',0.5)

            % get sessions
            [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
            session_dates = unique(monkey_evals.date);

            for sessionnum = 1:length(session_dates)
                [~,session_evals] = getNTidx(monkey_evals,'date',session_dates{sessionnum});
                pr2_winners = compareEncoderMetrics(session_evals,struct(...
                    'bonferroni_correction',6,...
                    'models',{models_to_plot},...
                    'model_pairs',{model_pairs},...
                    'postfix','_eval'));

                [~,avg_pR2] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});
                % scatter filled circles if there's a winner, empty circles if not
                no_winner =  cellfun(@isempty,pr2_winners(pairnum,:));
                scatterlims(...
                    [-0.4 0.6],...
                    [-0.4 0.6],...
                    avg_pR2.(strcat(model_pairs{pairnum,1},'_eval'))(no_winner),...
                    avg_pR2.(strcat(model_pairs{pairnum,2},'_eval'))(no_winner),...
                    [],session_colors(sessionnum,:))
                scatterlims(...
                    [-0.4 0.6],...
                    [-0.4 0.6],...
                    avg_pR2.(strcat(model_pairs{pairnum,1},'_eval'))(~no_winner),...
                    avg_pR2.(strcat(model_pairs{pairnum,2},'_eval'))(~no_winner),...
                    [],session_colors(sessionnum,:),'filled')
            end

            % make axes pretty
            set(gca,'box','off','tickdir','out')
            axis image
            if monkeynum ~= 1 || pairnum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,1})))
            ylabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,2})))
        end
    end

    % Plot within condition vs across condition pR2 for each neuron in all sessions
    conds = {'act','pas'};
    model_pairs = {'ext','handelbow'};
    for modelnum = 2:length(models_to_plot)
        figure('defaultaxesfontsize',18)
        for monkeynum = 1:length(monkey_names)
            for condnum = 1:2
                % set subplot
                subplot(2,length(monkey_names),(condnum-1)*length(monkey_names)+monkeynum)
                plot([-0.7 0.7],[-0.7 0.7],'k--','linewidth',0.5)
                hold on
                plot([0 0],[-0.7 0.7],'k-','linewidth',0.5)
                plot([-0.7 0.7],[0 0],'k-','linewidth',0.5)

                % get sessions
                [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
                session_dates = unique(monkey_evals.date);

                % plot out each session
                for sessionnum = 1:length(session_dates)
                    [~,avg_pR2] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});
                    
                    [~,session_evals] = getNTidx(monkey_evals,'date',session_dates{sessionnum});
                    pr2_winners = compareEncoderMetrics(session_evals,struct(...
                        'bonferroni_correction',6,...
                        'models',{models_to_plot},...
                        'model_pairs',{model_pairs},...
                        'postfix','_eval'));

                    % fill by whether separable or not?
                    sig_seps = avg_pR2.S1_FR_indiv_sep_CI_lo > 0.5;
                    
                    % fill by whether winner of comparison with ext or not?
                    no_winner =  cellfun(@isempty,pr2_winners(1,:));

                    scatterlims(...
                        [-0.7 0.7],...
                        [-0.7 0.7],...
                        avg_pR2.(sprintf('%s_%s_eval',models_to_plot{modelnum},conds{condnum}))(no_winner),...
                        avg_pR2.(sprintf('%s_train_%s_eval',models_to_plot{modelnum},conds{condnum}))(no_winner),...
                        [],session_colors(sessionnum,:),'filled')
                    scatterlims(...
                        [-0.7 0.7],...
                        [-0.7 0.7],...
                        avg_pR2.(sprintf('%s_%s_eval',models_to_plot{modelnum},conds{condnum}))(~no_winner),...
                        avg_pR2.(sprintf('%s_train_%s_eval',models_to_plot{modelnum},conds{condnum}))(~no_winner),...
                        [],session_colors(sessionnum,:),'filled')
                end
                % make axes pretty
                set(gca,'box','off','tickdir','out',...
                    'xlim',[-0.7 0.7],'ylim',[-0.7 0.7])
                axis equal
                % if monkeynum ~= 1 || condnum ~= 1
                %     set(gca,'box','off','tickdir','out',...
                %         'xtick',[],'ytick',[])
                % end
                xlabel(sprintf('%s pR2, trained full, tested %s',getModelTitles(models_to_plot{modelnum}),conds{condnum}))
                ylabel(sprintf('%s pR2, trained %s, tested %s',getModelTitles(models_to_plot{modelnum}),conds{condnum},conds{condnum}))
            end
        end
        % suptitle('Full pR^2 vs within condition pR^2')
    end

    % plot separability against full and within condition pR2
    conds = {'','act_','pas_'};
    for modelnum = 2:length(models_to_plot)
        figure('defaultaxesfontsize',18)
        for monkeynum = 1:length(monkey_names)
            for condnum = 1:length(conds)
                % set subplot
                subplot(length(monkey_names),length(conds),(monkeynum-1)*length(conds)+condnum)
                plot([0 0],[0 1],'k-','linewidth',0.5)
                hold on
                plot([-0.7 0.7],[0 0],'k-','linewidth',0.5)
                plot([-0.7 0.7],[0.5 0.5],'k--','linewidth',0.5)

                % get sessions
                [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
                session_dates = unique(monkey_evals.date);

                % plot out each session
                for sessionnum = 1:length(session_dates)
                    [~,session_evals] = getNTidx(monkey_evals,'date',session_dates{sessionnum});
                    [~,avg_pR2] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});

                    sig_seps = avg_pR2.S1_FR_indiv_sep_CI_lo > 0.5;

                    scatterlims(...
                        [-0.7 0.7],...
                        [0.4 1],...
                        avg_pR2.(sprintf('%s_%seval',models_to_plot{modelnum},conds{condnum}))(~sig_seps),...
                        avg_pR2.S1_FR_indiv_sep(~sig_seps),...
                        [],session_colors(sessionnum,:),'filled')
                    scatterlims(...
                        [-0.7 0.7],...
                        [0.4 1],...
                        avg_pR2.(sprintf('%s_%seval',models_to_plot{modelnum},conds{condnum}))(sig_seps),...
                        avg_pR2.S1_FR_indiv_sep(sig_seps),...
                        [],session_colors(sessionnum,:),'filled')

                    % fit quick linear model to plot fit line
                    % lm = fitlm(...
                    %     avg_pR2.(sprintf('%s_%seval',models_to_plot{modelnum},conds{condnum})),...
                    %     avg_pR2.S1_FR_indiv_sep);
                    % plot([-1;1],lm.predict([-1;1]),...
                    %     '--','color',session_colors(sessionnum,:),'linewidth',1)
                end
                % make axes pretty
                set(gca,'box','off','tickdir','out','xtick',[-0.7 0.7])
                if condnum ~= 1
                    set(gca,'ytick',[])
                else
                    ylabel('Neural Separability')
                end
                if monkeynum == length(monkey_names)
                    xlabel(sprintf('%s %s pR2',getModelTitles(models_to_plot{modelnum}),conds{condnum}))
                end
                axis image
                set(gca,'ylim',[0.4 1])
            end
        end
        suptitle('Neural separability vs pR^2')
    end

    % get correlation values for each crossval
    keycols = {'monkey','date','task','crossvalID'};
    keyTable = unique(neuron_eval(:,keycols));
    corr_cell = cell(height(keyTable),1);
    for key_idx = 1:height(keyTable)
        key = keyTable(key_idx,:);
        cond_idx = ismember(neuron_eval(:,keycols),key);
        neuron_eval_select = neuron_eval(cond_idx,:);

        % get correlations
        model_corr = cell(1,length(models_to_plot)-1);
        for modelnum = 2:length(models_to_plot)
            corr_pr2_neuronsep = corr(neuron_eval_select.S1_FR_indiv_sep,neuron_eval_select.(sprintf('%s_eval',models_to_plot{modelnum})),'rows','complete');
            corr_actpr2_neuronsep = corr(neuron_eval_select.S1_FR_indiv_sep,neuron_eval_select.(sprintf('%s_act_eval',models_to_plot{modelnum})),'rows','complete');
            corr_paspr2_neuronsep = corr(neuron_eval_select.S1_FR_indiv_sep,neuron_eval_select.(sprintf('%s_pas_eval',models_to_plot{modelnum})),'rows','complete');

            model_corr{modelnum-1} = table(...
                corr_pr2_neuronsep,...
                corr_actpr2_neuronsep,...
                corr_paspr2_neuronsep,...
                'VariableNames',strcat(models_to_plot{modelnum},{...
                    '_corr_pr2_neuronsep',...
                    '_corr_actpr2_neuronsep',...
                    '_corr_paspr2_neuronsep'}));
        end

        % put together in table
        corr_cell{key_idx} = horzcat(model_corr{:});
    end
    neuron_corr_table = horzcat(keyTable,vertcat(corr_cell{:}));

    % make figure for correlations of pR2 handelbow model with separability
    alpha = 0.05;
    xvals = [2 5 8]/10;
    for modelnum = 2:length(models_to_plot)
        figure('defaultaxesfontsize',18)
        for monkeynum = 1:length(monkey_names)
            subplot(1,length(monkey_names),monkeynum)
            plot([min(xvals)-0.2 max(xvals)+0.2],[0 0],'-k','linewidth',2)
            hold on
            
            % figure out what sessions we have for this monkey
            [~,monkey_corrs] = getNTidx(neuron_corr_table,'monkey',monkey_names{monkeynum});
            session_dates = unique(monkey_corrs.date);

            for sessionnum = 1:length(session_dates)
                [~,session_corrs] = getNTidx(monkey_corrs,'date',session_dates{sessionnum});

                % estimate error bars
                [~,cols] = ismember(...
                    strcat(models_to_plot{modelnum},{'_corr_pr2_neuronsep','_corr_actpr2_neuronsep','_corr_paspr2_neuronsep'}),...
                    session_corrs.Properties.VariableNames);
                num_repeats = double(max(session_corrs.crossvalID(:,1)));
                num_folds = double(max(session_corrs.crossvalID(:,2)));
                crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
                yvals = mean(session_corrs{:,cols});
                var_corrs = var(session_corrs{:,cols});
                upp = tinv(1-alpha/2,num_folds*num_repeats-1);
                low = tinv(alpha/2,num_folds*num_repeats-1);
                CI_lo = yvals + low * sqrt(crossval_correction*var_corrs);
                CI_hi = yvals + upp * sqrt(crossval_correction*var_corrs);
                
                % plot dots and lines
                plot(repmat(xvals,2,1),[CI_lo;CI_hi],'-','linewidth',2,'color',session_colors(sessionnum,:))
                scatter(xvals(:),yvals(:),50,session_colors(sessionnum,:),'filled')
            end
            ylabel('Correlation with active/passive separability')
            title(sprintf('Monkey %s',monkey_names{monkeynum}))
            set(gca,'box','off','tickdir','out',...
                'xlim',[min(xvals)-0.2 max(xvals)+0.2],...
                'xtick',xvals,'xticklabel',{'Full pR^2','Active pR^2','Passive pR^2'},...
                'ylim',[-1 1],'ytick',[-1 -0.5 0 0.5 1])
        end
        suptitle(models_to_plot{modelnum})
    end

