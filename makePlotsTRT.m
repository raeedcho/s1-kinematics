%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes plots from results given by analyzeTRT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
results = analyzeTRT(trial_data)

%%%%%%%%%%%%%%%%% Main line figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1 - Task and classic analysis methods
    % 1a - Monkey using manipulandum (in illustrator)

    % 1b - example movements in two workspaces (with neural firing dots?)
    % First, get example trials
    [~,td_pm_ex] = getTDidx(trial_data,'spaceNum',1,'result','R','rand',1);
    [~,td_dl_ex] = getTDidx(trial_data,'spaceNum',2,'result','R','rand',1);
    % trim to just go from target start to end
    td_ex = trimTD([td_pm_ex td_dl_ex],{'idx_ctHoldTime',0},{'idx_endTime',0});
    % plot the example trials
    figure
    plotTRTTrials(td_ex);
    % plot neural firing?
    unit_idx = 1
    plotSpikesOnHandle(td_ex,struct('unit_idx',unit_idx,'spikespec','b.','spikesize',10));
    % plot of same muscle movement given different Jacobians?

    % 1c - example directional rasters and tuning curves?

%% Figure 2 - Analysis block diagrams
    % 2a - Block diagram of three different models

    % 2b - Breaking up the data into training and testing sets

    % 2c - Example neural predictions for each model
    % neuron 1 has decent pseudo R2
    for neuron_idx = 1:length(td_eval{2,1})
        h = figure;
        temp_vel = cat(1,td_test{2}.vel);
        temp_spikes = cat(1,td_test{2}.(neural_signals));
        temp_pred_ext = cat(1,td_test{2}.(model_names{1}));
        temp_pred_musc = cat(1,td_test{2}.(model_names{3}));

        clf
        ax1 = subplot(2,1,1);
        plot(temp_vel(:,1),'b','linewidth',2)
        hold on
        plot(temp_vel(:,2),'g','linewidth',2)
        set(gca,'box','off','tickdir','out')

        ax2 = subplot(2,1,2);
        plot(temp_spikes(:,neuron_idx),'k','linewidth',2)
        hold on
        plot(temp_pred_ext(:,neuron_idx),'r','linewidth',2)
        plot(temp_pred_musc(:,neuron_idx),'c','linewidth',2)
        set(gca,'box','off','tickdir','out')

        linkaxes([ax1 ax2],'x')
        waitfor(h)
    end
    clearvars neuron_idx temp_* ax1 ax2

%% Figure 3 - Comparison of modeled tuning curves
    % 3a - a few example flat tuning curves (models on top of each other) DL under PM

    % 3b - DNA plot for each model/real tuning for each monkey

%% Figure 4 - Summaries
    % 3c - Summary over neurons: PD shift plots for each model and each monkey
    % 3d - Summary over monkeys/sessions: Bar chart of variance explained for each model by 1:1 line

%%%%%%%%%%%%%%%%% Supplementary Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary A - pR2s against each other for muscle vs endpoint

%% Supplementary B - Figure showing why one monkey didn't have as large a change in tuning

%% Plot handle positions
    if verbose
        figure
        pos_dl = cat(1,results.td_test{2}.pos);
        plot(pos_dl(:,1),pos_dl(:,2),'r')
        hold on
        pos_pm = cat(1,results.td_test{1}.pos);
        plot(pos_pm(:,1),pos_pm(:,2),'b')
        axis equal

        % clean up
        clearvars pos_*
    end

%% Plot example neuron model fit

%% Plot pR2s against each other
    % av_pR2_ext_dl = mean(td_ext_dl_eval,2);
    % av_pR2_musc_dl = mean(td_musc_dl_eval,2);
    % av_pR2_ext_pm = mean(td_ext_pm_eval,2);
    % av_pR2_musc_pm = mean(td_musc_pm_eval,2);
    av_pR2_ext = mean(td_eval{end,1},2);
    av_pR2_musc = mean(td_eval{end,3},2);
    
    % good_neurons = td_ext_dl_eval(:,1) > 0 & td_ext_pm_eval(:,1) > 0 & td_musc_pm_eval(:,1) > 0 & td_musc_dl_eval(:,1) > 0;
    % 
    % figure
    % plot(av_pR2_ext_dl(good_neurons),av_pR2_musc_dl(good_neurons),'ro','linewidth',2)
    % hold on
    % plot([-1 1],[-1 1],'k--','linewidth',2)
    % plot([0 0],[-1 1],'k-','linewidth',2)
    % plot([-1 1],[0 0],'k-','linewidth',2)
    % set(gca,'box','off','tickdir','out','xlim',[-0.1 0.5],'ylim',[-0.1 0.5])
    % 
    % figure
    % plot(av_pR2_ext_pm(good_neurons),av_pR2_musc_pm(good_neurons),'bo','linewidth',2)
    % hold on
    % plot([-1 1],[-1 1],'k--','linewidth',2)
    % plot([0 0],[-1 1],'k-','linewidth',2)
    % plot([-1 1],[0 0],'k-','linewidth',2)
    % set(gca,'box','off','tickdir','out','xlim',[-0.1 0.5],'ylim',[-0.1 0.5])

    % good_neurons = td_ext_eval(:,1) > 0 & td_musc_eval(:,1) > 0;
    good_neurons = isTuned{4};
    figure
    plot(repmat(av_pR2_ext(good_neurons)',2,1),td_eval{end,3}(good_neurons,:)','b-','linewidth',2)
    hold on
    plot(td_eval{end,1}(good_neurons,:)',repmat(av_pR2_musc(good_neurons)',2,1),'b-','linewidth',2)
    plot([-1 1],[-1 1],'k--','linewidth',2)
    plot([0 0],[-1 1],'k-','linewidth',2)
    plot([-1 1],[0 0],'k-','linewidth',2)
    set(gca,'box','off','tickdir','out','xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
    xlabel 'Hand-based pR2'
    ylabel 'Muscle-based pR2'

    % clean up
    clearvars av_pR2* good_neurons

%% Plot comparison of actual tuning curves with various modeled tuning curves
    % Get PDs and tuning curves for the modeled and actual neurons
    tuning_curves = cell(2,4); % PM is first row, DL is second. Column order is Ext, Ego, Musc, Real
    pdTables = cell(2,4);
    tuning_params = cell(1,4);
    
    % get PDs
    pdConvertedTable = getPDsFromWeights(results.tuningTable);

    num_bins = 8;
    % get PDs and tuning curves
    for modelnum = 1:4
        tuning_params{modelnum} = struct('num_bins',num_bins,'out_signals',{model_names(modelnum)},'out_signal_names',results.td_train(1).S1_unit_guide);
    
        for spacenum = 1:2
            % move converted PD table into separated cells with table selections
            [~,temp] = getNTidx(pdConvertedTable,'spaceNum',spacenum);
            % select only columns corresponding to the key or the model in question
            key_cols = ~contains(temp.Properties.VariableNames,'baseline') & ~contains(temp.Properties.VariableNames,'vel');
            model_cols = contains(temp.Properties.VariableNames,model_names{modelnum});
            temp = temp(:,key_cols | model_cols);
            % rename columns to get rid of model names
            temp.Properties.VariableNames = strrep(temp.Properties.VariableNames,[model_names{modelnum} '_'],'');
            pdTables{spacenum,modelnum} = temp;

            tuning_curves{spacenum,modelnum} = getTuningCurves(td_test{spacenum},tuning_params{modelnum});
        end
        % isTuned{modelnum} = checkIsTuned(pdTables{1,modelnum}) & checkIsTuned(pdTables{2,modelnum});
    end

    % first compare PM and DL tuning for each model
    for modelnum = 1:4
        figure;compareTuning(tuning_curves(:,modelnum),pdTables(:,modelnum))
        % figure;compareTuning(tuning_curves(:,modelnum),pdTables(:,modelnum),find(results.isTuned{4}))
    end

    % then compare PM and DL tuning for each model
    % reorder for color consistency..
    for spacenum = 1:2
        figure;compareTuning(tuning_curves(spacenum,[3,1,2,4]),pdTables(spacenum,[3,1,2,4]),find(results.isTuned{4}))
    end

%% Plot DL vs PM just for neurons actually tuned to velocity
    figure
    comparePDs(pm_real_pdTable(isTuned_real,:),dl_real_pdTable(isTuned_real,:),struct('move_corr','vel'),'ko','linewidth',2)
    hold on
    comparePDs(pm_ext_pdTable(isTuned_real,:),dl_ext_pdTable(isTuned_real,:),struct('move_corr','vel'),'ro','linewidth',2)
    comparePDs(pm_ego_pdTable(isTuned_real,:),dl_ego_pdTable(isTuned_real,:),struct('move_corr','vel'),'go','linewidth',2)
    comparePDs(pm_musc_pdTable(isTuned_real,:),dl_musc_pdTable(isTuned_real,:),struct('move_corr','vel'),'bo','linewidth',2)
    xlabel 'PM preferred direction'
    ylabel 'DL preferred direction'
    
%% Plot PD shift clouds for each neuron individually
    colors = {'r','g','b'};
    titles = {'Hand-based model PD shift vs Actual PD shift','Egocentric model PD shift vs Actual PD shift','Muscle-based model PD shift vs Actual PD shift'};
    for modelnum = 1:3
        comparePDClouds(results.shift_tables{4},results.shift_tables{modelnum},struct('filter_tuning',1),colors{modelnum},'facealpha',0.5)
        xlabel 'Actual PD Shift'
        ylabel 'Modeled PD Shift'
        title(titles{modelnum})
    end

    % clean up
    clearvars colors titles

%% Plot tuning weight clouds
    tuningHull = getTuningHull(results.tuningTable);
    % loop over each unit in one workspace
    % n_rows = ceil(sqrt(height(signalIDs)+1));
    cloud_fig = figure;
    surf_fig = figure;
    for neuron_idx = 10%1:height(signalIDs)
        close_fig = figure;

        figure(cloud_fig)
        clf
        plotMWTuningCloud(tuningHull,neuron_idx)

        figure(surf_fig)
        clf
        plotMultiworkspaceTuning(results,neuron_idx)

        waitfor(close_fig)
    end
