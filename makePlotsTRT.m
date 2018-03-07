%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes plots from results given by analyzeTRT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot handle positions
    if verbose
        figure
        pos_dl = cat(1,td_test{2}.pos);
        plot(pos_dl(:,1),pos_dl(:,2),'r')
        hold on
        pos_pm = cat(1,td_test{1}.pos);
        plot(pos_pm(:,1),pos_pm(:,2),'b')
        axis equal

        % clean up
        clearvars pos_*
    end

%% Plot example neuron model fit
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

    % clean up
    clearvars neuron_idx temp_* ax1 ax2

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
    % first compare PM and DL tuning for each model
    for modelnum = 1:4
        figure;compareTuning(tuning_curves(:,modelnum),pdTables(:,modelnum),find(isTuned{4}))
    end

    % then compare PM and DL tuning for each model
    % reorder for color consistency..
    for spacenum = 1:2
        figure;compareTuning(tuning_curves(spacenum,[3,1,2,4]),pdTables(spacenum,[3,1,2,4]),find(isTuned{4}))
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
        comparePDClouds(shift_tables{4},shift_tables{modelnum},struct('filter_tuning',1),colors{modelnum},'facealpha',0.5)
        xlabel 'Actual PD Shift'
        ylabel 'Modeled PD Shift'
        title(titles{modelnum})
    end

    % clean up
    clearvars colors titles
