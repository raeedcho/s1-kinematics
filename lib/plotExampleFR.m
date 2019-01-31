function plotExampleFR(trial_data,params)

    neuron_idx = 0;
    models = {'ext','ego','musc','handelbow'};
    assignParams(who,params)

    model_colors = getModelColors(models);

    %% Example predictions
    % trial_data = smoothSignals(trial_data,struct('signals','S1_FR','kernel_SD',0.1));
    trial_to_plot = randperm(length(trial_data),5);

    trial_endtime = 0;
    for trialnum = 1:length(trial_to_plot)
        temp_spikes = get_vars(trial_data(trial_to_plot(trialnum)),{'S1_FR',neuron_idx});
        bin_size = trial_data(trial_to_plot(trialnum)).bin_size;
        trial_starttime = trial_endtime;
        trial_endtime = trial_starttime + length(temp_spikes)*bin_size;
        plot(trial_starttime:bin_size:(trial_endtime-bin_size),...
            temp_spikes','k','linewidth',2)
        hold on

        for modelnum = 1:length(models)
            temp_pred = get_vars(trial_data(trial_to_plot(trialnum)),{sprintf('glm_%s_model',models{modelnum}),neuron_idx});
            plot(trial_starttime:bin_size:trial_endtime-bin_size,...
                temp_pred','color',model_colors(modelnum,:),'linewidth',2)

            % ax1 = subplot(2,1,1);
            % plot(temp_vel(:,1),'b','linewidth',2)
            % hold on
            % plot(temp_vel(:,2),'g','linewidth',2)
            % set(gca,'box','off','tickdir','out')
            % 
            % ax2 = subplot(2,1,2);
        end
        ylims = get(gca,'ylim');
        plot([trial_starttime trial_starttime],ylims,'--k','linewidth',2)
        plot([trial_endtime trial_endtime],ylims,'--k','linewidth',2)
    end
    set(gca,'box','off','tickdir','out')
end
