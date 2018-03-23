function plotMultiworkspaceTuning(results,neuron_idx)
%% Get tuning surfaces
    cm_viridis = viridis(200);
    model_names = {'glm_ext_model','glm_ego_model','glm_musc_model','S1_FR'};
    % model_names = {'glm_ext_model','glm_ego_model','glm_musc_model','S1_spikes'};
    model_alias = {'Extrinsic','Egocentric','Muscle','Real'};
    space_alias = {'PM','DL'};
    % for each model and space
    num_neurons = length(results.td_train(1).S1_unit_guide);
    fr_total = [];
    for modelnum = 1:4
        for spacenum = 1:2
            % position plots
            td = results.td_test{spacenum};
            poss = cat(1,td.pos);
            fr = get_vars(td,{model_names{modelnum},neuron_idx});

            ax_pos(modelnum) = subplot(4,3,(modelnum-1)*3+1);
            % tesselateFR(poss(:,1),poss(:,2),fr);
            plotSmoothFR(poss(:,1),poss(:,2),fr);
            hold on

            % get velocies and firing rates
            vels = cat(1,td.vel);
            fr = get_vars(td,{model_names{modelnum},neuron_idx});
            fr_total = [fr_total; fr];

            % Plot tuning surface
            ax_vel(modelnum,spacenum) = subplot(4,3,(modelnum-1)*3+spacenum+1);
            % tesselateFR(vels(:,1),vels(:,2),fr);
            plotSmoothFR(vels(:,1),vels(:,2),fr,[-40 40]);

            % Plot preferred direction in white
            % pd = results.pdTables{spacenum,modelnum}.velPD(neuron_idx);
            % pdCI = results.pdTables{spacenum,modelnum}.velPDCI(neuron_idx,:);
            % hold on
            % plot([0 20*cos(pd)],[0 20*sin(pd)],'-w','linewidth',2)
            % patch([0 20*cos(pdCI(1)) 20*cos(pd) 20*cos(pdCI(2))],[0 20*sin(pdCI(1)) 20*sin(pd) 20*sin(pdCI(2))],'w','edgecolor','none','facealpha',0.6)
        end
    end
    linkaxes(ax_pos,'xy')
    linkaxes(ax_vel,'xy')

    clims = prctile(fr_total,[5 95]);
    for modelnum = 1:4
        subplot(4,3,(modelnum-1)*3+1)
        caxis(clims)
        axis tight
        xlabel(sprintf('%s Pos tuning',model_alias{modelnum}))
        for spacenum = 1:2
            subplot(4,3,(modelnum-1)*3+spacenum+1)
            caxis(clims)
            xlim([-40 40])
            ylim([-40 40])
            xlabel(sprintf('%s Vel tuning, Workspace: %s',model_alias{modelnum},space_alias{spacenum}))

        end
    end

    % subplot(4,3,2)
    % label = ['Neuron ' num2str(signalIDs.signalID(neuron_idx))];
    % title(label)

    colormap(cm_viridis)
