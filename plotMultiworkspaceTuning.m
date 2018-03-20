function plotMultiworkspaceTuning(results,neuron_idx)
%% Get tuning surfaces
    cm_viridis = viridis(200);
    % model_names = {'glm_ext_model','glm_ego_model','glm_musc_model','S1_FR'};
    model_names = {'glm_ext_model','glm_ego_model','glm_musc_model','S1_spikes'};
    model_alias = {'Extrinsic','Egocentric','Muscle','Real'};
    space_alias = {'PM','DL'};
    % for each model and space
    num_neurons = length(results.td_train(1).S1_unit_guide);
    fr_total = [];
    for modelnum = 1:4
        for spacenum = 1:2
            % get velocies and firing rates
            ax(modelnum,spacenum) = subplot(4,2,(modelnum-1)*2+spacenum);
            td = results.td_test{spacenum};
            vels = cat(1,td.vel);
            fr = get_vars(td,{model_names{modelnum},neuron_idx});
            fr_total = [fr_total; fr];

            % Plot tuning surface
            % tesselateFR(vels(:,1),vels(:,2),fr);
            plotSmoothFR(vels(:,1),vels(:,2),fr,[-40 40]);

            % Plot preferred direction in white
            pd = results.pdTables{spacenum,modelnum}.velPD(neuron_idx);
            pdCI = results.pdTables{spacenum,modelnum}.velPDCI(neuron_idx,:);
            hold on
            plot([0 20*cos(pd)],[0 20*sin(pd)],'-w','linewidth',2)
            patch([0 20*cos(pdCI(1)) 20*cos(pd) 20*cos(pdCI(2))],[0 20*sin(pdCI(1)) 20*sin(pd) 20*sin(pdCI(2))],'w','edgecolor','none','facealpha',0.6)
        end
    end
    linkaxes(ax,'xy')

    clims = prctile(fr_total,[5 95]);
    for modelnum = 1:4
        for spacenum = 1:2
            subplot(4,2,(modelnum-1)*2+spacenum)
            caxis(clims)
            xlim([-40 40])
            ylim([-40 40])
            xlabel(sprintf('%s, Workspace: %s',model_alias{modelnum},space_alias{spacenum}))
        end
    end

    colormap(cm_viridis)
