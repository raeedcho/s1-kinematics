function plotMWTuningSurfaces(td_test,pdTables,neuron_idx,model_aliases)
%% Get tuning surfaces
    cm_viridis = viridis(200);
    model_type = 'glm';
    model_names = [strcat(model_type,'_',model_aliases,'_model') {'S1_FR'}];
    num_models = numel(model_names);
    model_titles = cell(1,num_models);
    for modelnum = 1:num_models-1
        switch model_aliases{modelnum}
        case 'musc'
            model_titles{modelnum} = 'Muscle-based';
        case 'ext'
            model_titles{modelnum} = 'Hand-based';
        case 'ego'
            model_titles{modelnum} = 'Egocentric';
        case 'cyl'
            model_titles{modelnum} = 'Cylindrical ego';
        case 'joint'
            model_titles{modelnum} = 'Joint-based';
        case 'markers'
            model_titles{modelnum} = 'Marker-based';
        end
    end
    model_titles{end} = 'Real';
    % model_alias = {'Extrinsic','Egocentric','Muscle','Real'};
    space_alias = {'PM','DL'};
    % for each model and space
    % num_neurons = length(td_test{1}(1).S1_unit_guide);
    fr_total = [];
    ax_pos = zeros(1,num_models);
    ax_vel = zeros(1,num_models);
    for modelnum = 1:num_models
        for spacenum = 1:2
            % position plots
            td = td_test{spacenum};
            fr = get_vars(td,{model_names{modelnum},neuron_idx});

            % position plots
            poss = cat(1,td.pos);
            ax_pos(modelnum) = subplot(num_models,3,(modelnum-1)*3+1);
            % tesselateFR(poss(:,1),poss(:,2),fr);
            plotSmoothFR(poss(:,1),poss(:,2),fr);
            hold on

            % velocity plots
            vels = cat(1,td.vel);
            % Plot tuning surface
            ax_vel(modelnum,spacenum) = subplot(num_models,3,(modelnum-1)*3+spacenum+1);
            % tesselateFR(vels(:,1),vels(:,2),fr);
            plotSmoothFR(vels(:,1),vels(:,2),fr,[-40 40]);

            % collect firing rates to get color limits later
            fr_total = [fr_total; fr];

            % Plot preferred direction in white (if available)
            if ~isempty(pdTables)
                pd = pdTables{spacenum,modelnum}.velPD(neuron_idx);
                pdCI = pdTables{spacenum,modelnum}.velPDCI(neuron_idx,:);
                hold on
                plot([0 20*cos(pd)],[0 20*sin(pd)],'-r','linewidth',2)
                % plot([0 20*cos(pdCI(1)) 20*cos(pd) 20*cos(pdCI(2)) 0],[0 20*sin(pdCI(1)) 20*sin(pd) 20*sin(pdCI(2)) 0],'r','linewidth',1.74)
                patch([0 20*cos(pdCI(1)) 20*cos(pd) 20*cos(pdCI(2))],[0 20*sin(pdCI(1)) 20*sin(pd) 20*sin(pdCI(2))],'r','edgecolor','none','facealpha',0.6)
            end
        end
    end
    linkaxes(ax_pos,'xy')
    linkaxes(ax_vel,'xy')

    clims = prctile(fr_total,[5 95]);
    for modelnum = 1:num_models
        subplot(num_models,3,(modelnum-1)*3+1)
        caxis(clims)
        axis tight
        xlabel(sprintf('%s Pos tuning',model_titles{modelnum}))
        for spacenum = 1:2
            subplot(num_models,3,(modelnum-1)*3+spacenum+1)
            caxis(clims)
            xlim([-40 40])
            ylim([-40 40])
            xlabel(sprintf('%s Vel tuning, Workspace: %s',model_titles{modelnum},space_alias{spacenum}))

        end
    end

    subplot(num_models,3,2)
    label = ['Neuron ' num2str(td_test{1}(1).S1_unit_guide(neuron_idx,:))];
    title(label)

    colormap(cm_viridis)
