function plotMWTuningCloud(tuningHull,neuron_idx)
    signalIDs = unique(tuningHull(:,{'monkey','date','signalID'}));

    plot([0 0],0.05*[-1 1],'-k','linewidth',2)
    hold on
    plot(0.05*[-1 1],[0 0],'-k','linewidth',2)

    for spacenum = 1:2
        ID = tuningHull(:,{'monkey','date','signalID'});
        unit_idx = ismember(ID,signalIDs(neuron_idx,:));
        space_idx = tuningHull.spaceNum == spacenum;
        tuningHull_unit = tuningHull(unit_idx & space_idx,:);

        % skip if real neuron is not tuned
        real_hull = tuningHull_unit.S1_FR_velWeight{1};
        if inpolygon(0,0,real_hull(:,1),real_hull(:,2))
            continue
        end

        % plot patched cloud for each neuron
        patch(real_hull(:,1),real_hull(:,2),'k','facealpha',0.25*spacenum)

        % plot modeled patches
        ext_hull = tuningHull_unit.glm_ext_model_velWeight{1};
        patch(ext_hull(:,1),ext_hull(:,2),'r','facealpha',0.25*spacenum)
        ego_hull = tuningHull_unit.glm_ego_model_velWeight{1};
        patch(ego_hull(:,1),ego_hull(:,2),'g','facealpha',0.25*spacenum)
        musc_hull = tuningHull_unit.glm_musc_model_velWeight{1};
        patch(musc_hull(:,1),musc_hull(:,2),'b','facealpha',0.25*spacenum)
    end

    if isnumeric(signalIDs.signalID(neuron_idx))
        label = ['Neuron ' num2str(signalIDs.signalID(neuron_idx))];
    else
        label = ['Neuron ' signalIDs.signalID(neuron_idx)];
    end
    
    title(label)
    set(gca,'box','off','tickdir','out')

    axis equal
