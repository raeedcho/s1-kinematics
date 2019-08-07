%% Plot margins
    S1_margin = sepResults.trial_table.S1_FR_self_margin;
    ext_margin = sepResults.trial_table.ext_self_margin;
    handelbow_margin = sepResults.trial_table.handelbow_self_margin;
    ext_actpasbaseline_margin = sepResults.trial_table.ext_actpasbaseline_self_margin;

    ext_input_margin = sepResults.trial_table.ext_input_margin;
    handelbow_input_margin = sepResults.trial_table.handelbow_input_margin;
    ext_actpasbaseline_input_margin = sepResults.trial_table.ext_actpasbaseline_input_margin;
    
    % S1_margin = sepResults.trial_table.S1_FR_true_margin;
    % ext_margin = sepResults.trial_table.ext_true_margin;
    % handelbow_margin = sepResults.trial_table.handelbow_true_margin;
    % ext_actpasbaseline_margin = sepResults.trial_table.ext_actpasbaseline_true_margin;
    isPassive = sepResults.trial_table.isPassive;

    figure('defaultaxesfontsize',18)
    subplot(1,3,1)
    % scatter(S1_margin(~isPassive),handelbow_margin(~isPassive),[],'k','filled')
    scatter(S1_margin(~isPassive),handelbow_input_margin(~isPassive),[],'k','filled')
    hold on
    % scatter(S1_margin(isPassive),handelbow_margin(isPassive),[],'r','filled')
    scatter(S1_margin(isPassive),handelbow_input_margin(isPassive),[],'r','filled')
    plot(xlim,[0 0],'-k')
    plot([0 0],ylim,'-k')
    plot(xlim,xlim,'--k')
    set(gca,'box','off','tickdir','out','xlim',xlim,'ylim',xlim)
    % axis equal
    title('Hand/Elbow vs. S1 Margin comparison')
    ylabel('Hand/Elbow model margin')
    xlabel('S1 FR margin')

    subplot(1,3,2)
    scatter(S1_margin(~isPassive),ext_input_margin(~isPassive),[],'k','filled')
    % scatter(S1_margin(~isPassive),ext_margin(~isPassive),[],'k','filled')
    hold on
    scatter(S1_margin(isPassive),ext_input_margin(isPassive),[],'r','filled')
    % scatter(S1_margin(isPassive),ext_margin(isPassive),[],'r','filled')
    plot(xlim,[0 0],'-k')
    plot([0 0],ylim,'-k')
    plot(xlim,xlim,'--k')
    set(gca,'box','off','tickdir','out','xlim',xlim,'ylim',xlim)
    % axis equal
    title('Extrinsic vs. S1 Margin comparison')
    ylabel('Extrinsic model margin')
    xlabel('S1 FR margin')

    subplot(1,3,3)
    scatter(S1_margin(~isPassive),ext_actpasbaseline_input_margin(~isPassive),[],'k','filled')
    % scatter(S1_margin(~isPassive),ext_actpasbaseline_margin(~isPassive),[],'k','filled')
    hold on
    scatter(S1_margin(isPassive),ext_actpasbaseline_input_margin(isPassive),[],'r','filled')
    % scatter(S1_margin(isPassive),ext_actpasbaseline_margin(isPassive),[],'r','filled')
    plot(xlim,[0 0],'-k')
    plot([0 0],ylim,'-k')
    plot(xlim,xlim,'--k')
    set(gca,'box','off','tickdir','out','xlim',xlim,'ylim',xlim)
    % axis equal
    title('Hand+Act/Pas baseline param vs. S1 Margin comparison')
    ylabel('Hand+baseline param model margin')
    xlabel('S1 FR margin')
