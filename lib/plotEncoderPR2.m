function plotEncoderPR2(encoderResults,x_model,y_model)
% makes scatter plot of pR2, given encoder results and model names

    % Set up variables
    % aliases
    switch(x_model)
    case 'ext'
        x_model_alias = 'Hand';
    case 'ego'
        x_model_alias = 'Egocentric';
    case 'musc'
        x_model_alias = 'Muscle';
    case 'cyl'
        x_model_alias = 'Cylindrical Hand';
    case 'joint'
        x_model_alias = 'Joint';
    case 'markers'
        x_model_alias = 'Marker';
    end
    switch(y_model)
    case 'ext'
        y_model_alias = 'Hand';
    case 'ego'
        y_model_alias = 'Egocentric';
    case 'musc'
        y_model_alias = 'Muscle';
    case 'cyl'
        y_model_alias = 'Cylindrical Hand';
    case 'joint'
        y_model_alias = 'Joint';
    case 'markers'
        y_model_alias = 'Marker';
    end
    
    avgEval = neuronAverage(encoderResults.crossEval,contains(encoderResults.crossEval.Properties.VariableDescriptions,'meta'));
    av_pR2_x = avgEval.(sprintf('glm_%s_model_eval',x_model));
    av_pR2_y = avgEval.(sprintf('glm_%s_model_eval',y_model));
    
    % get stats on pR2 diff between musc model and markers model
    alpha = 0.05;
    diffstat = encoderResults.crossEval.(sprintf('glm_%s_model_eval',y_model))-encoderResults.crossEval.(sprintf('glm_%s_model_eval',x_model));
    correction = 1/(encoderResults.params.num_folds*encoderResults.params.num_repeats) + 1/(encoderResults.params.num_folds-1);
    alphaup = 1-alpha/2;
    alphalow = alpha/2;
    
    dpR2CI = zeros(height(avgEval),2);
    for i = 1:height(avgEval)
        sigID = avgEval.signalID(i,:);
        idx = getNTidx(encoderResults.crossEval,'signalID',sigID);
        mudiff = mean(diffstat(idx));
        vardiff = var(diffstat(idx));
        upp = tinv(alphaup,encoderResults.params.num_folds*encoderResults.params.num_repeats-1);
        low = tinv(alphalow,encoderResults.params.num_folds*encoderResults.params.num_repeats-1);
    
        dpR2CI(i,1) = mudiff + low * sqrt(correction*vardiff);
        dpR2CI(i,2) = mudiff + upp * sqrt(correction*vardiff);
    end
    
    % classify neurons
    y_neurons = dpR2CI(:,1)>0;
    x_neurons = dpR2CI(:,2)<0;
    
    % Scatterplot of pR2 plotted against each other
    plot([-1 1],[-1 1],'k--','linewidth',2)
    hold on
    plot([0 0],[-1 1],'k-','linewidth',2)
    plot([-1 1],[0 0],'k-','linewidth',2)
    color = zeros(height(avgEval),3);
    color(y_neurons,:) = repmat(getModelColors(y_model),sum(y_neurons),1);
    color(x_neurons,:) = repmat(getModelColors(x_model),sum(x_neurons),1);
    scatter(av_pR2_x(:),av_pR2_y(:),50,color,'filled')
    set(gca,'box','off','tickdir','out','xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
    axis square
    xlabel(sprintf('%s-based pR2',x_model_alias))
    ylabel(sprintf('%s-based pR2',y_model_alias))
    
    % make CI plot
    % figure
    % for i = 1:height(avgEval)
    %     if y_neurons(i)
    %         color = getModelColors(y_model);
    %     elseif x_neurons(i)
    %         color = getModelColors(x_model);
    %     else
    %         color = [0 0 0];
    %     end
    %     plot(dpR2CI(i,:),[i i],'-','color',color,'linewidth',2);
    %     hold on
    %     plot(dpR2CI(i,:),[i i],'.','color',color,'markersize',30);
    % end
    % plot([0 0],[0 1+height(avgEval)],'k-','linewidth',2);
    % axis ij
    % axis([-0.15 0.15 -inf inf])
    % set(gca,'box','off','tickdir','out','ylim',[0 1+height(avgEval)],'ytick',[1 height(avgEval)])

