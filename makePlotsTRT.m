%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes plots from results given by analyzeTRT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    unit_idx = 1;
    plotSpikesOnHandle(td_ex,struct('unit_idx',unit_idx,'spikespec','b.','spikesize',10));
    % plot of same muscle movement given different Jacobians?

    % Switch to classical PDs?
    % 1c - example directional rasters and tuning curves?

%% Figure 2 - Analysis block diagrams
    % 2a - Block diagram of three different models

    % 2b - Breaking up the data into training and testing sets

    % 2c - Example neural predictions for each model

%% Figure 3 - Comparison of modeled tuning curves
    % 3a - a few example flat tuning curves (models on top of each other) DL under PM

    % 3b - DNA plot for each model/real tuning for each monkey

%% Figure 4 - Summaries
    % 3c - Summary over neurons: PD shift plots for each model and each monkey
    % 3d - Summary over monkeys/sessions: Bar chart of variance explained for each model by 1:1 line

%%%%%%%%%%%%%%%%% Supplementary Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary A - pR2s against each other for muscle vs endpoint

%% Supplementary B - Figure showing why one monkey didn't have as large a change in tuning

%% Set up plotting variables
    model_aliases = {'ext','ego','musc','markers'};
    num_models = length(model_aliases)+1;
    model_titles = cell(num_models-1,1);
    model_colors = zeros(num_models-1,3);
    for modelnum = 1:num_models-1
        switch model_aliases{modelnum}
        case 'musc'
            model_colors(modelnum,:) = [0, 174, 239]/255;
            model_titles{modelnum} = 'Muscle-based';
        case 'ext'
            model_colors(modelnum,:) = [247, 148, 30]/255;
            model_titles{modelnum} = 'Hand-based';
        case 'ego'
            model_colors(modelnum,:) = [105, 189, 69]/255;
            model_titles{modelnum} = 'Egocentric';
        case 'cyl'
            model_colors(modelnum,:) = [113, 191, 110]/255;
            model_titles{modelnum} = 'Cylindrical ego';
        case 'joint'
            model_colors(modelnum,:) = [38, 34, 98]/255;
            model_titles{modelnum} = 'Joint-based';
        case 'markers'
            model_colors(modelnum,:) = [193, 25, 47]/255;
            model_titles{modelnum} = 'Hand/Elbow-based';
        end
    end

    % colors for pm, dl conditions
    cond_colors = [0.6,0.5,0.7;...
        1,0,0];

%% Plot PD shifts
    % get shifts from weights
    shift_tables = cell(1,num_models);
    mean_shifts = cell(1,num_models);
    for modelnum = 1:num_models
        % select tables for each space
        [~,pm_tuningTable] = getNTidx(encoderResults.crossTuning,'spaceNum',1);
        [~,dl_tuningTable] = getNTidx(encoderResults.crossTuning,'spaceNum',2);

        % compose shift table for this model/bootstrap sample
        key_cols = contains(pm_tuningTable.Properties.VariableDescriptions,'meta');
        shift_tables{modelnum} = pm_tuningTable(:,key_cols);

        % remove spaceNum from columns
        shift_tables{modelnum}.spaceNum = [];

        % get PDs from pm and dl
        pm_PDs = pm_tuningTable.([encoderResults.params.model_names{modelnum} '_velPD']);
        dl_PDs = dl_tuningTable.([encoderResults.params.model_names{modelnum} '_velPD']);
        dPDs = minusPi2Pi(dl_PDs-pm_PDs);

        tab_append = table(dPDs,'VariableNames',{'velPD'});
        tab_append.Properties.VariableDescriptions = {'circular'};
        shift_tables{modelnum} = [shift_tables{modelnum} tab_append];

        mean_shifts{modelnum} = neuronAverage(shift_tables{modelnum},contains(shift_tables{modelnum}.Properties.VariableDescriptions,'meta'));
    end

    figure
    for modelnum = 1:num_models-1
        [~,real_shifts] = getNTidx(mean_shifts{end},'signalID',encoderResults.tunedNeurons);
        [~,model_shifts] = getNTidx(mean_shifts{modelnum},'signalID',encoderResults.tunedNeurons);

        subplot(1,num_models-1,modelnum)
        plot([-180 180],[0 0],'-k','linewidth',2)
        hold on
        plot([0 0],[-180 180],'-k','linewidth',2)
        plot([-180 180],[-180 180],'--k','linewidth',2)
        axis equal
        set(gca,'box','off','tickdir','out','xtick',[-180 180],'ytick',[-180 180],'xlim',[-180 180],'ylim',[-180 180])
        scatter(180/pi*real_shifts.velPD,180/pi*model_shifts.velPD,50,model_colors(modelnum,:),'filled')

        % labels
        xlabel 'Actual PD Shift'
        ylabel 'Modeled PD Shift'
        title(sprintf('%s model PD shift vs Actual PD shift',model_titles{modelnum}))
    end

%% Plot out tuning curves
    % compare PM and DL tuning for each model
    for modelnum = 1:num_models
        figure
        title(encoderResults.params.model_names{modelnum})
        compareTuning(encoderResults.tuning_curves(:,modelnum),encoderResults.pdTables(:,modelnum),struct('which_units',find(encoderResults.isTuned),'cond_colors',cond_colors))
        % compareTuning(tuning_curves(:,modelnum),pdTables(:,modelnum),struct('which_units',find(encoderResults.isTuned),'cond_colors',cond_colors,'maxFR',1))
    end

%% Make histogram plots of PD changes
    figure
    subplot(num_models,1,1)
    h = histogram(gca,mean_shifts{end}.velPD*180/pi,'BinWidth',10,'DisplayStyle','stair');
    set(h,'edgecolor','k')
    set(gca,'box','off','tickdir','out','xlim',[-180 180],'xtick',[-180 0 180],'ylim',[0 20],'ytick',[0 20])
    for modelnum = 1:num_models-1
        subplot(num_models,1,modelnum+1)
        h = histogram(gca,mean_shifts{modelnum}.velPD*180/pi,'BinWidth',10,'DisplayStyle','stair');
        set(h,'edgecolor',model_colors(modelnum,:))
        set(gca,'box','off','tickdir','out','xlim',[-180 180],'xtick',[-180 0 180],'ylim',[0 20],'ytick',[0 20])
    end
    xlabel 'Change in Preferred Direction'

%% Calculate mean error on shifts
    err = zeros(100,num_models-1);
    for modelnum = 1:num_models-1
        [~,real_shifts] = getNTidx(shift_tables{end},'signalID',encoderResults.tunedNeurons);
        [~,model_shifts] = getNTidx(shift_tables{modelnum},'signalID',encoderResults.tunedNeurons);
        err_arr = model_shifts.velPD-real_shifts.velPD;
        for i = 1:100
            err_idx = 1:length(encoderResults.tunedNeurons);
            err_idx = err_idx + (i-1)*length(encoderResults.tunedNeurons);
            % use a 1-cos style error because of circular data
            % This value will range between 0 and 2
            err(i,modelnum) = mean(1-cos(err_arr(err_idx)));
        end
    end

    % plot histograms
    % figure
    % for modelnum = 1:num_models-1
    %     h = histogram(gca,err(:,modelnum),'DisplayStyle','bar');
    %     set(h,'edgecolor',model_colors(modelnum,:),'facecolor',model_colors(modelnum,:))
    %     hold on
    % end

    % compute statistics
    models_to_compare = [find(contains(model_aliases,'musc')) find(contains(model_aliases,'marker'))];
    alpha = 0.05/2; % bonferroni correction for multiple comparisons...?
    diffstat = err(:,models_to_compare(1))-err(:,models_to_compare(2)); % musc - ext
    mudiff = mean(diffstat);
    vardiff = var(diffstat);
    correction = 1/100 + 1/4;
    % alphaup = 1-alpha;
    upp = tinv(0.975,99);
    low = tinv(0.025,99);
    errCIhi = mudiff + upp * sqrt(correction*vardiff)
    errCIlo = mudiff + low * sqrt(correction*vardiff)

    diffstat = err(:,models_to_compare(1))-err(:,models_to_compare(2)); % musc - ego
    mudiff = mean(diffstat);
    vardiff = var(diffstat);
    correction = 1/100 + 1/4;
    alphaup = 1-alpha;
    upp = tinv(alphaup,99);
    errCIhigh_ego = mudiff + upp * sqrt(correction*vardiff)

    % plot errors
    figure
    for modelnum = 1:num_models-1
        scatter(err(:,modelnum),repmat(modelnum/10,size(err,1),1),50,model_colors(modelnum,:),'filled')
        hold on
        plot(mean(err(:,modelnum)),modelnum/10,'k.','linewidth',3,'markersize',40)
    end
    set(gca,'tickdir','out','box','off','ytick',(1:(num_models-1))/10,'yticklabel',model_titles,'xtick',[0 1])
    axis equal
    axis ij
    xlim([0 1])
    ylim([0 num_models/10])
    xlabel('Cosine error of model')

%% Plot pR2s against each other
    % setup
    x_model = 'ext';
    y_model = 'markers';

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

    % make plot of pR2 of muscle against markers
    % figure
    % plot([-1 1],[-1 1],'k--','linewidth',2)
    % hold on
    % plot([0 0],[-1 1],'k-','linewidth',2)
    % plot([-1 1],[0 0],'k-','linewidth',2)
    % for i = 1:height(avgEval)
    %     sigID = avgEval.signalID(i,:);
    %     idx = getNTidx(crossEval,'signalID',sigID);
    %     if y_neurons(i)
    %         color = model_colors(contains(model_names,y_model),:);
    %     elseif x_neurons(i)
    %         color = model_colors(contains(model_names,x_model),:);
    %     else
    %         color = [0 0 0];
    %     end
    %     scatter(crossEval.(sprintf('glm_%s_model_eval',x_model))(idx),crossEval.(sprintf('glm_%s_model_eval',y_model))(idx),25,color,'filled')
    % end
    % clear i idx color sigID;
    % scatter(av_pR2_x(:),av_pR2_y(:),50,'r','filled')
    % set(gca,'box','off','tickdir','out','xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
    % axis square
    % xlabel(sprintf('%s-based pR2',x_model_alias))
    % ylabel(sprintf('%s-based pR2',y_model_alias))

    % Scatterplot of pR2 plotted against each other
    figure
    plot([-1 1],[-1 1],'k--','linewidth',2)
    hold on
    plot([0 0],[-1 1],'k-','linewidth',2)
    plot([-1 1],[0 0],'k-','linewidth',2)
    color = zeros(height(avgEval),3);
    color(y_neurons,:) = repmat(model_colors(contains(model_aliases,y_model),:),sum(y_neurons),1);
    color(x_neurons,:) = repmat(model_colors(contains(model_aliases,x_model),:),sum(x_neurons),1);
    scatter(av_pR2_x(:),av_pR2_y(:),50,color,'filled')
    set(gca,'box','off','tickdir','out','xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
    axis square
    xlabel(sprintf('%s-based pR2',x_model_alias))
    ylabel(sprintf('%s-based pR2',y_model_alias))

    % make plot
    figure
    for i = 1:height(avgEval)
        if y_neurons(i)
            color = model_colors(contains(model_aliases,y_model),:);
        elseif x_neurons(i)
            color = model_colors(contains(model_aliases,x_model),:);
        else
            color = [0 0 0];
        end
        plot(dpR2CI(i,:),[i i],'-','color',color,'linewidth',2);
        hold on
        plot(dpR2CI(i,:),[i i],'.','color',color,'markersize',30);
    end
    plot([0 0],[0 1+height(avgEval)],'k-','linewidth',2);
    axis ij
    axis([-0.15 0.15 -inf inf])
    set(gca,'box','off','tickdir','out','ylim',[0 1+height(avgEval)],'ytick',[1 height(avgEval)])

%% Plot error on all monkeys (requires already saved errors)
    % get mean and standard error of errors
    err = [han_err chips_err lando_err];
    mean_err = mean(err);
    var_err = var(err);
    correction = 1/100 + 1/4;
    std_err_err = sqrt(correction*var_err);

    num_monks = 1;
    % num_models = 3;
    num_cols = num_monks*(num_models-1); % this should be the number columns in err now
    model_colors_rep = repmat(model_colors,num_monks,1);
    % model_x = [1.5 2 2.5 4.5 5 5.5 7.5 8 8.5]/10;
    monk_x = (2:3:((num_monks-1)*3+2))/10;
    template_x = linspace(-0.5,0.5,num_models-1)/10;
    model_spacing = mode(diff(template_x));
    model_x = [];
    for i = 1:length(monk_x)
        model_x = [model_x template_x+monk_x(i)];
    end
    figure
    for colnum = 1:num_cols
        % scatter(repmat(colnum/10,size(err,1),1),err(:,colnum),50,model_colors_rep(colnum,:),'filled')
        % hold on
        % plot(colnum/10,mean(err(:,colnum)),'.','color',model_colors_rep(colnum,:),'linewidth',3,'markersize',40)
        bar(model_x(colnum),mean(err(:,colnum)),model_spacing,'facecolor',model_colors_rep(colnum,:))
        hold on
        % errorbar(colnum/10,mean_err(:,colnum),std_err_err(:,colnum),'color',model_colors_rep(colnum,:),'linewidth',3)
    end
    errorbar(model_x,mean_err,std_err_err,'k.','linewidth',3)
    set(gca,'tickdir','out','box','off','xtick',monk_x,'xticklabel',{'Monkey H','Monkey C','Monkey L'},'ytick',[0 1])
    axis equal
    ylim([0 1])
    xlim([0 1])
    ylabel('Error of model')

%% Tuning curve covariances
    tuning_covar = zeros(2,num_models,height(encoderResults.tuning_curves{1,1}));
    for neuron_idx = 1:height(encoderResults.tuning_curves{1,1})
        for spacenum = 1:2
            tuning_curve_mat = zeros(8,num_models);
            for modelnum = 1:num_models
                tuning_curve_mat(:,modelnum) = encoderResults.tuning_curves{spacenum,modelnum}(neuron_idx,:).velCurve';
            end
            covar_mat = cov(tuning_curve_mat);
            true_tuning_idx = contains(model_names,'S1');
            tuning_covar(spacenum,:,neuron_idx) = covar_mat(:,true_tuning_idx)/covar_mat(true_tuning_idx,true_tuning_idx);
        end
    end

%% Example predictions
    for neuron_idx = 23% 1:height(avgEval)
        h = figure;
        temp_vel = cat(1,encoderResults.td_tuning{2}.vel);
        temp_spikes = get_vars(encoderResults.td_tuning{2},{'S1_FR',neuron_idx});
        temp_pred_ext = get_vars(encoderResults.td_tuning{2},{'glm_ext_model',neuron_idx});
        temp_pred_musc = get_vars(encoderResults.td_tuning{2},{'glm_musc_model',neuron_idx});

        clf
        ax1 = subplot(2,1,1);
        plot(temp_vel(:,1),'b','linewidth',2)
        hold on
        plot(temp_vel(:,2),'g','linewidth',2)
        set(gca,'box','off','tickdir','out')

        ax2 = subplot(2,1,2);
        plot(temp_spikes,'k','linewidth',2)
        hold on
        plot(temp_pred_ext,'color',model_colors(contains(model_names,'ext'),:),'linewidth',2)
        plot(temp_pred_musc,'color',model_colors(contains(model_names,'musc'),:),'linewidth',2)
        title(sprintf('Hand-based pR^2: %f, Musc-based pR^2: %f',av_pR2_ext(neuron_idx),av_pR2_musc(neuron_idx)))
        set(gca,'box','off','tickdir','out')

        linkaxes([ax1 ax2],'x')
        waitfor(h)
    end
    clearvars neuron_idx temp_* ax1 ax2

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

%% Plot tuning weight clouds
    % tuningHull = getTuningHull(results.tuningTable);
    % loop over each unit in one workspace
    % n_rows = ceil(sqrt(height(signalIDs)+1));
    % cloud_fig = figure;
    % surf_fig = figure;
    % neuron_idx = getNTidx(pdTables{1,4},'signalID',[95 2]);
    for neuron_idx = 1:length(encoderResults.isTuned)
        % close_fig = figure;

        % figure(cloud_fig)
        % clf
        % plotMWTuningCloud(tuningHull,neuron_idx)

        surf_fig = figure;
        clf
        plotMWTuningSurfaces(encoderResults.td_tuning,pdTables,neuron_idx,model_aliases)

        waitfor(surf_fig)
    end

%% Decoder crap
    %% load data

    %% make decoder performance scatter plots
    barwidth = 0.4;
    figure
    bar([1 3 2 4]*0.20,mean([hand_decoder_vaf neur_decoder_vaf]),barwidth,'facecolor',[0.5 0.5 0.5],'edgecolor','none')
    hold on
    plot([1 2]*0.20,[hand_decoder_vaf(:,1) neur_decoder_vaf(:,1)]','-k','linewidth',1)
    plot([3 4]*0.20,[hand_decoder_vaf(:,2) neur_decoder_vaf(:,2)]','-k','linewidth',1)
    plot([1 3 2 4]*0.20, [hand_decoder_vaf neur_decoder_vaf]','.k','markersize',30)
    xtick = (1:4)*0.20;
    xticklabel = {'Hand-only pos','Hand+Neuron pos','Hand vel','Hand+Neuron vel'};
    set(gca,'box','off','tickdir','out','xtick',xtick,...
        'xticklabel',xticklabel,...
        'xlim',[0 3],'ylim',[0 1])
    ylabel 'Fraction VAF'
    xlabel 'Model'
    title 'Decoding performance'

