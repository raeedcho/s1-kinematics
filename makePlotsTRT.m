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
    figure('defaultaxesfontsize',18)
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
    datadir = '/home/raeed/Projects/limblab/data-td/MultiWorkspace/Results/Encoding';
    filename = {'Han_20171101_TRT_encodingResults_run20180809.mat','Chips_20170915_TRT_encodingResults_run20180809.mat','Lando_20170802_encodingResults_run20180809.mat'};
    num_monks = length(filename);
    err = cell(num_monks,1);
    hyp = cell(num_monks,1);
    p_val = cell(num_monks,1);

    model_aliases = {'ext','ego','musc','markers'};
    num_models = length(model_aliases)+1;
    model_titles = getModelTitles(model_aliases);
    model_colors = getModelColors(model_aliases);
    % colors for pm, dl conditions
    cond_colors = [0.6,0.5,0.7;...
        1,0,0];

%% Loop over all monkeys for encoder figures and errors
    for monkeynum = 1:num_monks
        clear encoderResults

        % load data
        load(fullfile(datadir,filename{monkeynum}))
    
        %% Plot PD shifts
            % get shifts from weights
            shift_tables = calculatePDShiftTables(encoderResults);
        
            mean_shifts = cell(num_models,1);
            for modelnum = 1:num_models
                mean_shifts{modelnum} = neuronAverage(shift_tables{modelnum},contains(shift_tables{modelnum}.Properties.VariableDescriptions,'meta'));
            end

            figure('defaultaxesfontsize',18)
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
                title({sprintf('%s model PD shift vs Actual PD shift',model_titles{modelnum});filename{monkeynum}},'interpreter','none')
            end
        
        %% Plot out tuning curves
            % compare PM and DL tuning for each model
            % for modelnum = 1:num_models
            %     figure('defaultaxesfontsize',18)
            %     title(encoderResults.params.model_names{modelnum})
            %     compareTuning(encoderResults.tuning_curves(:,modelnum),encoderResults.pdTables(:,modelnum),struct('which_units',find(encoderResults.isTuned),'cond_colors',cond_colors))
            %     % compareTuning(tuning_curves(:,modelnum),pdTables(:,modelnum),struct('which_units',find(encoderResults.isTuned),'cond_colors',cond_colors,'maxFR',1))
            % end
        
        %% Make histogram plots of PD changes
            figure('defaultaxesfontsize',18)
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
            title(filename(monkeynum),'interpreter','none')
        
        %% Calculate mean error on shifts
            err{monkeynum} = calculateEncoderPDShiftErr(encoderResults);

            % % plot errors
            % figure('defaultaxesfontsize',18)
            % for modelnum = 1:num_models-1
            %     scatter(err{:,modelnum},repmat(modelnum/10,size(err,1),1),50,model_colors(modelnum,:),'filled')
            %     hold on
            %     plot(mean(err{:,modelnum}),modelnum/10,'k.','linewidth',3,'markersize',40)
            % end
            % set(gca,'tickdir','out','box','off','ytick',(1:(num_models-1))/10,'yticklabel',model_titles,'xtick',[0 1])
            % axis equal
            % axis ij
            % xlim([0 1])
            % ylim([0 num_models/10])
            % xlabel('Cosine error of model')

            modelcompare = {'ext','ego';...
                'ext','musc';...
                'ext','markers';...
                'ego','musc';...
                'ego','markers';
                'musc','markers'};
            tails = {'both';'right';'right';'right';'right';'right'};
            [hyp{monkeynum},p_val{monkeynum}] = stattestPDShiftErr(err{monkeynum},modelcompare,tails,encoderResults.params.num_repeats,encoderResults.params.num_folds);
        
        %% Plot pR2s against each other
            % setup
            figure('defaultaxesfontsize',18)
            subplot(2,3,1)
            plotEncoderPR2(encoderResults,'ext','ego')
            subplot(2,3,2)
            plotEncoderPR2(encoderResults,'ext','markers')
            title(filename{monkeynum},'interpreter','none') % centered title
            subplot(2,3,3)
            plotEncoderPR2(encoderResults,'ego','markers')
            subplot(2,3,4)
            plotEncoderPR2(encoderResults,'musc','markers')
            subplot(2,3,5)
            plotEncoderPR2(encoderResults,'ext','musc')
            subplot(2,3,6)
            plotEncoderPR2(encoderResults,'ego','musc')
    end
    
%% Histogram of PD shift for all monkeys
    figure('defaultaxesfontsize',18)
    for monkeynum = 1:num_monks
        % load data
        load(fullfile(datadir,filename{monkeynum}))

        shift_tables = calculatePDShiftTables(encoderResults);
        mean_shifts = cell(num_models,1);
        for modelnum = 1:num_models
            mean_shifts{modelnum} = neuronAverage(shift_tables{modelnum},contains(shift_tables{modelnum}.Properties.VariableDescriptions,'meta'));
        end

        subplot(num_monks,1,monkeynum)
        h = histogram(gca,mean_shifts{end}.velPD*180/pi,'BinWidth',10,'DisplayStyle','stair');
        set(h,'edgecolor','k')
        set(gca,'box','off','tickdir','out','xlim',[-180 180],'xtick',[-180 0 180],'ylim',[0 20],'ytick',[0 20])
        title(filename(monkeynum),'interpreter','none')
    end
    xlabel 'Change in Preferred Direction'

%% Plot error on all monkeys
    num_monks = 3;
    correction = 1/100 + 1/4;
    models_to_plot = {'ext','ego','musc','markers'};
    % x coordinate of individual monkey bars
    monk_x = (2:3:((num_monks-1)*3+2))/10;
    % template for within monkey bars separation
    template_x = linspace(-0.5,0.5,length(models_to_plot))/10;
    model_spacing = mode(diff(template_x));

    % make plot
    figure('defaultaxesfontsize',18)
    for monkeynum = 1:num_monks
        model_idx = find(contains(err{monkeynum}.Properties.VariableNames,models_to_plot));
        mean_err = mean(err{monkeynum}{:,model_idx});
        var_err = var(err{monkeynum}{:,model_idx});
        std_err_err = sqrt(correction*var_err);

        for modelnum = 1:length(models_to_plot)
            xval = monk_x(monkeynum) + template_x(modelnum);
            bar(xval,mean_err(modelnum),model_spacing,'facecolor',model_colors(model_idx(modelnum),:),'edgecolor','none')
            hold on
            plot([xval xval],[mean_err(modelnum)-std_err_err(modelnum) mean_err(modelnum)+std_err_err(modelnum)],'k','linewidth',3)
        end
    end
    set(gca,'tickdir','out','box','off','xtick',monk_x,...
        'xticklabel',filename,'ytick',[0 1],'ticklabelinterpreter','none')
    % axis equal
    ylim([0 0.7])
    % xlim([0 1])
    ylabel('Error of model')

%% Plot pR2 of all monkeys
    num_monks = 3;
    correction = 1/100 + 1/4;
    models_to_plot = {'ego','ext','musc','markers'};
    % x coordinate of individual monkey bars
    monk_x = (2:3:((num_monks-1)*3+2))/10;
    % template for within monkey bars separation
    template_x = linspace(-0.5,0.5,length(models_to_plot))/10;
    model_spacing = mode(diff(template_x));

    % make plot
    figure('defaultaxesfontsize',18)
    for monkeynum = 1:num_monks
        % load data
        load(fullfile(datadir,filename{monkeynum}))

        avgEval = neuronAverage(encoderResults.crossEval,contains(encoderResults.crossEval.Properties.VariableDescriptions,'meta'));

        % model_idx = find(contains(err{monkeynum}.Properties.VariableNames,models_to_plot));
        % mean_err = mean(err{monkeynum}{:,model_idx});
        % var_err = var(err{monkeynum}{:,model_idx});
        % std_err_err = sqrt(correction*var_err);

        avg_pR2 = zeros(height(avgEval),length(models_to_plot));
        for modelnum = 1:length(models_to_plot)
            xval = monk_x(monkeynum) + template_x(modelnum);
            mean_pR2 = mean(avgEval.(sprintf('glm_%s_model_eval',models_to_plot{modelnum})));
            avg_pR2(:,modelnum) = avgEval.(sprintf('glm_%s_model_eval',models_to_plot{modelnum}));
            bar(xval,mean_pR2,model_spacing,'facecolor',getModelColors(models_to_plot{modelnum}),'edgecolor','none')
            hold on
            % plot([xval xval],[mean_err(modelnum)-std_err_err(modelnum) mean_err(modelnum)+std_err_err(modelnum)],'k','linewidth',3)
        end
        xval = repmat(monk_x(monkeynum)+template_x,length(avg_pR2),1);
        scatter(xval(:),avg_pR2(:),[],'k','filled')
        plot(xval',avg_pR2','-k','linewidth',1)
    end
    set(gca,'tickdir','out','box','off','xtick',monk_x,...
        'xticklabel',filename,'ytick',[0 0.25 0.5],'ticklabelinterpreter','none')
    % axis equal
    ylim([0 0.6])
    % xlim([0 1])
    ylabel('Model pseudo-R^2')

%% Get example tuning curves for all models
    for monkeynum = 1%:num_monks
        clear encoderResults

        % load data
        load(fullfile(datadir,filename{monkeynum}))

        %% Plot out tuning curves
            % compare PM and DL tuning for each model
            for modelnum = 1:num_models
                figure('defaultaxesfontsize',18)
                % figure
                compareTuning(encoderResults.tuning_curves(:,modelnum),encoderResults.pdTables(:,modelnum),struct('which_units',find(encoderResults.isTuned),'cond_colors',cond_colors))
                % compareTuning(encoderResults.tuning_curves(:,modelnum),encoderResults.pdTables(:,modelnum),struct('which_units',find(encoderResults.isTuned),'cond_colors',cond_colors,'maxFR',1))
                title(encoderResults.params.model_names{modelnum},'interpreter','none')
            end
    end

%% Decoder crap
    datadir = '/home/raeed/Projects/limblab/data-td/FullWS/Results/Decoding';
    filename = {'Han_20160325_RWhold_decodingResults_run20180813.mat','Chips_20151211_RW_decodingResults_run20180813.mat'};
    num_monks = length(filename);
    barwidth = 0.4;
    markerstyle = '^o';
    barfig = figure('defaultaxesfontsize',18);
    scatfig = figure('defaultaxesfontsize',18);
    for monkeynum = 1:num_monks
        %% load data
        load(fullfile(datadir,filename{monkeynum}))

        %% make decoder performance scatter plots
        x_offset = monkeynum-1;
        % first monkey
        figure(barfig)
        bar([1 3 2 4]*0.20+x_offset,mean([decoderResults.hand_decoder_vaf decoderResults.neur_decoder_vaf]),barwidth,'facecolor',[0.5 0.5 0.5],'edgecolor','none')
        hold on
        plot([1 2]*0.20+x_offset,[decoderResults.hand_decoder_vaf(:,1) decoderResults.neur_decoder_vaf(:,1)]','-k','linewidth',1)
        plot([3 4]*0.20+x_offset,[decoderResults.hand_decoder_vaf(:,2) decoderResults.neur_decoder_vaf(:,2)]','-k','linewidth',1)
        plot([1 3 2 4]*0.20+x_offset, [decoderResults.hand_decoder_vaf decoderResults.neur_decoder_vaf]','.k','markersize',30)

        figure(scatfig)
        subplot(1,2,1)
        scatter(decoderResults.hand_decoder_vaf(:,1),decoderResults.neur_decoder_vaf(:,1),markerstyle(monkeynum),'filled')
        hold on
        plot([0 1],[0 1],'--k','linewidth',2)

        subplot(1,2,2)
        scatter(decoderResults.hand_decoder_vaf(:,2),decoderResults.neur_decoder_vaf(:,2),markerstyle(monkeynum),'filled')
        hold on
        plot([0 1],[0 1],'--k','linewidth',2)
    end
    figure(barfig)
    xtick = [(1:4)*0.20 1+(1:4)*0.20];
    xticklabel = repmat({'Hand-only pos','Hand+Neuron pos','Hand vel','Hand+Neuron vel'},1,2);
    set(gca,'box','off','tickdir','out','xtick',xtick,...
        'xticklabel',xticklabel,...
        'xlim',[0 num_monks],'ylim',[0 1])
    ylabel 'Fraction VAF'
    xlabel 'Model'
    title([{'Decoding performance'};filename'],'interpreter','none')

    figure(scatfig)
    subplot(1,2,1)
    axis equal
    subplot(1,2,2)
    axis equal

%% Extra stuff/in progress...
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
        td_tuning = encoderResults.td_tuning{2};
        td_tuning = smoothSignals(td_tuning,struct('signals','S1_FR','kernel_SD',0.1));
        trial_to_plot = randperm(length(td_tuning),5);
        num_neurons = length(td_tuning(1).S1_unit_guide);
        for neuron_idx = 23
            h = figure('defaultaxesfontsize',18);
            trial_to_plot = randperm(length(td_tuning),5);
            temp_vel = get_vars(td_tuning(trial_to_plot),{'vel',1:2});
            temp_spikes = get_vars(td_tuning(trial_to_plot),{'S1_FR',neuron_idx});
            temp_pred_ext = get_vars(td_tuning(trial_to_plot),{'glm_ext_model',neuron_idx});
            temp_pred_ego = get_vars(td_tuning(trial_to_plot),{'glm_ego_model',neuron_idx});
            temp_pred_musc = get_vars(td_tuning(trial_to_plot),{'glm_musc_model',neuron_idx});
            temp_pred_markers = get_vars(td_tuning(trial_to_plot),{'glm_markers_model',neuron_idx});

            clf
%             ax1 = subplot(2,1,1);
%             plot(temp_vel(:,1),'b','linewidth',2)
%             hold on
%             plot(temp_vel(:,2),'g','linewidth',2)
%             set(gca,'box','off','tickdir','out')
%     
%             ax2 = subplot(2,1,2);
            plot(temp_spikes,'k','linewidth',2)
            hold on
            plot(temp_pred_ext,'color',model_colors(contains(model_aliases,'ext'),:),'linewidth',2)
            plot(temp_pred_ego,'color',model_colors(contains(model_aliases,'ego'),:),'linewidth',2)
            plot(temp_pred_musc,'color',model_colors(contains(model_aliases,'musc'),:),'linewidth',2)
%             plot(temp_pred_markers,'color',model_colors(contains(model_aliases,'markers'),:),'linewidth',2)
%             title(sprintf('Hand-based pR^2: %f, Musc-based pR^2: %f',av_pR2_ext(neuron_idx),av_pR2_musc(neuron_idx)))
            set(gca,'box','off','tickdir','out')
    
%             linkaxes([ax1 ax2],'x')
            waitfor(h)
        end
        clearvars neuron_idx temp_* ax1 ax2
    
    %% Plot handle positions
        if verbose
            figure('defaultaxesfontsize',18)
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
        % cloud_fig = figure('defaultaxesfontsize',18);
        % surf_fig = figure('defaultaxesfontsize',18);
        % neuron_idx = getNTidx(pdTables{1,4},'signalID',[95 2]);
        for neuron_idx = 1:length(encoderResults.isTuned)
            % close_fig = figure('defaultaxesfontsize',18);
    
            % figure(cloud_fig)
            % clf
            % plotMWTuningCloud(tuningHull,neuron_idx)
    
            surf_fig = figure('defaultaxesfontsize',18);
            clf
            plotMWTuningSurfaces(encoderResults.td_tuning,pdTables,neuron_idx,model_aliases)
    
            waitfor(surf_fig)
        end
