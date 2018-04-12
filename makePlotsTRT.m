%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes plots from results given by analyzeTRT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get results

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

    % 1c - example directional rasters and tuning curves?

%% Figure 2 - Analysis block diagrams
    % 2a - Block diagram of three different models

    % 2b - Breaking up the data into training and testing sets

    % 2c - Example neural predictions for each model
    % neuron 1 has decent pseudo R2
    for neuron_idx = 52%1:length(td_eval{2,1})
        h = figure;
        temp_vel = cat(1,results.td_test{2}.vel);
        temp_spikes = cat(1,results.td_test{2}.S1_FR);
        temp_pred_ext = cat(1,results.td_test{2}.(model_names{1}));
        temp_pred_musc = cat(1,results.td_test{2}.(model_names{3}));

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
    clearvars neuron_idx temp_* ax1 ax2

%% Figure 3 - Comparison of modeled tuning curves
    % 3a - a few example flat tuning curves (models on top of each other) DL under PM

    % 3b - DNA plot for each model/real tuning for each monkey

%% Figure 4 - Summaries
    % 3c - Summary over neurons: PD shift plots for each model and each monkey
    % 3d - Summary over monkeys/sessions: Bar chart of variance explained for each model by 1:1 line

%%%%%%%%%%%%%%%%% Supplementary Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Supplementary A - pR2s against each other for muscle vs endpoint

%% Supplementary B - Figure showing why one monkey didn't have as large a change in tuning

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

%% Plot example neuron model fit

%% Plot pR2s against each other
    % av_pR2_ext_dl = mean(td_ext_dl_eval,2);
    % av_pR2_musc_dl = mean(td_musc_dl_eval,2);
    % av_pR2_ext_pm = mean(td_ext_pm_eval,2);
    % av_pR2_musc_pm = mean(td_musc_pm_eval,2);
    avgEval = neuronAverage(crossEval,contains(crossEval.Properties.VariableDescriptions,'meta'));
    av_pR2_ext = avgEval.glm_ext_model_eval;
    av_pR2_musc = avgEval.glm_musc_model_eval;

    good_neurons = isTuned{4};
    figure
    scatter(av_pR2_ext(good_neurons),av_pR2_musc(good_neurons),50,'k','filled')
    hold on
    plot([-1 1],[-1 1],'k--','linewidth',2)
    plot([0 0],[-1 1],'k-','linewidth',2)
    plot([-1 1],[0 0],'k-','linewidth',2)
    set(gca,'box','off','tickdir','out','xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
    xlabel 'Hand-based pR2'
    ylabel 'Muscle-based pR2'

    % clean up
    clearvars av_pR2* good_neurons

%% Plot comparison of actual tuning curves with various modeled tuning curves
    % TODO: Redo this block so that PDs and tuning curves are calculated from whole trial_data w/ bootstrapping
    % Get PDs and tuning curves for the modeled and actual neurons
    tuning_curves = cell(2,4); % PM is first row, DL is second. Column order is Ext, Ego, Musc, Real
    pdTables = cell(2,4);
    tuning_params = cell(1,4);
    model_type = 'glm';
    model_names = [strcat(model_type,{'_ext','_ego','_musc'},'_model') {'S1_FR'}];
    isTuned = cell(1,4);
    tunedNeurons = cell(1,4);
    
    % get PDs
    pdConvertedTable = getPDsFromWeights(crossTuning);
    curveTable = getCurvesFromCrossval(crossTuning);

    % get PDs and tuning curves
    for modelnum = 1:4
        for spacenum = 1:2
            % move converted PD table into separated cells with table selections
            [~,temp] = getNTidx(pdConvertedTable,'spaceNum',spacenum);
            % select only columns corresponding to the key or the model in question
            key_cols = contains(temp.Properties.VariableDescriptions,'meta');
            model_cols = contains(temp.Properties.VariableNames,model_names{modelnum});
            temp = temp(:,key_cols | model_cols);
            % rename columns to get rid of model names
            temp.Properties.VariableNames = strrep(temp.Properties.VariableNames,[model_names{modelnum} '_'],'');
            pdTables{spacenum,modelnum} = temp;

            % tuning_curves{spacenum,modelnum} = getTuningCurves(results.td_test{spacenum},tuning_params{modelnum});
            [~,temp] = getNTidx(curveTable,'spaceNum',spacenum);
            % select only columns corresponding to the key or the model in question
            key_cols = contains(temp.Properties.VariableDescriptions,'meta');
            bins_cols = endsWith(temp.Properties.VariableNames,'bins');
            model_cols = contains(temp.Properties.VariableNames,model_names{modelnum});
            temp = temp(:,key_cols | bins_cols | model_cols);
            % rename columns to get rid of model names
            temp.Properties.VariableNames = strrep(temp.Properties.VariableNames,[model_names{modelnum} '_'],'');
            tuning_curves{spacenum,modelnum} = temp;
        end
    end

    % check whether each model is tuned
    tuningHull = getTuningHull(crossTuning);
    for modelnum = 1:4
        % set up isTuned cell
        isTuned{modelnum} = true(height(tuningHull)/2,1);
        tunedNeurons{modelnum} = [];
        for spacenum = 1:2
            % get only entries in given space
            [~,tuningHull_space] = getNTidx(tuningHull,'spaceNum',spacenum);
            for neuron_idx = 1:height(tuningHull_space)
                hull = tuningHull_space(neuron_idx,:).([model_names{modelnum} '_velWeight']){1};
                isTuned{modelnum}(neuron_idx) = isTuned{modelnum}(neuron_idx) & ~inpolygon(0,0,hull(:,1),hull(:,2));
            end
        end
        tunedNeurons{modelnum} = tuningHull_space.signalID(isTuned{modelnum},:);
    end

    % compare PM and DL tuning for each model
    for modelnum = 1:4
        % figure;compareTuning(tuning_curves(:,modelnum),pdTables(:,modelnum))
        figure;compareTuning(tuning_curves(:,modelnum),pdTables(:,modelnum),find(isTuned{4}))
    end

    % then compare PM and DL tuning for each model
    % reorder for color consistency..
    % for spacenum = 1:2
    %     figure;compareTuning(tuning_curves(spacenum,[3,1,2,4]),pdTables(spacenum,[3,1,2,4]),find(results.isTuned{4}))
    % end

%% Make iris and dna plots
    % use only "tuned" neurons
    num_models = 4;
    f1 = figure;
    f2 = figure;
    for modelnum = 1:num_models
        figure(f1)
        subplot(2,2,modelnum)
        irisPlot(pdTables{1,modelnum}(isTuned{4},:),pdTables{2,modelnum}(isTuned{4},:));
        title(model_names{modelnum})

        figure(f2)
        subplot(2,2,modelnum)
        dnaPlot(pdTables{1,modelnum}(isTuned{4},:),pdTables{2,modelnum}(isTuned{4},:));
        title(model_names{modelnum})
    end
    
%% Cross-validate models of neural data
    % Get crossval info
    [crossEval,crossTuning] = analyzeTRT(trial_data);

    % get shifts from weights
    shift_tables = cell(1,num_models);
    mean_shifts = cell(1,num_models);
    for modelnum = 1:num_models
        % select tables for each space
        [~,pm_tuningTable] = getNTidx(crossTuning,'spaceNum',1);
        [~,dl_tuningTable] = getNTidx(crossTuning,'spaceNum',2);

        % compose shift table for this model/bootstrap sample
        key_cols = contains(pm_tuningTable.Properties.VariableDescriptions,'meta');
        shift_tables{modelnum} = pm_tuningTable(:,key_cols);

        % remove spaceNum from columns
        shift_tables{modelnum}.spaceNum = [];

        % get PDs from pm and dl
        weights = pm_tuningTable.([model_names{modelnum} '_velWeight']);
        [pm_PDs,pm_moddepth] = cart2pol(weights(:,1),weights(:,2));
        weights = dl_tuningTable.([model_names{modelnum} '_velWeight']);
        [dl_PDs,dl_moddepth] = cart2pol(weights(:,1),weights(:,2));
        dPDs = minusPi2Pi(dl_PDs-pm_PDs);
        % use log for moddepth difference because of glm link?
        dMod = log(dl_moddepth)-log(pm_moddepth);

        tab_append = table(dPDs,dMod,'VariableNames',{'velPD','velModdepth'});
        tab_append.Properties.VariableDescriptions = {'circular','linear'};
        shift_tables{modelnum} = [shift_tables{modelnum} tab_append];

        mean_shifts{modelnum} = neuronAverage(shift_tables{modelnum},contains(shift_tables{modelnum}.Properties.VariableDescriptions,'meta'));
    end

    colors = {'r','g','b'};
    markers = {'x','+','.'};
    markersize = [15,15,40];
    titles = {'Hand-based model PD shift vs Actual PD shift','Egocentric model PD shift vs Actual PD shift','Muscle-based model PD shift vs Actual PD shift'};
    for modelnum = 1:3
        % [~,real_shifts] = getNTidx(shift_tables{4},'signalID',tunedNeurons{4});
        % [~,model_shifts] = getNTidx(shift_tables{modelnum},'signalID',tunedNeurons{4});
        % real_shifts = shift_tables{4};
        % model_shifts = shift_tables{modelnum};
        % comparePDClouds(real_shifts,model_shifts,struct('filter_tuning',[1]),colors{modelnum},'linewidth',1.85)
        % comparePDClouds(real_shifts,model_shifts,struct('filter_tuning',[]),colors{modelnum},'facealpha',0.1)
        [~,real_shifts] = getNTidx(mean_shifts{4},'signalID',tunedNeurons{4});
        [~,model_shifts] = getNTidx(mean_shifts{modelnum},'signalID',tunedNeurons{4});
        figure
        plot([-pi pi],[0 0],'-k','linewidth',2)
        hold on
        plot([0 0],[-pi pi],'-k','linewidth',2)
        plot([-pi pi],[-pi pi],'--k','linewidth',2)
        axis equal
        set(gca,'box','off','tickdir','out','xtick',[-pi pi],'ytick',[-pi pi],'xlim',[-pi pi],'ylim',[-pi pi],...
            'xticklabel',{'-\pi','\pi'},'yticklabel',{'-\pi','\pi'})
        scatter(real_shifts.velPD,model_shifts.velPD,50,colors{modelnum},'filled')
        xlabel 'Actual PD Shift'
        ylabel 'Modeled PD Shift'
        title(titles{modelnum})
    end

    % clean up
    clearvars colors titles

%% Calculate mean error on shifts
    err = zeros(100,3);
    for modelnum = 1:3
        [~,real_shifts] = getNTidx(shift_tables{4},'signalID',tunedNeurons{4});
        [~,model_shifts] = getNTidx(shift_tables{modelnum},'signalID',tunedNeurons{4});
        err_arr = model_shifts.velPD-real_shifts.velPD;
        for i = 1:100
            err_idx = 1:length(tunedNeurons{4});
            err_idx = err_idx + (i-1)*length(tunedNeurons{4});
            % use a 1-cos style error because of circular data
            % This value will range between 0 and 2
            err(i,modelnum) = mean(1-cos(err_arr(err_idx)));
        end
    end

    % plot histograms
    figure
    for modelnum = 1:3
        h = histogram(gca,err(:,modelnum),'DisplayStyle','bar');
        set(h,'edgecolor',colors{modelnum},'facecolor',colors{modelnum})
        hold on
    end

    % plot errors
    figure
    for modelnum = 1:3
        scatter(err(:,modelnum),repmat(modelnum/10,size(err,1),1),50,colors{modelnum},'filled')
        hold on
        plot(mean(err(:,modelnum)),modelnum/10,'k.','linewidth',3,'markersize',40)
    end
    set(gca,'tickdir','out','box','off','ytick',(1:3)/10,'yticklabel',{'Hand-based','Egocentric','Muscle-based'},'xtick',[0 1])
    axis equal
    axis ij
    xlim([0 1])
    ylim([0 4/10])
    xlabel('Cosine error of model')

    % compute statistics
    alpha = 0.05/3; % bonferroni correction for multiple comparisons...?
    diffstat = err(:,3)-err(:,1);
    mudiff = mean(diffstat);
    vardiff = var(diffstat);
    correction = 1/100 + 1/4;
    alphaup = 1-alpha/2;
    alphalow = alpha/2;
    upp = tinv(alphaup,99);
    low = tinv(alphalow,99);
    CIhigh = mudiff + upp * sqrt(correction*vardiff);
    CIlow = mudiff + low * sqrt(correction*vardiff);

%% Plot tuning weight clouds
    tuningHull = getTuningHull(results.tuningTable);
    % loop over each unit in one workspace
    % n_rows = ceil(sqrt(height(signalIDs)+1));
    cloud_fig = figure;
    surf_fig = figure;
    % neuron_idx = getNTidx(pdTables{1,4},'signalID',[95 2]);
    for neuron_idx = 1:height(tuningHull)
        close_fig = figure;

        figure(cloud_fig)
        clf
        plotMWTuningCloud(tuningHull,neuron_idx)

        figure(surf_fig)
        clf
        plotMWTuningSurfaces(results.td_test,[],neuron_idx)

        waitfor(close_fig)
    end

%% Plot DL vs PM just for neurons actually tuned to velocity
    % figure
    % comparePDs(pm_real_pdTable(isTuned_real,:),dl_real_pdTable(isTuned_real,:),struct('move_corr','vel'),'ko','linewidth',2)
    % hold on
    % comparePDs(pm_ext_pdTable(isTuned_real,:),dl_ext_pdTable(isTuned_real,:),struct('move_corr','vel'),'ro','linewidth',2)
    % comparePDs(pm_ego_pdTable(isTuned_real,:),dl_ego_pdTable(isTuned_real,:),struct('move_corr','vel'),'go','linewidth',2)
    % comparePDs(pm_musc_pdTable(isTuned_real,:),dl_musc_pdTable(isTuned_real,:),struct('move_corr','vel'),'bo','linewidth',2)
    % xlabel 'PM preferred direction'
    % ylabel 'DL preferred direction'

%% Plot tuning pR2 distribution for each neuron
    neuron_list = trial_data(1).S1_unit_guide;
    figure
    for neuron_idx = 1:length(neuron_list)
        clf
        for spacenum = 1:2
            [~,temp] = getNTidx(crossTuning,'signalID',neuron_list(neuron_idx,:),'spaceNum',spacenum);

            subplot(2,1,spacenum)
            hist(temp.S1_FR_eval)
            title(sprintf('Neuron %d %d, Space %d',neuron_list(neuron_idx,:),spacenum))
        end
        if ~isTuned{4}(neuron_idx)
            xlabel 'This is not tuned'
        end
        waitforbuttonpress
    end
