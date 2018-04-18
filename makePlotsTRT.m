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

%% Split up trial data and preprocess
    % prep trial data by getting only rewards and trimming to only movements
    [~,td] = getTDidx(trial_data,'result','R');
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    % for bumps
    % [~,td] = getTDidx(trial_data,'result','R');
    % td = td(~isnan(cat(1,td.idx_bumpTime)));
    % td = trimTD(td,{'idx_bumpTime',0},{'idx_bumpTime',15});

    % bin data at 50ms
    td = binTD(td,5);
    % add in spherical coordinates
    td = addSphereHand2TD(td);
    % add in cylindrical coordinates
    td = addCylHand2TD(td);
    % add firing rates rather than spike counts
    td = addFiringRates(td,struct('array','S1'));

    %% Do PCA on muscle space
    % do PCA on muscles, training on only the training set
    % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    %                     'do_plot',true);
    PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len'))}}, 'do_plot',true);
    [td,~] = getPCA(td,PCAparams);
    % temporary hack to allow us to do PCA on velocity too
    for i=1:length(td)
        td(i).opensim_len_pca = td(i).opensim_pca;
    end
    % get rid of superfluous PCA
    td = rmfield(td,'opensim_pca');
    % get velocity PCA
    % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    %                     'do_plot',true);
    PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel'))}}, 'do_plot',true);
    [td,~] = getPCA(td,PCAparams_vel);
    % temporary hack to allow us to save into something useful
    for i=1:length(td)
        td(i).opensim_muscVel_pca = td(i).opensim_pca;
    end
    % get rid of superfluous PCA
    td = rmfield(td,'opensim_pca');

    % Get PCA for neural space
    % PCAparams = struct('signals',{{'S1_spikes'}}, 'do_plot',true,'pca_recenter_for_proj',true,'sqrt_transform',true);
    % [td,~] = getPCA(td,PCAparams);

    % Split td into different workspaces (workspace 1 is PM and workspace 2 is DL)
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    [~,td_pm] = getTDidx(td,'spaceNum',1);
    [~,td_dl] = getTDidx(td,'spaceNum',2);
    minsize = min(length(td_pm),length(td_dl));
    td_pm = td_pm(1:minsize);
    td_dl = td_dl(1:minsize);

    % recombine for later...
    td = [td_pm td_dl];

%% Set up model variables
    num_folds = 5; % 5 is default number of folds, no need to pass in
    num_repeats = 20; % 20 is default number of repeats, no need to pass in
    model_type = 'glm';
    model_names = [strcat(model_type,'_',{'musc','ext','ego','cyl','joint'},'_model') {'S1_FR'}];
    num_models = length(model_names);
    % colors for models
    model_colors = [  0,174,239;...
                    247,148, 30;...
                    105,189, 69;...
                    113,191,110;...
                     38, 34, 98]/255;

    % colors for pm, dl conditions
    cond_colors = [0.6,0.5,0.7;...
                   1,0,0];
    
    % set up glm parameters
    num_musc_pcs = 5;
    neural_signals = 'S1_FR';
    % indices for cartesian hand coordinates
    opensim_hand_idx = find(contains(td(1).opensim_names,'_handPos') | contains(td(1).opensim_names,'_handVel'));
    opensim_joint_idx = find(contains(td(1).opensim_names,'_ang') | contains(td(1).opensim_names,'_vel'));
    glm_params{1} = struct('model_type',model_type,...
                            'model_name','musc_model',...
                            'in_signals',{{'opensim_len_pca',1:num_musc_pcs;'opensim_muscVel_pca',1:num_musc_pcs}},...
                            'out_signals',neural_signals);
    glm_params{2} = struct('model_type',model_type,...
                            'model_name','ext_model',...
                            'in_signals',{{'opensim',opensim_hand_idx}},...
                            'out_signals',neural_signals);
    glm_params{3} = struct('model_type',model_type,...
                            'model_name','ego_model',...
                            'in_signals',{{'sphere_hand_pos';'sphere_hand_vel'}},...
                            'out_signals',neural_signals);
    glm_params{4} = struct('model_type',model_type,...
                            'model_name','cyl_model',...
                            'in_signals',{{'cyl_hand_pos';'cyl_hand_vel'}},...
                            'out_signals',neural_signals);
    glm_params{5} = struct('model_type',model_type,...
                            'model_name','joint_model',...
                            'in_signals',{{'opensim',opensim_joint_idx}},...
                            'out_signals',neural_signals);

%% Plot comparison of actual tuning curves with various modeled tuning curves
    % use K-fold crossvalidation to get neural predictions from each model for tuning curves and PDs
    indices = crossvalind('Kfold',length(td_pm),num_folds);
    td_test = cell(2,num_folds);

    % do the crossval
    for foldctr = 1:num_folds
        % split into testing and training
        test_idx = (indices==foldctr);
        train_idx = ~test_idx;
        td_train = [td_pm(train_idx) td_dl(train_idx)];
        td_test{1,foldctr} = td_pm(test_idx);
        td_test{2,foldctr} = td_dl(test_idx);

        % Fit models on training data
        for modelnum = 1:num_models-1
            [~,glm_info] = getModel(td_train,glm_params{modelnum});

            % predict firing rates for td_test
            for spacenum = 1:2
                td_test{spacenum,foldctr} = getModel(td_test{spacenum,foldctr},glm_info);
            end
        end
    end
    td_tuning = {horzcat(td_test{1,:}); horzcat(td_test{2,:})};
    
    % get PDs and tuning curves
    pdTables = cell(2,num_models);
    tuning_curves = cell(2,num_models);
    for modelnum = 1:num_models
        for spacenum = 1:2
            % First PDs
            pd_params = struct('out_signals',model_names{modelnum},'out_signal_names',td(1).S1_unit_guide,'do_plot',false,'meta',struct('spaceNum',spacenum));
            % pdTables{spacenum,modelnum} = getTDClassicalPDs(td_tuning{spacenum},pd_params);
            pdTables{spacenum,modelnum} = getTDPDs(td_tuning{spacenum},pd_params);

            tuning_params = struct('out_signals',model_names{modelnum},'out_signal_names',td(1).S1_unit_guide,'meta',struct('spaceNum',spacenum));
            tuning_curves{spacenum,modelnum} = getTuningCurves(td_tuning{spacenum},tuning_params);
        end
    end
    isTuned = pdTables{1,end}.velTuned & pdTables{2,end}.velTuned;

    % compare PM and DL tuning for each model
    for modelnum = 1:num_models
        figure;compareTuning(tuning_curves(:,modelnum),pdTables(:,modelnum),struct('which_units',find(isTuned),'cond_colors',cond_colors))
    end

%% Make iris and dna plots
    % get names of tuned neurons
    signalIDs = td_pm(1).S1_unit_guide;
    tunedNeurons = signalIDs(isTuned,:);

    f1 = figure;
    f2 = figure;
    for modelnum = 1:num_models
        figure(f1)
        subplot(2,num_models/2,modelnum)
        irisPlot(pdTables{1,modelnum}(isTuned,:),pdTables{2,modelnum}(isTuned,:));
        title(model_names{modelnum})

        figure(f2)
        subplot(2,num_models/2,modelnum)
        dnaPlot(pdTables{1,modelnum}(isTuned,:),pdTables{2,modelnum}(isTuned,:));
        title(model_names{modelnum})
    end
    
%% Cross-validate models of neural data
    % Get crossval info
    crossval_params = struct('model_names',{model_names},'glm_params',{glm_params});
    [crossEval,crossTuning] = analyzeTRT(td,crossval_params);

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
        pm_PDs = pm_tuningTable.([model_names{modelnum} '_velPD']);
        dl_PDs = dl_tuningTable.([model_names{modelnum} '_velPD']);
        dPDs = minusPi2Pi(dl_PDs-pm_PDs);

        tab_append = table(dPDs,'VariableNames',{'velPD'});
        tab_append.Properties.VariableDescriptions = {'circular'};
        shift_tables{modelnum} = [shift_tables{modelnum} tab_append];

        mean_shifts{modelnum} = neuronAverage(shift_tables{modelnum},contains(shift_tables{modelnum}.Properties.VariableDescriptions,'meta'));
    end

    markers = {'x','+','.'};
    markersize = [15,15,40];
    titles = {'Muscle-based model PD shift vs Actual PD shift',...
        'Hand-based model PD shift vs Actual PD shift',...
        'Egocentric model PD shift vs Actual PD shift',...
        'Cylindrical ego model PD shift vs Actual PD shift',...
        'Joint-based model PD shift vs Actual PD shift'};
    for modelnum = 1:num_models-1
        [~,real_shifts] = getNTidx(mean_shifts{end},'signalID',tunedNeurons);
        [~,model_shifts] = getNTidx(mean_shifts{modelnum},'signalID',tunedNeurons);
        figure
        plot([-pi pi],[0 0],'-k','linewidth',2)
        hold on
        plot([0 0],[-pi pi],'-k','linewidth',2)
        plot([-pi pi],[-pi pi],'--k','linewidth',2)
        axis equal
        set(gca,'box','off','tickdir','out','xtick',[-pi pi],'ytick',[-pi pi],'xlim',[-pi pi],'ylim',[-pi pi],...
            'xticklabel',{'-\pi','\pi'},'yticklabel',{'-\pi','\pi'})
        scatter(real_shifts.velPD,model_shifts.velPD,50,model_colors(modelnum,:),'filled')
        xlabel 'Actual PD Shift'
        ylabel 'Modeled PD Shift'
        title(titles{modelnum})
    end

%% Calculate mean error on shifts
    err = zeros(100,num_models-1);
    for modelnum = 1:num_models-1
        [~,real_shifts] = getNTidx(shift_tables{end},'signalID',tunedNeurons);
        [~,model_shifts] = getNTidx(shift_tables{modelnum},'signalID',tunedNeurons);
        err_arr = model_shifts.velPD-real_shifts.velPD;
        for i = 1:100
            err_idx = 1:length(tunedNeurons);
            err_idx = err_idx + (i-1)*length(tunedNeurons);
            % use a 1-cos style error because of circular data
            % This value will range between 0 and 2
            err(i,modelnum) = mean(1-cos(err_arr(err_idx)));
        end
    end

    % plot histograms
    figure
    for modelnum = 1:num_models-1
        h = histogram(gca,err(:,modelnum),'DisplayStyle','bar');
        set(h,'edgecolor',model_colors(modelnum,:),'facecolor',model_colors(modelnum,:))
        hold on
    end

    % plot errors
    figure
    for modelnum = 1:num_models-1
        scatter(err(:,modelnum),repmat(modelnum/10,size(err,1),1),50,model_colors(modelnum,:),'filled')
        hold on
        plot(mean(err(:,modelnum)),modelnum/10,'k.','linewidth',3,'markersize',40)
    end
    set(gca,'tickdir','out','box','off','ytick',(1:(num_models-1))/10,'yticklabel',{'Muscle-based','Hand-based','Egocentric','Cylindrical','Joint-based'},'xtick',[0 1])
    axis equal
    axis ij
    xlim([0 1])
    ylim([0 num_models/10])
    xlabel('Cosine error of model')

    % compute statistics
    alpha = 0.05; % bonferroni correction for multiple comparisons...?
    diffstat = err(:,1)-err(:,2); % musc - ext
    mudiff = mean(diffstat);
    vardiff = var(diffstat);
    correction = 1/100 + 1/4;
    alphaup = 1-alpha/2;
    alphalow = alpha/2;
    upp = tinv(alphaup,99);
    low = tinv(alphalow,99);
    errCIhigh = mudiff + upp * sqrt(correction*vardiff);
    errCIlow = mudiff + low * sqrt(correction*vardiff);

%% Plot pR2s against each other
    % av_pR2_ext_dl = mean(td_ext_dl_eval,2);
    % av_pR2_musc_dl = mean(td_musc_dl_eval,2);
    % av_pR2_ext_pm = mean(td_ext_pm_eval,2);
    % av_pR2_musc_pm = mean(td_musc_pm_eval,2);
    avgEval = neuronAverage(crossEval,contains(crossEval.Properties.VariableDescriptions,'meta'));
    av_pR2_ext = avgEval.glm_ext_model_eval;
    av_pR2_musc = avgEval.glm_musc_model_eval;

    % make plot of pR2 of muscle against ext
    figure
    plot([-1 1],[-1 1],'k--','linewidth',2)
    hold on
    plot([0 0],[-1 1],'k-','linewidth',2)
    plot([-1 1],[0 0],'k-','linewidth',2)
    scatter(crossEval.glm_ext_model_eval(:),crossEval.glm_musc_model_eval(:),25,'k','filled')
    scatter(av_pR2_ext(:),av_pR2_musc(:),50,'r','filled')
    set(gca,'box','off','tickdir','out','xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
    xlabel 'Hand-based pR2'
    ylabel 'Muscle-based pR2'

    % plot histogram of dpR2
    diffstat = crossEval.glm_musc_model_eval-crossEval.glm_ext_model_eval;
    figure
    histogram(diffstat)
    hold on
    plot([0 0],get(gca,'ylim'),'--k','linewidth',3)

    % get stats on pR2 diff between musc model and ext model
    alpha = 0.05; % bonferroni correction for multiple comparisons...?
    mudiff = mean(diffstat);
    vardiff = var(diffstat);
    correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
    alphaup = 1-alpha/2;
    alphalow = alpha/2;
    upp = tinv(alphaup,num_folds*num_repeats-1);
    low = tinv(alphalow,num_folds*num_repeats-1);
    dpR2CIhigh = mudiff + upp * sqrt(correction*vardiff);
    dpR2CIlow = mudiff + low * sqrt(correction*vardiff);

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

%% Plot tuning pR2 distribution for each neuron
    % neuron_list = trial_data(1).S1_unit_guide;
    % figure
    % for neuron_idx = 1:length(neuron_list)
    %     clf
    %     for spacenum = 1:2
    %         [~,temp] = getNTidx(crossTuning,'signalID',neuron_list(neuron_idx,:),'spaceNum',spacenum);

    %         subplot(2,1,spacenum)
    %         hist(temp.S1_FR_eval)
    %         title(sprintf('Neuron %d %d, Space %d',neuron_list(neuron_idx,:),spacenum))
    %     end
    %     if ~isTuned(neuron_idx)
    %         xlabel 'This is not tuned'
    %     end
    %     waitforbuttonpress
    % end

