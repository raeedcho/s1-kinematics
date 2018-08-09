%% Split up trial data and preprocess
    % prep trial data by getting only rewards and trimming to only movements
    % first process marker data
    td = trial_data;
    td = smoothSignals(td,struct('signals','markers'));
    td = getDifferential(td,struct('signal','markers','alias','marker_vel'));
    [~,td] = getTDidx(td,'result','R');
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

    % Split td into different workspaces (workspace 1 is PM and workspace 2 is DL)
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    [~,td_pm] = getTDidx(td,'spaceNum',1);
    [~,td_dl] = getTDidx(td,'spaceNum',2);
    minsize = min(length(td_pm),length(td_dl));
    td_pm = td_pm(1:minsize);
    td_dl = td_dl(1:minsize);

    % recombine for later...
    td = [td_pm td_dl];

    %% Do PCA on muscle space
    % do PCA on muscles, training on only the training set
    % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    %                     'do_plot',true);
    PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len'))}}, 'do_plot',true);
    [td,~] = dimReduce(td,PCAparams);
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
    [td_temp,~] = dimReduce(td,PCAparams_vel);
    % temporary hack to allow us to save into something useful
    for i=1:length(td)
        td(i).opensim_muscVel_pca = td(i).opensim_pca;
    end
    % get rid of superfluous PCA
    td = rmfield(td,'opensim_pca');

    % Get PCA for neural space
    % PCAparams = struct('signals',{{'S1_spikes'}}, 'do_plot',true,'pca_recenter_for_proj',true,'sqrt_transform',true);
    % [td,~] = dimReduce(td,PCAparams);

    % Get PCA for marker space
    td = dimReduce(td,struct('signals','markers'));
    td = dimReduce(td,struct('signals','marker_vel'));

%% Set up model variables
    num_folds = 5; % 5 is default number of folds, no need to pass in
    num_repeats = 20; % 20 is default number of repeats, no need to pass in
    model_type = 'glm';
    model_aliases = {'ext','ego','musc','markers'};
    model_names = [strcat(model_type,'_',model_aliases,'_model') {'S1_FR'}];
    num_models = length(model_names);
    model_titles = cell(num_models-1,1);
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
    % colors for models
    model_colors = zeros(num_models-1,3);

    % colors for pm, dl conditions
    cond_colors = [0.6,0.5,0.7;...
                   1,0,0];
    
    % set up glm parameters
    num_musc_pcs = 5;
    neural_signals = 'S1_FR';
    % indices for cartesian hand coordinates
    opensim_hand_idx = find(contains(td(1).opensim_names,'_handPos') | contains(td(1).opensim_names,'_handVel'));
    opensim_joint_idx = find(contains(td(1).opensim_names,'_ang') | contains(td(1).opensim_names,'_vel'));
    glm_params = cell(num_models-1,1);
    for modelnum = 1:num_models-1
        switch model_aliases{modelnum}
        case 'musc'
            model_colors(modelnum,:) = [0, 174, 239]/255;
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','musc_model',...
                                    'in_signals',{{'opensim_len_pca',1:num_musc_pcs;'opensim_muscVel_pca',1:num_musc_pcs}},...
                                    'out_signals',neural_signals);
        case 'ext'
            model_colors(modelnum,:) = [247, 148, 30]/255;
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','ext_model',...
                                    'in_signals',{{'opensim',opensim_hand_idx}},...
                                    'out_signals',neural_signals);
        case 'ego'
            model_colors(modelnum,:) = [105, 189, 69]/255;
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','ego_model',...
                                    'in_signals',{{'sphere_hand_pos';'sphere_hand_vel'}},...
                                    'out_signals',neural_signals);
        case 'cyl'
            model_colors(modelnum,:) = [113, 191, 110]/255;
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','cyl_model',...
                                    'in_signals',{{'cyl_hand_pos';'cyl_hand_vel'}},...
                                    'out_signals',neural_signals);
        case 'joint'
            model_colors(modelnum,:) = [38, 34, 98]/255;
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','joint_model',...
                                    'in_signals',{{'opensim',opensim_joint_idx}},...
                                    'out_signals',neural_signals);
        case 'markers'
            model_colors(modelnum,:) = [193, 25, 47]/255;
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','markers_model',...
                                    'in_signals',{{'markers_pca',1:num_musc_pcs;'marker_vel_pca',1:num_musc_pcs}},...
                                    'out_signals',neural_signals);
        end
    end

%% Plot comparison of actual tuning curves with various modeled tuning curves
    % Split td into different workspaces (workspace 1 is PM and workspace 2 is DL)
    [~,td_pm] = getTDidx(td,'spaceNum',1);
    [~,td_dl] = getTDidx(td,'spaceNum',2);

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
        figure
        title(model_names{modelnum})
        compareTuning(tuning_curves(:,modelnum),pdTables(:,modelnum),struct('which_units',find(isTuned),'cond_colors',cond_colors))
        % compareTuning(tuning_curves(:,modelnum),pdTables(:,modelnum),struct('which_units',find(isTuned),'cond_colors',cond_colors,'maxFR',1))
    end

%% Make iris and dna plots
    % get names of tuned neurons
    signalIDs = td_pm(1).S1_unit_guide;
    tunedNeurons = signalIDs(isTuned,:);

    f1 = figure;
    f2 = figure;
    for modelnum = 1:num_models
        figure(f1)
        subplot(2,ceil(num_models/2),modelnum)
        irisPlot(pdTables{1,modelnum}(isTuned,:),pdTables{2,modelnum}(isTuned,:));
        title(model_names{modelnum})

        figure(f2)
        subplot(2,ceil(num_models/2),modelnum)
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

    figure
    for modelnum = 1:num_models-1
        [~,real_shifts] = getNTidx(mean_shifts{end},'signalID',tunedNeurons);
        [~,model_shifts] = getNTidx(mean_shifts{modelnum},'signalID',tunedNeurons);

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
    

