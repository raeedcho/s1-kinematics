function results = mwEncoders(td,params)
%% Split up trial data and preprocess
    % Split td into different workspaces (workspace 1 is PM and workspace 2 is DL)
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    [~,td_pm] = getTDidx(td,'spaceNum',1);
    [~,td_dl] = getTDidx(td,'spaceNum',2);
    minsize = min(length(td_pm),length(td_dl));
    td_pm = td_pm(1:minsize);
    td_dl = td_dl(1:minsize);

    % recombine for later...
    td = [td_pm td_dl];

    % Get PCA for neural space
    % PCAparams = struct('signals',{{'S1_spikes'}}, 'do_plot',true,'pca_recenter_for_proj',true,'sqrt_transform',true);
    % [td,~] = dimReduce(td,PCAparams);

%% Set up model variables
    num_folds = 5; % 5 is default number of folds, no need to pass in
    num_repeats = 20; % 20 is default number of repeats, no need to pass in
    num_tuning_bins = 16;
    model_type = 'glm';
    model_aliases = {'ext','ego','musc','markers'};
    neural_signals = 'S1_FR';
    model_names = [strcat(model_type,'_',model_aliases,'_model') {neural_signals}];
    num_models = length(model_names);
    
    % set up glm parameters
    num_musc_pcs = 5;
    glm_params = cell(num_models-1,1);
    for modelnum = 1:num_models-1
        switch model_aliases{modelnum}
        case 'musc'
            %% Do PCA on muscle space
            % do PCA on muscles, training on only the training set
            % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
            % PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
            %                     'do_plot',true);
            PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len'))}}, 'do_plot',false);
            [td_temp,~] = dimReduce(td,PCAparams);
            % temporary hack to allow us to do PCA on velocity too
            for i=1:length(td)
                td(i).opensim_len_pca = td_temp(i).opensim_pca;
            end
            % get velocity PCA
            % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
            % PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
            %                     'do_plot',true);
            PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel'))}}, 'do_plot',false);
            [td_temp,~] = dimReduce(td,PCAparams_vel);
            % temporary hack to allow us to save into something useful
            for i=1:length(td)
                td(i).opensim_muscVel_pca = td_temp(i).opensim_pca;
            end
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','musc_model',...
                                    'in_signals',{{'opensim_len_pca',1:num_musc_pcs;'opensim_muscVel_pca',1:num_musc_pcs}},...
                                    'out_signals',neural_signals);
        case 'ext'
            marker_hand_idx = 1:3;
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','ext_model',...
                                    'in_signals',{{'markers',marker_hand_idx;'marker_vel',marker_hand_idx}},...
                                    'out_signals',neural_signals);
            % opensim_hand_idx = find(contains(td(1).opensim_names,'_handPos') | contains(td(1).opensim_names,'_handVel'));
            % glm_params{modelnum} = struct('model_type',model_type,...
            %                         'model_name','ext_model',...
            %                         'in_signals',{{'opensim',opensim_hand_idx}},...
            %                         'out_signals',neural_signals);
        case 'ego'
            % add in spherical coordinates
            td = addSphereHand2TD(td);
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','ego_model',...
                                    'in_signals',{{'sphere_hand_pos';'sphere_hand_vel'}},...
                                    'out_signals',neural_signals);
        case 'cyl'
            % add in cylindrical coordinates
            td = addCylHand2TD(td);
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','cyl_model',...
                                    'in_signals',{{'cyl_hand_pos';'cyl_hand_vel'}},...
                                    'out_signals',neural_signals);
        case 'joint'
            opensim_joint_idx = find(contains(td(1).opensim_names,'_ang') | contains(td(1).opensim_names,'_vel'));
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','joint_model',...
                                    'in_signals',{{'opensim',opensim_joint_idx}},...
                                    'out_signals',neural_signals);
        case 'markers'
            % indices for cartesian hand coordinates
            marker_hand_idx = 1:3;
            marker_elbow_idx = 28:30;
            glm_params{modelnum} = struct('model_type',model_type,...
                                    'model_name','markers_model',...
                                    'in_signals',{{'markers',[marker_hand_idx marker_elbow_idx];'marker_vel',[marker_hand_idx marker_elbow_idx]}},...
                                    'out_signals',neural_signals);

            % Get PCA for marker space
            % td = dimReduce(td,struct('signals','markers'));
            % td = dimReduce(td,struct('signals','marker_vel'));
            % glm_params{modelnum} = struct('model_type',model_type,...
            %                         'model_name','markers_model',...
            %                         'in_signals',{{'markers_pca',1:num_musc_pcs;'marker_vel_pca',1:num_musc_pcs}},...
            %                         'out_signals',neural_signals);
        end
    end
    assignParams(who,params);

%% Get comparison of actual tuning curves with various modeled tuning curves
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
            pdTables{spacenum,modelnum} = getTDClassicalPDs(td_tuning{spacenum},pd_params);
            % pdTables{spacenum,modelnum} = getTDPDs(td_tuning{spacenum},pd_params);

            tuning_params = struct('out_signals',model_names{modelnum},'out_signal_names',td(1).S1_unit_guide,...
                'num_bins',num_tuning_bins,'meta',struct('spaceNum',spacenum));
            tuning_curves{spacenum,modelnum} = getTuningCurves(td_tuning{spacenum},tuning_params);
        end
    end

%% Cross-validate models of neural data
    % Get crossval info
    crossval_params = struct('model_names',{model_names},'glm_params',{glm_params},'num_folds',num_folds,'num_repeats',num_repeats);
    [crossEval,crossTuning] = analyzeTRT(td,crossval_params);
    
%% Make iris and dna plots
    % f1 = figure;
    % f2 = figure;
    % for modelnum = 1:num_models
    %     figure(f1)
    %     subplot(2,ceil(num_models/2),modelnum)
    %     irisPlot(pdTables{1,modelnum}(isTuned,:),pdTables{2,modelnum}(isTuned,:));
    %     title(model_names{modelnum})

    %     figure(f2)
    %     subplot(2,ceil(num_models/2),modelnum)
    %     dnaPlot(pdTables{1,modelnum}(isTuned,:),pdTables{2,modelnum}(isTuned,:));
    %     title(model_names{modelnum})
    % end

%% create return struct
    % for cross validation plots
    results.tuning_curves = tuning_curves;
    results.pdTables = pdTables;
    results.crossEval = crossEval;
    results.crossTuning = crossTuning;

    % for showing predictive capability
    results.td_tuning = td_tuning;

    % get names of tuned neurons
    signalIDs = td(1).S1_unit_guide;
    results.isTuned = pdTables{1,end}.velTuned & pdTables{2,end}.velTuned;
    results.tunedNeurons = signalIDs(results.isTuned,:);

    % get parameters
    results.params.num_folds = num_folds;
    results.params.num_repeats = num_repeats;
    results.params.num_tuning_bins = num_tuning_bins;
    results.params.model_type = model_type;
    results.params.model_aliases = model_aliases;
    results.params.model_names = model_names;
    results.params.num_musc_pcs = num_musc_pcs;
    results.params.neural_signals = neural_signals;
    results.params.glm_params = glm_params;