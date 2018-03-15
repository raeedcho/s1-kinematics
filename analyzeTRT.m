%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function results = analyzeTRT(trial_data,params)
% 
% For multiworkspace files, with dl and pm workspaces:
%   * Fits three different coordinate frame models to data from both workspaces
%   * Predicts firing rates from models on test data in both workspaces
%   * Compare predicted firing rates to actual firing rates in test data
%
% INPUTS:
%   trial_data : the struct
%   params     : parameter struct
%       .neural_signals  : which signals to calculate PDs for
%                           default: 'S1_spikes'
%       .model_type     :   type of model to fit
%                           default; 'glm'
%       .glm_distribution : distribution to use for GLM
%                           default: 'poisson'
%       .model_eval_metric : Evaluation metric for models
%                           default: 'pr2'
%       .num_musc_pcs   :   Number of muscle principle components to use
%                           default: 5
%       .num_boots    : # bootstrap iterations to use
%                           default: 100
%       .verbose      :     Print diagnostic plots and information
%
% OUTPUTS:
%   result  : results of analysis
%       .td_eval        : evaluation metrics for various models, in cell array
%                           rows: PM, DL, full
%                           cols: ext, ego, musc
%       .pdTables       : PD tables in cell array
%                           Rows: PM, DL
%                           Cols: ext_model, ego_model, musc_model, real
%       .tuning_curves  : tuning curves in cell array, same as pdTables
%       .shift_tables   : PD shift tables in cell array, cols same as above
%       .glm_info       : output from fitting GLMs to different models
%                           Cols: ext, ego, musc
%       .isTuned        : cell array of isTuned, one for each model, as above
%       .td_train       : trial data structure used to train models
%       .td_test        : trial data structures used to test models
%                           PM is first, DL is second
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function results = analyzeTRT(trial_data,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up
    % default parameters
    neural_signals = 'S1_spikes';
    glm_distribution = 'poisson';
    model_eval_metric = 'pr2';
    model_type = 'glm';
    num_musc_pcs = 5;
    num_boots = 100;
    verbose = false;
    if nargin > 1, assignParams(who,params); end % overwrite parameters

%% Compile training and test sets
    % prep trial data by getting only rewards and trimming to only movements
    [~,td] = getTDidx(trial_data,'result','R');
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    % add in spherical coordinates
    td = addSphereHand2TD(td);

    % add firing rates rather than spike counts
    td = addFiringRates(td,struct('array','S1'));

    % bin data at 50ms
    td = binTD(td,5);
    
    % Split td into different workspaces (workspace 1 is PM and workspace 2 is DL)
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    [~,td_pm] = getTDidx(td,'spaceNum',1);
    [~,td_dl] = getTDidx(td,'spaceNum',2);
    minsize = min(length(td_pm),length(td_dl));
    td_pm = td_pm(1:minsize);
    td_dl = td_dl(1:minsize);
    
    % get training and test set indices
    [train_idx,test_idx] = crossvalind('HoldOut',minsize,0.2);
    td_train = [td_pm(train_idx) td_dl(train_idx)];
    td_test = cell(2,1);
    td_test{1} = td_pm(test_idx);
    td_test{2} = td_dl(test_idx);

%% Get PCA for muscle space
    % do PCA on muscles, training on only the training set
    % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    %                     'do_plot',true);
    PCAparams = struct('signals',{{'opensim',find(contains(td_train(1).opensim_names,'_len'))}}, 'do_plot',verbose);
    [td_train,pca_info] = getPCA(td_train,PCAparams);
    for spacenum = 1:2
        td_test{spacenum} = getPCA(td_test{spacenum},pca_info);
    end
    
    % temporary hack to allow us to do PCA on velocity too
    for i=1:length(td_train)
        td_train(i).opensim_len_pca = td_train(i).opensim_pca;
    end
    for i=1:length(td_test{1})
        for spacenum = 1:2
            td_test{spacenum}(i).opensim_len_pca = td_test{spacenum}(i).opensim_pca;
        end
    end
    
    % get rid of superfluous PCA
    td_train = rmfield(td_train,'opensim_pca');
    for spacenum = 1:2
        td_test{spacenum} = rmfield(td_test{spacenum},'opensim_pca');
    end

    % get velocity PCA
    % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    %                     'do_plot',true);
    PCAparams_vel = struct('signals',{{'opensim',find(contains(td_train(1).opensim_names,'_muscVel'))}}, 'do_plot',verbose);
    [td_train,pca_info_vel] = getPCA(td_train,PCAparams_vel);
    for spacenum = 1:2
        td_test{spacenum} = getPCA(td_test{spacenum},pca_info_vel);
    end
    
    % temporary hack to allow us to save into something useful
    for i=1:length(td_train)
        td_train(i).opensim_muscVel_pca = td_train(i).opensim_pca;
    end
    for i=1:length(td_test{1})
        td_test{1}(i).opensim_muscVel_pca = td_test{1}(i).opensim_pca;
        td_test{2}(i).opensim_muscVel_pca = td_test{2}(i).opensim_pca;
    end
    
    % get rid of superfluous PCA
    td_train = rmfield(td_train,'opensim_pca');
    for spacenum = 1:2
        td_test{spacenum} = rmfield(td_test{spacenum},'opensim_pca');
    end

%% Get PCA for neural space
    PCAparams = struct('signals',{{'S1_spikes'}}, 'do_plot',verbose,'pca_recenter_for_proj',true,'sqrt_transform',true);
    [td_train,pca_info] = getPCA(td_train,PCAparams);
    for spacenum = 1:2
        td_test{spacenum} = getPCA(td_test{spacenum},pca_info);
    end

%% Fit models
    % set up parameters for models
    glm_params = cell(1,3);
    glm_info = cell(1,3);
    % cartesian hand coordinates for position and velocity
    model_names = {[model_type '_ext_model'],[model_type '_ego_model'],[model_type '_musc_model'],neural_signals};
    opensim_hand_idx = find(contains(td_train(1).opensim_names,'_handPos') | contains(td_train(1).opensim_names,'_handVel'));
    glm_params{1} = struct('model_type',model_type,...
                            'model_name','ext_model',...
                            'in_signals',{{'opensim',opensim_hand_idx}},...
                            'out_signals',{{neural_signals}},'glm_distribution',glm_distribution);
    glm_params{2} = struct('model_type',model_type,...
                            'model_name','ego_model',...
                            'in_signals',{{'sphere_hand_pos';'sphere_hand_vel'}},...
                            'out_signals',{{neural_signals}},'glm_distribution',glm_distribution);
    glm_params{3} = struct('model_type',model_type,...
                            'model_name','musc_model',...
                            'in_signals',{{'opensim_len_pca',1:num_musc_pcs;'opensim_muscVel_pca',1:num_musc_pcs}},...
                            'out_signals',{{neural_signals}},'glm_distribution',glm_distribution);
    for modelnum = 1:3
        [~,glm_info{modelnum}] = getModel(td_train,glm_params{modelnum});
    end

    % Predict firing rates
    for modelnum = 1:3
        for spacenum = 1:2
            td_test{spacenum} = getModel(td_test{spacenum},glm_info{modelnum});
        end
    end

    % Evaluate model fits
    % set up eval array
    % row 1 is PM, row 2 is DL, and row 3 is PM and DL together
    % col 1 is ext, col 2 is ego, col 3 is musc
    td_eval = cell(3,3);
    eval_params = glm_info;
    for modelnum = 1:3
        eval_params{modelnum}.eval_metric = model_eval_metric;
        for spacenum = 1:2
            td_eval{spacenum,modelnum} = squeeze(evalModel(td_test{spacenum},eval_params{modelnum}));
        end
        td_eval{end,modelnum} = squeeze(evalModel([td_test{2} td_test{1}],eval_params{modelnum}));
    end

%% Get PDs and tuning curves for the modeled and actual neurons
    % set up outputs
    pdTables = cell(2,4); % PM is first row, DL is second. Column order is Ext, Ego, Musc, Real
    tuning_curves = cell(2,4); % PM is first row, DL is second. Column order is Ext, Ego, Musc, Real
    isTuned = cell(1,4);
    pd_params = cell(1,4);
    tuning_params = cell(1,4);
    
    % set up parameters
    % PDs
    
    num_bins = 8;
    % get PDs and tuning curves
    for modelnum = 1:4
        pd_params{modelnum} = struct('num_boots',num_boots,'out_signals',{model_names(modelnum)},'out_signal_names',td_train(1).S1_unit_guide,'disp_times',verbose,'distribution',glm_distribution);
        tuning_params{modelnum} = struct('num_bins',num_bins,'out_signals',{model_names(modelnum)},'out_signal_names',td_train(1).S1_unit_guide);
    
        for spacenum = 1:2
            pdTables{spacenum,modelnum} = getTDPDs(td_test{spacenum},pd_params{modelnum});
            tuning_curves{spacenum,modelnum} = getTuningCurves(td_test{spacenum},tuning_params{modelnum});
        end
        isTuned{modelnum} = checkIsTuned(pdTables{1,modelnum}) & checkIsTuned(pdTables{2,modelnum});
    end

%% Bootstrap on PD shifts
    shift_tables = cell(1,length(model_names));
    num_internal_boots = 1;
    trial_idx = randi(length(td_test{1}),length(td_test{1}),num_boots);
    if verbose
        tic
    end
    for bootctr = 1:num_boots
        for modelctr = 1:length(model_names)
            pd_params = struct('num_boots',num_internal_boots,'out_signals',{model_names(modelctr)},'out_signal_names',td_train(1).S1_unit_guide,'distribution',glm_distribution);

            pm_pdTable = getTDPDs(td_test{1}(trial_idx(:,bootctr)),pd_params);
            dl_pdTable = getTDPDs(td_test{2}(trial_idx(:,bootctr)),pd_params);

            % compose shift table for this model/bootstrap sample
            temp_shift_table = pm_pdTable;
            temp_shift_table.velPD = minusPi2Pi(dl_pdTable.velPD-pm_pdTable.velPD);
            temp_shift_table.velPDCI = minusPi2Pi(dl_pdTable.velPDCI-pm_pdTable.velPDCI); % this doesn't actually mean anything with one internal boot sample
            % don't care about moddepth currently

            if bootctr == 1
                % slot new table into cell array
                shift_tables{modelctr} = temp_shift_table;
            else
                % append to old table
                shift_tables{modelctr} = [shift_tables{modelctr};temp_shift_table];
            end

        end
        if verbose
            disp(['Bootstrap sample ' num2str(bootctr) ', ending at ' num2str(toc) 's'])
        end
    end

%% Package up outputs
    results = struct('td_eval',{td_eval},...
                        'pdTables',{pdTables},...
                        'tuning_curves',{tuning_curves},...
                        'shift_tables',{shift_tables},...
                        'glm_info',{glm_info},...
                        'isTuned',{isTuned},...
                        'td_train',{td_train},...
                        'td_test',{td_test});
