% function fitPredictTRT(td)
% For multiworkspace files, with dl and pm workspaces:
%   * Fits three different coordinate frame models to data from both workspaces
%   * Predicts firing rates from models on test data in both workspaces
%   * Compare predicted firing rates to actual firing rates in test data

%% Compile training and test sets
    % prep trial data by getting only rewards and trimming to only movements
    [~,td] = getTDidx(trial_data,'result','R');
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    % add in spherical coordinates
    td = addSphereHand2TD(td);

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
    
    % clean up
    clearvars -except trial_data td_train td_test

%% Train models on training set composed of both workspaces together
    % do PCA on muscles, training on only the training set
    % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    %                     'do_plot',true);
    PCAparams = struct('signals',{{'opensim',find(contains(td_train(1).opensim_names,'_len'))}}, 'do_plot',true);
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
    PCAparams_vel = struct('signals',{{'opensim',find(contains(td_train(1).opensim_names,'_muscVel'))}}, 'do_plot',true);
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
    
    % clean up
    clearvars -except trial_data td_train td_test

%% Fit models
    % set up parameters for models
    glm_params = cell(1,3);
    glm_info = cell(1,3);
    % cartesian hand coordinates for position and velocity
    opensim_hand_idx = find(contains(td_train(1).opensim_names,'_handPos') | contains(td_train(1).opensim_names,'_handVel'));
    glm_params{1} = struct('model_type','glm',...
                            'model_name','ext_model',...
                            'in_signals',{{'opensim',opensim_hand_idx}},...
                            'out_signals',{{'S1_spikes'}});
    glm_params{2} = struct('model_type','glm',...
                            'model_name','ego_model',...
                            'in_signals',{{'sphere_hand_pos';'sphere_hand_vel'}},...
                            'out_signals',{{'S1_spikes'}});
    glm_params{3} = struct('model_type','glm',...
                            'model_name','musc_model',...
                            'in_signals',{{'opensim_len_pca',1:5;'opensim_muscVel_pca',1:5}},...
                            'out_signals',{'S1_spikes'});
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
        eval_params{modelnum}.eval_metric = 'pr2';
        for spacenum = 1:2
            td_eval{spacenum,modelnum} = squeeze(evalModel(td_test{spacenum},eval_params{modelnum}));
        end
        td_eval{end,modelnum} = squeeze(evalModel([td_test{2} td_test{1}],eval_params{modelnum}));
    end
    
    % clean up
    clearvars opensim_hand_idx modelnum spacenum 

%% Get PDs and tuning curves for the modeled and actual neurons
    % set up outputs
    pdTables = cell(2,4); % PM is first row, DL is second. Column order is Ext, Ego, Musc, Real
    tuning_curves = cell(2,4); % PM is first row, DL is second. Column order is Ext, Ego, Musc, Real
    isTuned = cell(1,4);
    pd_params = cell(1,4);
    tuning_params = cell(1,4);
    
    % set up parameters
    % PDs
    num_boots = 10;
    pd_params{1} = struct('num_boots',num_boots,'out_signals',{{'glm_ext_model'}},'out_signal_names',td_train(1).S1_unit_guide,'disp_times',true);
    pd_params{2} = struct('num_boots',num_boots,'out_signals',{{'glm_ego_model'}},'out_signal_names',td_train(1).S1_unit_guide,'disp_times',true);
    pd_params{3} = struct('num_boots',num_boots,'out_signals',{{'glm_musc_model'}},'out_signal_names',td_train(1).S1_unit_guide,'disp_times',true);
    pd_params{4} = struct('num_boots',num_boots,'out_signals',{{'S1_spikes'}},'out_signal_names',td_train(1).S1_unit_guide,'disp_times',true);
    % tuning curves
    num_bins = 8;
    tuning_params{1} = struct('num_bins',num_bins,'out_signals',{{'glm_ext_model'}},'out_signal_names',td_train(1).S1_unit_guide);
    tuning_params{2} = struct('num_bins',num_bins,'out_signals',{{'glm_ego_model'}},'out_signal_names',td_train(1).S1_unit_guide);
    tuning_params{3} = struct('num_bins',num_bins,'out_signals',{{'glm_musc_model'}},'out_signal_names',td_train(1).S1_unit_guide);
    tuning_params{4} = struct('num_bins',num_bins,'out_signals',{{'S1_spikes'}},'out_signal_names',td_train(1).S1_unit_guide);
    
    % get PDs and tuning curves
    for modelnum = 1:4
        for spacenum = 1:2
            pdTables{spacenum,modelnum} = getTDPDs(td_test{spacenum},pd_params{modelnum});
            tuning_curves{spacenum,modelnum} = getTuningCurves(td_test{spacenum},tuning_params{modelnum});
        end
        isTuned{modelnum} = checkIsTuned(pdTables{1,modelnum}) & checkIsTuned(pdTables{2,modelnum});
    end
    
    % clean up
    clearvars num_boots num_bins modelnum spacenum

%% Bootstrap on PD shifts
    model_names = {'glm_ext_model','glm_ego_model','glm_musc_model','S1_spikes'};
    shift_tables = cell(length(model_names),1);
    num_internal_boots = 1;
    num_outer_boots = 1000;
    trial_idx = randi(length(td_test{1}),length(td_test{1}),num_outer_boots);
    tic
    for bootctr = 1:num_outer_boots
        disp(['Bootstrap sample ' num2str(bootctr) ', starting at ' num2str(toc) 's'])
        for modelctr = 1:length(model_names)
            pd_params = struct('num_boots',num_internal_boots,'out_signals',{model_names(modelctr)},'out_signal_names',td_train(1).S1_unit_guide);

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
    end

    % clean up
    clearvars num_*_boots trial_idx bootctr modelctr pd_params *_pdTable temp_shift_table

%% Look at dot products of model parameters from simulated neural tuning to actual neural tuning

%% Plot handle positions
    figure
    pos_dl = cat(1,td_test{2}.pos);
    plot(pos_dl(:,1),pos_dl(:,2),'r')
    hold on
    pos_pm = cat(1,td_test{1}.pos);
    plot(pos_pm(:,1),pos_pm(:,2),'b')
    axis equal

    % clean up
    clearvars pos_*

%% Plot example neuron model fit
    % neuron 1 has decent pseudo R2
    figure
    for neuron_idx = 1:length(td_eval{2,1})
        temp_vel = cat(1,td_test{2}.vel);
        temp_spikes = cat(1,td_test{2}.S1_spikes);
        temp_pred_ext = cat(1,td_test{2}.glm_ext_model);
        temp_pred_musc = cat(1,td_test{2}.glm_musc_model);

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
        waitforbuttonpress
    end

    % clean up
    clearvars neuron_idx temp_* ax1 ax2

%% Plot pR2s against each other
    % av_pR2_ext_dl = mean(td_ext_dl_eval,2);
    % av_pR2_musc_dl = mean(td_musc_dl_eval,2);
    % av_pR2_ext_pm = mean(td_ext_pm_eval,2);
    % av_pR2_musc_pm = mean(td_musc_pm_eval,2);
    av_pR2_ext = mean(td_eval{end,1},2);
    av_pR2_musc = mean(td_eval{end,3},2);
    
    % good_neurons = td_ext_dl_eval(:,1) > 0 & td_ext_pm_eval(:,1) > 0 & td_musc_pm_eval(:,1) > 0 & td_musc_dl_eval(:,1) > 0;
    % 
    % figure
    % plot(av_pR2_ext_dl(good_neurons),av_pR2_musc_dl(good_neurons),'ro','linewidth',2)
    % hold on
    % plot([-1 1],[-1 1],'k--','linewidth',2)
    % plot([0 0],[-1 1],'k-','linewidth',2)
    % plot([-1 1],[0 0],'k-','linewidth',2)
    % set(gca,'box','off','tickdir','out','xlim',[-0.1 0.5],'ylim',[-0.1 0.5])
    % 
    % figure
    % plot(av_pR2_ext_pm(good_neurons),av_pR2_musc_pm(good_neurons),'bo','linewidth',2)
    % hold on
    % plot([-1 1],[-1 1],'k--','linewidth',2)
    % plot([0 0],[-1 1],'k-','linewidth',2)
    % plot([-1 1],[0 0],'k-','linewidth',2)
    % set(gca,'box','off','tickdir','out','xlim',[-0.1 0.5],'ylim',[-0.1 0.5])

    % good_neurons = td_ext_eval(:,1) > 0 & td_musc_eval(:,1) > 0;
    good_neurons = isTuned{4};
    figure
    plot(repmat(av_pR2_ext(good_neurons)',2,1),td_eval{end,3}(good_neurons,:)','b-','linewidth',2)
    hold on
    plot(td_eval{end,1}(good_neurons,:)',repmat(av_pR2_musc(good_neurons)',2,1),'b-','linewidth',2)
    plot([-1 1],[-1 1],'k--','linewidth',2)
    plot([0 0],[-1 1],'k-','linewidth',2)
    plot([-1 1],[0 0],'k-','linewidth',2)
    set(gca,'box','off','tickdir','out','xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
    xlabel 'Hand-based pR2'
    ylabel 'Muscle-based pR2'

    % clean up
    clearvars av_pR2* good_neurons

%% Plot comparison of actual tuning curves with various modeled tuning curves
    % first compare PM and DL tuning for each model
    for modelnum = 1:4
        figure;compareTuning(tuning_curves(:,modelnum),pdTables(:,modelnum),find(isTuned{4}))
    end

    % then compare PM and DL tuning for each model
    % reorder for color consistency..
    for spacenum = 1:2
        figure;compareTuning(tuning_curves(spacenum,[3,1,2,4]),pdTables(spacenum,[3,1,2,4]),find(isTuned{4}))
    end

%% Plot DL vs PM just for neurons actually tuned to velocity
    figure
    comparePDs(pm_real_pdTable(isTuned_real,:),dl_real_pdTable(isTuned_real,:),struct('move_corr','vel'),'ko','linewidth',2)
    hold on
    comparePDs(pm_ext_pdTable(isTuned_real,:),dl_ext_pdTable(isTuned_real,:),struct('move_corr','vel'),'ro','linewidth',2)
    comparePDs(pm_ego_pdTable(isTuned_real,:),dl_ego_pdTable(isTuned_real,:),struct('move_corr','vel'),'go','linewidth',2)
    comparePDs(pm_musc_pdTable(isTuned_real,:),dl_musc_pdTable(isTuned_real,:),struct('move_corr','vel'),'bo','linewidth',2)
    xlabel 'PM preferred direction'
    ylabel 'DL preferred direction'
    
%% Plot PD shift clouds for each neuron individually
    colors = {'r','g','b'};
    titles = {'Hand-based model PD shift vs Actual PD shift','Egocentric model PD shift vs Actual PD shift','Muscle-based model PD shift vs Actual PD shift'};
    for modelnum = 1:3
        comparePDClouds(shift_tables{4},shift_tables{modelnum},struct('filter_tuning',1),colors{modelnum},'facealpha',0.5)
        xlabel 'Actual PD Shift'
        ylabel 'Modeled PD Shift'
        title(titles{modelnum})
    end

    % clean up
    clearvars colors titles

