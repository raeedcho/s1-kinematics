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
    td_pm_test = td_pm(test_idx);
    td_dl_test = td_dl(test_idx);

%% Train models on training set composed of both workspaces together
    % do PCA on muscles, training on only the training set
    % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    %                     'do_plot',true);
    PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len'))}}, 'do_plot',true);
    [td_train,pca_info] = getPCA(td_train,PCAparams);
    td_pm_test = getPCA(td_pm_test,pca_info);
    td_dl_test = getPCA(td_dl_test,pca_info);
    
    % temporary hack to allow us to do PCA on velocity too
    for i=1:length(td_train)
        td_train(i).opensim_len_pca = td_train(i).opensim_pca;
    end
    for i=1:length(td_pm_test)
        td_pm_test(i).opensim_len_pca = td_pm_test(i).opensim_pca;
        td_dl_test(i).opensim_len_pca = td_dl_test(i).opensim_pca;
    end
    
    % get rid of superfluous PCA
    td_train = rmfield(td_train,'opensim_pca');
    td_pm_test = rmfield(td_pm_test,'opensim_pca');
    td_dl_test = rmfield(td_dl_test,'opensim_pca');

    % get velocity PCA
    % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    %                     'do_plot',true);
    PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel'))}}, 'do_plot',true);
    [td_train,pca_info_vel] = getPCA(td_train,PCAparams_vel);
    td_pm_test = getPCA(td_pm_test,pca_info_vel);
    td_dl_test = getPCA(td_dl_test,pca_info_vel);
    
    % temporary hack to allow us to save into something useful
    for i=1:length(td_train)
        td_train(i).opensim_muscVel_pca = td_train(i).opensim_pca;
    end
    for i=1:length(td_pm_test)
        td_pm_test(i).opensim_muscVel_pca = td_pm_test(i).opensim_pca;
        td_dl_test(i).opensim_muscVel_pca = td_dl_test(i).opensim_pca;
    end
    
    % get rid of superfluous PCA
    td_train = rmfield(td_train,'opensim_pca');
    td_pm_test = rmfield(td_pm_test,'opensim_pca');
    td_dl_test = rmfield(td_dl_test,'opensim_pca');

%% Plot handle positions
    figure
    pos_dl = cat(1,td_dl_test.pos);
    plot(pos_dl(:,1),pos_dl(:,2),'r')
    mean(pos_dl)
    hold on
    pos_pm = cat(1,td_pm_test.pos);
    plot(pos_pm(:,1),pos_pm(:,2),'b')
    axis equal
    clear pos_*

%% Fit models
    % set up parameters for models
    % cartesian hand coordinates for position and velocity
    opensim_hand_idx = find(contains(td(1).opensim_names,'_handPos') | contains(td(1).opensim_names,'_handVel'));
    glm_ext_params = struct('model_type','glm',...
                            'model_name','ext_model',...
                            'in_signals',{{'opensim',opensim_hand_idx}},...
                            'out_signals',{{'S1_spikes'}});
    glm_ego_params = struct('model_type','glm',...
                            'model_name','ego_model',...
                            'in_signals',{{'sphere_hand_pos';'sphere_hand_vel'}},...
                            'out_signals',{{'S1_spikes'}});
    glm_musc_params = struct('model_type','glm',...
                            'model_name','musc_model',...
                            'in_signals',{{'opensim_len_pca',1:5;'opensim_muscVel_pca',1:5}},...
                            'out_signals',{'S1_spikes'});
    [~,glm_info_ext] = getModel(td_train,glm_ext_params);
    [~,glm_info_ego] = getModel(td_train,glm_ego_params);
    [~,glm_info_musc] = getModel(td_train,glm_musc_params);

%% Predict firing rates
    td_dl_test = getModel(td_dl_test,glm_info_ext);
    td_pm_test = getModel(td_pm_test,glm_info_ext);
    td_dl_test = getModel(td_dl_test,glm_info_ego);
    td_pm_test = getModel(td_pm_test,glm_info_ego);
    td_dl_test = getModel(td_dl_test,glm_info_musc);
    td_pm_test = getModel(td_pm_test,glm_info_musc);

%% Get PDs for the modeled and actual neurons
    num_boots = 100;
    pd_ext_params = struct('num_boots',num_boots,'out_signals',{{'glm_ext_model'}},'out_signal_names',td(1).S1_unit_guide,'disp_times',true);
    pd_ego_params = struct('num_boots',num_boots,'out_signals',{{'glm_ego_model'}},'out_signal_names',td(1).S1_unit_guide,'disp_times',true);
    pd_musc_params = struct('num_boots',num_boots,'out_signals',{{'glm_musc_model'}},'out_signal_names',td(1).S1_unit_guide,'disp_times',true);
    pd_real_params = struct('num_boots',num_boots,'out_signals',{{'S1_spikes'}},'out_signal_names',td(1).S1_unit_guide,'disp_times',true);
    pm_ext_pdTable = getTDPDs(td_pm_test,pd_ext_params);
    pm_ego_pdTable = getTDPDs(td_pm_test,pd_ego_params);
    pm_musc_pdTable = getTDPDs(td_pm_test,pd_musc_params);
    pm_real_pdTable = getTDPDs(td_pm_test,pd_real_params);
    dl_ext_pdTable = getTDPDs(td_dl_test,pd_ext_params);
    dl_ego_pdTable = getTDPDs(td_dl_test,pd_ego_params);
    dl_musc_pdTable = getTDPDs(td_dl_test,pd_musc_params);
    dl_real_pdTable = getTDPDs(td_dl_test,pd_real_params);

    % get PD shifts
    ext_pdShift = minusPi2Pi(dl_ext_pdTable.velPD-pm_ext_pdTable.velPD);
    ego_pdShift = minusPi2Pi(dl_ego_pdTable.velPD-pm_ego_pdTable.velPD);
    musc_pdShift = minusPi2Pi(dl_musc_pdTable.velPD-pm_musc_pdTable.velPD);
    real_pdShift = minusPi2Pi(dl_real_pdTable.velPD-pm_real_pdTable.velPD);

    % figure out which neurons are tuned
    isTuned_ext = checkIsTuned(pm_ext_pdTable) & checkIsTuned(dl_ext_pdTable);
    isTuned_ego = checkIsTuned(pm_ego_pdTable) & checkIsTuned(dl_ego_pdTable);
    isTuned_musc = checkIsTuned(pm_musc_pdTable) & checkIsTuned(dl_musc_pdTable);
    isTuned_real = checkIsTuned(pm_real_pdTable) & checkIsTuned(dl_real_pdTable);

%% Bootstrap on PD shifts
    model_names = {'glm_ext_model','glm_ego_model','glm_musc_model','S1_spikes'};
    num_internal_boots = 1;
    num_outer_boots = 1000;
    shift_tables = cell(length(model_names),1);
    trial_idx = randi(length(td_pm_test),length(td_pm_test),num_outer_boots);
    tic
    for bootctr = 1:num_outer_boots
        disp(['Bootstrap sample ' num2str(bootctr) ', starting at ' num2str(toc) 's'])
        for modelctr = 1:length(model_names)
            pd_params = struct('num_boots',num_internal_boots,'out_signals',{model_names(modelctr)},'out_signal_names',td(1).S1_unit_guide);

            pm_pdTable = getTDPDs(td_pm_test(trial_idx(:,bootctr)),pd_params);
            dl_pdTable = getTDPDs(td_dl_test(trial_idx(:,bootctr)),pd_params);

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

%% Plot DL vs PM just for neurons actually tuned to velocity
    figure
    comparePDs(pm_real_pdTable(isTuned_real,:),dl_real_pdTable(isTuned_real,:),struct('move_corr','vel'),'ko','linewidth',2)
    hold on
    comparePDs(pm_ext_pdTable(isTuned_real,:),dl_ext_pdTable(isTuned_real,:),struct('move_corr','vel'),'ro','linewidth',2)
    comparePDs(pm_ego_pdTable(isTuned_real,:),dl_ego_pdTable(isTuned_real,:),struct('move_corr','vel'),'go','linewidth',2)
    comparePDs(pm_musc_pdTable(isTuned_real,:),dl_musc_pdTable(isTuned_real,:),struct('move_corr','vel'),'bo','linewidth',2)
    xlabel 'PM preferred direction'
    ylabel 'DL preferred direction'
    
%% Plot PDs against real PD shift (one point per bootstrap shift)
    isTuned_rep = repmat(isTuned_real,num_outer_boots,1);
    figure
    comparePDs(shift_tables{4}(isTuned_rep,:),shift_tables{1}(isTuned_rep,:),struct('move_corr','vel'),'ro','linewidth',2)
    hold on
    comparePDs(shift_tables{4}(isTuned_rep,:),shift_tables{2}(isTuned_rep,:),struct('move_corr','vel'),'go','linewidth',2)
    comparePDs(shift_tables{4}(isTuned_rep,:),shift_tables{3}(isTuned_rep,:),struct('move_corr','vel'),'bo','linewidth',2)
    xlabel 'Actual PD shift'
    ylabel 'Modeled PD shift'
    title 'One point per bootstrap sample'

%% Plot PD shift clouds for each neuron individually
    % Loop through and get mean shifts out of bootstrapped shifts
    signalIDs = unique(shift_tables{1}.signalID,'rows');
    % only look at tuned neurons
    signalIDs = signalIDs(isTuned_real,:);
    model_cstrings = {'r','g','b'}
    model_figures = {figure;figure;figure};
    model_titles = {'Hand-based','Egocentric','Muscle-based'}
    for modelctr = 1:length(shift_tables)-1
        % set up figure axes
        figure(model_figures{modelctr})
        plot([-pi pi],[0 0],'-k','linewidth',2)
        hold on
        plot([0 0],[-pi pi],'-k','linewidth',2)
        plot([-pi pi],[-pi pi],'--k','linewidth',2)
        axis equal
        set(gca,'box','off','tickdir','out','xtick',[-pi pi],'ytick',[-pi pi],'xlim',[-pi pi],'ylim',[-pi pi],...
            'xticklabel',{'-\pi','\pi'},'yticklabel',{'-\pi','\pi'})
        xlabel 'Actual PD Shift'
        ylabel 'Modeled PD Shift'
        title([model_titles{modelctr} ' model PD shift vs Actual PD shift'])
    end
    for i = 1:size(signalIDs,1)
        % make new table set
        shift_tables_unit = cell(size(shift_tables));
        for modelctr = 1:length(shift_tables)
            % Populate new table with only that unit's PD shifts
            unit_idx = ismember(shift_tables{modelctr}.signalID,signalIDs(i,:),'rows');
            shift_tables_unit{modelctr} = shift_tables{modelctr}(unit_idx,:);
        end
        % now create hulls and plot
        for modelctr = 1:length(shift_tables)-1
            figure(model_figures{modelctr})

            % first plot cloud with center point
            %comparePDs(shift_tables_unit{4},shift_tables_unit{modelctr},struct('move_corr','vel'),[model_cstrings{modelctr} 'o'],'linewidth',2)
            %plot(circ_mean(shift_tables_unit{4}.velPD),circ_mean(shift_tables_unit{modelctr}.velPD),'ko','linewidth',2);

            % figure out 95% confidence interval by 'euclidean' distance to circular mean
            % should I be calculating geodesic distance somehow? seems like no, based on a preliminary google search
            clust = [shift_tables_unit{4}.velPD shift_tables_unit{modelctr}.velPD];
            means = [circ_mean(shift_tables_unit{4}.velPD) circ_mean(shift_tables_unit{modelctr}.velPD)];
            centered_clust = minusPi2Pi(clust-repmat(means,num_outer_boots,1));
            dists = sqrt(sum(centered_clust.^2,2));
            inliers = dists<prctile(dists,95);
            clust = clust(inliers,:);
            centered_clust = centered_clust(inliers,:);

            % generate convex hull of inliers
            hull_idx = convhull(centered_clust);

            % plot cloud
            % plot hull (have to figure out what to do about wraparound)
            patch(centered_clust(hull_idx,1)+means(1),centered_clust(hull_idx,2)+means(2),model_cstrings{modelctr},'facealpha',0.5);
            patch(centered_clust(hull_idx,1)+means(1)-2*pi,centered_clust(hull_idx,2)+means(2),model_cstrings{modelctr},'facealpha',0.5);
            patch(centered_clust(hull_idx,1)+means(1),centered_clust(hull_idx,2)+means(2)-2*pi,model_cstrings{modelctr},'facealpha',0.5);
            patch(centered_clust(hull_idx,1)+means(1)-2*pi,centered_clust(hull_idx,2)+means(2)-2*pi,model_cstrings{modelctr},'facealpha',0.5);
        end
    end

%% Plot PDs against each other (one point per mean PD shift)
    figure
    % axes
    plot([-pi pi],[0 0],'-k','linewidth',2)
    hold on
    plot([0 0],[-pi pi],'-k','linewidth',2)
    % unity line
    plot([-pi pi],[-pi pi],'--k','linewidth',2)
    % first plot ext vs. real
    plot(real_pdShift(isTuned_real),ext_pdShift(isTuned_real),'ro','linewidth',3,'markersize',5)
    % then ego vs real
    plot(real_pdShift(isTuned_real),ego_pdShift(isTuned_real),'go','linewidth',3,'markersize',5)
    % then musc vs real
    plot(real_pdShift(isTuned_real),musc_pdShift(isTuned_real),'bo','linewidth',3,'markersize',5)
    set(gca,'box','off','tickdir','out','xlim',[-pi pi],'ylim',[-pi pi],'xtick',[-pi 0 pi],'ytick',[-pi 0 pi])
    axis equal
    xlabel 'Actual PD shift'
    ylabel 'Modeled PD shift'
    title 'One point per mean PD shift'

%% Evaluate model fits
    eval_params = glm_info_ext;
    eval_params.eval_metric = 'pr2';
    td_ext_dl_eval = squeeze(evalModel(td_dl_test,eval_params));
    td_ext_pm_eval = squeeze(evalModel(td_pm_test,eval_params));
    td_ext_eval = squeeze(evalModel([td_dl_test td_pm_test],eval_params))
    
    eval_params = glm_info_ego;
    eval_params.eval_metric = 'pr2';
    td_ego_dl_eval = squeeze(evalModel(td_dl_test,eval_params));
    td_ego_pm_eval = squeeze(evalModel(td_pm_test,eval_params));
    td_ego_eval = squeeze(evalModel([td_dl_test td_pm_test],eval_params))
    
    eval_params = glm_info_musc;
    eval_params.eval_metric = 'pr2';
    td_musc_dl_eval = squeeze(evalModel(td_dl_test,eval_params));
    td_musc_pm_eval = squeeze(evalModel(td_pm_test,eval_params));
    td_musc_eval = squeeze(evalModel([td_dl_test td_pm_test],eval_params))

%% Plot example
    % neuron 1 has decent pseudo R2
    figure
    for neuron_idx = 1:length(td_ext_dl_eval)
        temp_vel = cat(1,td_dl_test.vel);
        temp_spikes = cat(1,td_dl_test.S1_spikes);
        temp_pred_ext = cat(1,td_dl_test.glm_ext_model);
        temp_pred_musc = cat(1,td_dl_test.glm_musc_model);

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

%% Plot pR2s against each other
    % av_pR2_ext_dl = mean(td_ext_dl_eval,2);
    % av_pR2_musc_dl = mean(td_musc_dl_eval,2);
    % av_pR2_ext_pm = mean(td_ext_pm_eval,2);
    % av_pR2_musc_pm = mean(td_musc_pm_eval,2);
    av_pR2_ext = mean(td_ext_eval,2);
    av_pR2_musc = mean(td_musc_eval,2);
    
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
    good_neurons = isTuned_real
    figure
    plot(repmat(av_pR2_ext(good_neurons)',2,1),td_musc_eval(good_neurons,:)','b-','linewidth',2)
    hold on
    plot(td_ext_eval(good_neurons,:)',repmat(av_pR2_musc(good_neurons)',2,1),'b-','linewidth',2)
    plot([-1 1],[-1 1],'k--','linewidth',2)
    plot([0 0],[-1 1],'k-','linewidth',2)
    plot([-1 1],[0 0],'k-','linewidth',2)
    set(gca,'box','off','tickdir','out','xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
    xlabel 'Hand-based pR2'
    ylabel 'Muscle-based pR2'
