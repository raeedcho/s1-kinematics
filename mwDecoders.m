%% decoder analysis - predict hand movement and arm movement from spikes
%% Setup
    % prep trial data by getting only rewards and trimming to only movements
    % first process marker data
    td = trial_data;
    td = smoothSignals(td,struct('signals','markers'));
    td = getDifferential(td,struct('signal','markers','alias','marker_vel'));

    % bin data at 50ms
    td = binTD(td,5);
    % add in spherical coordinates
    td = addSphereHand2TD(td);

    % duplicate and shift spike data before chopping data up
    % Get future 150 ms of S1 activity to predict current kinematics
    % non-overlapping bins...
    td = dupeAndShift(td,'S1_spikes',-3);

    % get only rewards
    [~,td] = getTDidx(td,'result','R');

    % only go from first go cue to trial end
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    % for bumps
    % [~,td] = getTDidx(trial_data,'result','R');
    % td = td(~isnan(cat(1,td.idx_bumpTime)));
    % td = trimTD(td,{'idx_bumpTime',0},{'idx_bumpTime',15});

    % clean up trial data by taking out bad trials
    pause;

    % Split td into different workspaces (workspace 1 is PM and workspace 2 is DL)
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    [~,td_pm] = getTDidx(td,'spaceNum',1);
    [~,td_dl] = getTDidx(td,'spaceNum',2);
    minsize = min(length(td_pm),length(td_dl));
    td_pm = td_pm(1:minsize);
    td_dl = td_dl(1:minsize);

    % recombine for later...
    td = [td_pm td_dl];

    % %% Do PCA on muscle space
    % % do PCA on muscles, training on only the training set
    % % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % % PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    % %                     'do_plot',true);
    % PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len'))}}, 'do_plot',true);
    % [td,~] = getPCA(td,PCAparams);
    % % temporary hack to allow us to do PCA on velocity too
    % for i=1:length(td)
    %     td(i).opensim_len_pca = td(i).opensim_pca;
    % end
    % % get rid of superfluous PCA
    % td = rmfield(td,'opensim_pca');
    % % get velocity PCA
    % % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
    % % PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
    % %                     'do_plot',true);
    % PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel'))}}, 'do_plot',true);
    % [td_temp,~] = getPCA(td,PCAparams_vel);
    % % temporary hack to allow us to save into something useful
    % for i=1:length(td)
    %     td(i).opensim_muscVel_pca = td(i).opensim_pca;
    % end
    % % get rid of superfluous PCA
    % td = rmfield(td,'opensim_pca');

    % Get PCA for neural space
    % PCAparams = struct('signals',{{'S1_spikes'}}, 'do_plot',true,'pca_recenter_for_proj',true,'sqrt_transform',true);
    % [td,~] = getPCA(td,PCAparams);

    % Get PCA for marker space
    td = getPCA(td,struct('signals','markers'));
    td = getPCA(td,struct('signals','marker_vel'));

%% train decoders and evaluate
    % split into folds
    num_repeats = 10;
    num_folds = 10;

    % preallocate vaf holders
    neur_decoder_vaf = zeros(num_repeats,num_folds,2);
    hand_decoder_vaf = zeros(num_repeats,num_folds,2);

    % reverse model, predicing hand from elbow
    reverse_neur_decoder_vaf = zeros(num_repeats,num_folds,2);
    reverse_hand_decoder_vaf = zeros(num_repeats,num_folds,2);

    % another control, predicting one hand marker from another
    other_neur_decoder_vaf = zeros(num_repeats,num_folds,2);
    other_hand_decoder_vaf = zeros(num_repeats,num_folds,2);

    % set model parameters
    % to predict elbow from both neurons and hand
    elbow_idx = 28:30;
    hand_idx = 1:3;
    other_hand_idx = 4:6;

    % decoder parameters
    neural_signals = {'S1_spikes','all';'S1_spikes_shift','all'};
    elbow_signals = {'markers',elbow_idx;'marker_vel',elbow_idx};
    hand_signals = {'markers',hand_idx;'marker_vel',hand_idx};
    other_hand_signals = {'markers',other_hand_idx;'marker_vel',other_hand_idx};

    % predict elbow from either just hand or hand and neurons
    neur_decoder_params = struct('model_type','linmodel','model_name','neur_decoder',...
        'in_signals',{[neural_signals;hand_signals]},...
        'out_signals',{elbow_signals});
    hand_decoder_params = struct('model_type','linmodel','model_name','hand_decoder',...
        'in_signals',{hand_signals},...
        'out_signals',{elbow_signals});
    % predict hand from either elbow or hand and elbow
    reverse_neur_decoder_params = struct('model_type','linmodel','model_name','reverse_neur_decoder',...
        'in_signals',{[neural_signals;elbow_signals]},...
        'out_signals',{hand_signals});
    reverse_hand_decoder_params = struct('model_type','linmodel','model_name','reverse_hand_decoder',...
        'in_signals',{elbow_signals},...
        'out_signals',{hand_signals});
    % predict hand from either hand or hand and neurons
    other_neur_decoder_params = struct('model_type','linmodel','model_name','other_neur_decoder',...
        'in_signals',{[neural_signals;hand_signals]},...
        'out_signals',{other_hand_signals});
    other_hand_decoder_params = struct('model_type','linmodel','model_name','other_hand_decoder',...
        'in_signals',{hand_signals},...
        'out_signals',{other_hand_signals});
    % define helpful function...
    square_sum = @(x) sum(x.^2,1)*[1 1 1 0 0 0; 0 0 0 1 1 1]';

    % do crossval
    for repeatnum = 1:num_repeats
        inds = crossvalind('kfold',length(td),num_folds);

        for foldnum = 1:num_folds
            % fit models
            [~,neur_decoder] = getModel(td(inds~=foldnum),neur_decoder_params);
            [~,hand_decoder] = getModel(td(inds~=foldnum),hand_decoder_params);
            [~,reverse_neur_decoder] = getModel(td(inds~=foldnum),reverse_neur_decoder_params);
            [~,reverse_hand_decoder] = getModel(td(inds~=foldnum),reverse_hand_decoder_params);
            [~,other_neur_decoder] = getModel(td(inds~=foldnum),other_neur_decoder_params);
            [~,other_hand_decoder] = getModel(td(inds~=foldnum),other_hand_decoder_params);

            % evaluate models and save into array
            td_test = td(inds==foldnum);
            td_test = getModel(td_test,neur_decoder);
            td_test = getModel(td_test,hand_decoder);
            td_test = getModel(td_test,reverse_neur_decoder);
            td_test = getModel(td_test,reverse_hand_decoder);
            td_test = getModel(td_test,other_neur_decoder);
            td_test = getModel(td_test,other_hand_decoder);

            % get error on elbow
            elbow_true = get_vars(td_test,{'markers',elbow_idx;'marker_vel',elbow_idx});
            elbow_mean = mean(elbow_true);
            elbow_pred_neur = get_vars(td_test,check_signals(td_test,'linmodel_neur_decoder'));
            elbow_pred_hand = get_vars(td_test,check_signals(td_test,'linmodel_hand_decoder'));

            % calculate fraction of variance explained
            SS_total = square_sum(elbow_true-elbow_mean);
            neur_decoder_vaf(repeatnum,foldnum,:) = 1 - square_sum(elbow_pred_neur-elbow_true)./SS_total;
            hand_decoder_vaf(repeatnum,foldnum,:) = 1 - square_sum(elbow_pred_hand-elbow_true)./SS_total;

            % reverse models -- get error on hand
            hand_true = get_vars(td_test,{'markers',hand_idx;'marker_vel',hand_idx});
            hand_mean = mean(hand_true);
            hand_pred_neur = get_vars(td_test,check_signals(td_test,'linmodel_reverse_neur_decoder'));
            hand_pred_hand = get_vars(td_test,check_signals(td_test,'linmodel_reverse_hand_decoder'));

            SS_total = square_sum(hand_true-hand_mean);
            reverse_neur_decoder_vaf(repeatnum,foldnum,:) = 1 - square_sum(hand_pred_neur-hand_true)./SS_total;
            reverse_hand_decoder_vaf(repeatnum,foldnum,:) = 1 - square_sum(hand_pred_hand-hand_true)./SS_total;

            % as control, get error on other hand marker from hand
            other_hand_true = get_vars(td_test,{'markers',other_hand_idx;'marker_vel',other_hand_idx});
            other_hand_mean = mean(other_hand_true);
            other_hand_pred_neur = get_vars(td_test,check_signals(td_test,'linmodel_other_neur_decoder'));
            other_hand_pred_hand = get_vars(td_test,check_signals(td_test,'linmodel_other_hand_decoder'));

            % calculate fraction of variance explained
            SS_total = square_sum(other_hand_true-other_hand_mean);
            other_neur_decoder_vaf(repeatnum,foldnum,:) = 1 - square_sum(other_hand_pred_neur-other_hand_true)./SS_total;
            other_hand_decoder_vaf(repeatnum,foldnum,:) = 1 - square_sum(other_hand_pred_hand-other_hand_true)./SS_total;
        end
    end
    
    neur_decoder_vaf = reshape(neur_decoder_vaf,num_repeats*num_folds,2);
    hand_decoder_vaf = reshape(hand_decoder_vaf,num_repeats*num_folds,2);
    reverse_neur_decoder_vaf = reshape(reverse_neur_decoder_vaf,num_repeats*num_folds,2);
    reverse_hand_decoder_vaf = reshape(reverse_hand_decoder_vaf,num_repeats*num_folds,2);
    other_neur_decoder_vaf = reshape(other_neur_decoder_vaf,num_repeats*num_folds,2);
    other_hand_decoder_vaf = reshape(other_hand_decoder_vaf,num_repeats*num_folds,2);
    
    vaf_diff = neur_decoder_vaf-hand_decoder_vaf;
    mean_vaf = mean(vaf_diff);
    correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
    % upp = tinv(0.975,99);
    low = tinv(0.01,num_folds*num_repeats-1);
    % vafcihi = mean_vaf + upp * sqrt(correction*vardiff);
    vafcilo = mean_vaf + low * sqrt(correction*var(vaf_diff));
    p_val = tcdf(-mean_vaf./sqrt(correction*var(vaf_diff)),num_folds*num_repeats-1);

    % plot vafs
    figure
    bar([1 3 2 4]*0.20,mean([hand_decoder_vaf neur_decoder_vaf]),0.20,'facecolor',[0.5 0.5 0.5],'edgecolor','none')
    hold on
    plot([1 2]*0.20,[hand_decoder_vaf(:,1) neur_decoder_vaf(:,1)]','-k','linewidth',1)
    plot([3 4]*0.20,[hand_decoder_vaf(:,2) neur_decoder_vaf(:,2)]','-k','linewidth',1)
    plot([1 3 2 4]*0.20, [hand_decoder_vaf neur_decoder_vaf]','.k','markersize',30)
    % set(gca,'box','off','tickdir','out','xtick',1:4,...
    %     'xticklabel',{'Hand-only pos','Hand+Neuron pos','Hand vel','Hand+Neuron vel'},...
    %     'xlim',[0 5],'ylim',[0 1])
    % ylabel 'Fraction VAF'
    % xlabel 'Model'
    % title 'Elbow decoding model performance'

    % plot reverse model
    figure
    bar([1 3 2 4]*0.20+1,mean([reverse_hand_decoder_vaf reverse_neur_decoder_vaf]),0.20,'facecolor',[0.5 0.5 0.5],'edgecolor','none')
    hold on
    plot([1 2]*0.20+1,[reverse_hand_decoder_vaf(:,1) reverse_neur_decoder_vaf(:,1)]','-k','linewidth',1)
    plot([3 4]*0.20+1,[reverse_hand_decoder_vaf(:,2) reverse_neur_decoder_vaf(:,2)]','-k','linewidth',1)
    plot([1 3 2 4]*0.20+1, [reverse_hand_decoder_vaf reverse_neur_decoder_vaf]','.k','markersize',30)
    % set(gca,'box','off','tickdir','out','xtick',1:4,...
    %     'xticklabel',{'Elbow-only pos','Elbow+Neuron pos','Elbow vel','Elbow+Neuron vel'},...
    %     'xlim',[0 5],'ylim',[0 1])
    % ylabel 'Fraction VAF'
    % xlabel 'Model'
    % title 'Decoding model performance'

    % plot control
    figure
    bar([1 3 2 4]*0.20+2,mean([other_hand_decoder_vaf other_neur_decoder_vaf]),0.20,'facecolor',[0.5 0.5 0.5],'edgecolor','none')
    hold on
    plot([1 2]*0.20+2,[other_hand_decoder_vaf(:,1) other_neur_decoder_vaf(:,1)]','-k','linewidth',1)
    plot([3 4]*0.20+2,[other_hand_decoder_vaf(:,2) other_neur_decoder_vaf(:,2)]','-k','linewidth',1)
    plot([1 3 2 4]*0.20+2, [other_hand_decoder_vaf other_neur_decoder_vaf]','.k','markersize',30)
    % set(gca,'box','off','tickdir','out','xtick',1:4,...
    %     'xticklabel',{'Hand-only pos','Hand+Neuron pos','Hand vel','Hand+Neuron vel'},...
    %     'xlim',[0 5],'ylim',[0 1])
    % ylabel 'Fraction VAF'
    % xlabel 'Model'
    % title 'Control decoding model performance'

    % plot niceness
    xtick = ones(3,1)*(1:4)*0.20 + (0:2)';
    xticklabel = {'Hand-only pos','Hand+Neuron pos','Hand vel','Hand+Neuron vel',...
        'Elbow-only pos','Elbow+Neuron pos','Elbow vel','Elbow+Neuron vel',...
        'HandControl-only pos','HandControl+Neuron pos','HandControl vel','HandControl+Neuron vel'};
    set(gca,'box','off','tickdir','out','xtick',xtick,...
        'xticklabel',xticklabel,...
        'xlim',[0 3],'ylim',[0 1])

    % plot example
    % for i = 1:length(td_test)
    %     f = figure;

    %     % plot3(td_test(1).markers(:,elbow_idx(1)),td_test(1).markers(:,elbow_idx(2)),td_test(1).markers(:,elbow_idx(3)),'-k','linewidth',3)
    %     % hold on
    %     % plot3(td_test(1).linmodel_neur_decoder(:,(1)),td_test(1).linmodel_neur_decoder(:,(2)),td_test(1).linmodel_neur_decoder(:,(3)),'-b','linewidth',2)
    %     % plot3(td_test(1).linmodel_hand_decoder(:,(1)),td_test(1).linmodel_hand_decoder(:,(2)),td_test(1).linmodel_hand_decoder(:,(3)),'-r','linewidth',2)
    %     % axis equal
    %     subplot(2,1,1)
    %     plot(td_test(i).markers(:,elbow_idx),'k-','linewidth',3)
    %     hold on
    %     plot(td_test(i).linmodel_neur_decoder(:,1:3),'b-','linewidth',1)
    %     plot(td_test(i).linmodel_hand_decoder(:,1:3),'r-','linewidth',1)

    %     subplot(2,1,2)
    %     plot(td_test(i).marker_vel(:,elbow_idx),'k-','linewidth',3)
    %     hold on
    %     plot(td_test(i).linmodel_neur_decoder(:,4:6),'b-','linewidth',1)
    %     plot(td_test(i).linmodel_hand_decoder(:,4:6),'r-','linewidth',1)


    %     waitfor(f)
    % end

