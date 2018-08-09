%% decoder analysis - predict hand movement and arm movement from spikes
%% Setup
    % Load data
    datadir = '/home/raeed/Projects/limblab/data-td/MultiWorkspace';
    load(sprintf('%s/Han_20171101_TD.mat',datadir))

    % prep trial data by getting only rewards and trimming to only movements
    % first process marker data
    td = trial_data;
    td = smoothSignals(td,struct('signals','markers'));
    td = getDifferential(td,struct('signal','markers','alias','marker_vel'));

    % get only rewards
    [~,td] = getTDidx(td,'result','R');

    % Split td into different workspaces (workspace 1 is PM and workspace 2 is DL)
    % also make sure we have balanced workspaces (slightly biases us towards early trials, but this isn't too bad)
    [~,td_pm] = getTDidx(td,'spaceNum',1);
    [~,td_dl] = getTDidx(td,'spaceNum',2);
    minsize = min(length(td_pm),length(td_dl));
    td_pm = td_pm(1:minsize);
    td_dl = td_dl(1:minsize);

    % recombine for later...
    td = [td_pm td_dl];

%% Set up helper variables
    % to predict elbow from both neurons and hand
    elbow_idx = 28:30;
    hand_idx = 4:6;
    other_hand_idx = 7:9;
    wrist1_idx = 10:12;
    wrist2_idx = 13:15;

%% Simulate neurons that depend on only hand and elbow
    num_sim_neurons = 50;
    weights = randn(13,num_sim_neurons);

    kinem = get_vars(td,{'markers',[hand_idx elbow_idx];'marker_vel',[hand_idx elbow_idx]});
    mean_kinem = mean(kinem);
    std_kinem = std(kinem);

    for trial_idx = 1:length(td)
        kinem = get_vars(td(trial_idx),{'markers',[hand_idx elbow_idx];'marker_vel',[hand_idx elbow_idx]});
        z_kinem = (kinem-mean_kinem)./std_kinem;
        simmed_neur = binornd(1,1./(1+exp([ones(length(z_kinem),1) z_kinem]*weights)));
        td(trial_idx).S1_sim = simmed_neur;
    end

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

    % another control, predicting one hand marker from another
    wrist_neur_decoder_vaf = zeros(num_repeats,num_folds,2);
    wrist_decoder_vaf = zeros(num_repeats,num_folds,2);

    % duplicate and shift spike data before chopping data up
    % Get future 150 ms of S1 activity to predict current kinematics
    % non-overlapping bins...
    neur_name = 'S1_spikes';
    td = dupeAndShift(td,neur_name,-3);

    % Trim td
    td = trimTD(td,{'idx_targetStartTime',0},{'idx_endTime',0});

    % set model parameters
    % decoder parameters
    neural_signals = {neur_name,'all';sprintf('%s_shift',neur_name),'all'};
    elbow_signals = {'markers',elbow_idx;'marker_vel',elbow_idx};
    hand_signals = {'markers',hand_idx;'marker_vel',hand_idx};
    other_hand_signals = {'markers',other_hand_idx;'marker_vel',other_hand_idx};
    wrist1_signals = {'markers',wrist1_idx;'marker_vel',wrist1_idx};
    wrist2_signals = {'markers',wrist2_idx;'marker_vel',wrist2_idx};

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
    wrist_neur_decoder_params = struct('model_type','linmodel','model_name','wrist_neur_decoder',...
        'in_signals',{[neural_signals;wrist1_signals]},...
        'out_signals',{wrist2_signals});
    wrist_decoder_params = struct('model_type','linmodel','model_name','wrist_decoder',...
        'in_signals',{wrist1_signals},...
        'out_signals',{wrist2_signals});
    % define helpful function...
    square_sum = @(x) sum(x.^2,1)*[1 1 1 0 0 0; 0 0 0 1 1 1]';

    % do crossval
    repeat_tic = tic;
    fprintf('Starting %dx%d-fold crossvalidation at time %f\n',num_repeats,num_folds,toc(repeat_tic));
    for repeatnum = 1:num_repeats
        fprintf('Starting cross-val repeat number %d at time %f\n',repeatnum,toc(repeat_tic))
        inds = crossvalind('kfold',length(td),num_folds);

        fold_tic = tic;
        for foldnum = 1:num_folds
            % fit models
            [~,neur_decoder] = getModel(td(inds~=foldnum),neur_decoder_params);
            [~,hand_decoder] = getModel(td(inds~=foldnum),hand_decoder_params);
            [~,reverse_neur_decoder] = getModel(td(inds~=foldnum),reverse_neur_decoder_params);
            [~,reverse_hand_decoder] = getModel(td(inds~=foldnum),reverse_hand_decoder_params);
            [~,other_neur_decoder] = getModel(td(inds~=foldnum),other_neur_decoder_params);
            [~,other_hand_decoder] = getModel(td(inds~=foldnum),other_hand_decoder_params);
            [~,wrist_neur_decoder] = getModel(td(inds~=foldnum),wrist_neur_decoder_params);
            [~,wrist_decoder] = getModel(td(inds~=foldnum),wrist_decoder_params);

            % evaluate models and save into array
            td_test = td(inds==foldnum);
            td_test = getModel(td_test,neur_decoder);
            td_test = getModel(td_test,hand_decoder);
            td_test = getModel(td_test,reverse_neur_decoder);
            td_test = getModel(td_test,reverse_hand_decoder);
            td_test = getModel(td_test,other_neur_decoder);
            td_test = getModel(td_test,other_hand_decoder);
            td_test = getModel(td_test,wrist_neur_decoder);
            td_test = getModel(td_test,wrist_decoder);

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

            % as control, get error on other hand marker from hand
            wrist2_true = get_vars(td_test,{'markers',wrist2_idx;'marker_vel',wrist2_idx});
            wrist2_mean = mean(wrist2_true);
            wrist2_pred_neur = get_vars(td_test,check_signals(td_test,'linmodel_wrist_neur_decoder'));
            wrist2_pred_hand = get_vars(td_test,check_signals(td_test,'linmodel_wrist_decoder'));

            % calculate fraction of variance explained
            SS_total = square_sum(wrist2_true-wrist2_mean);
            wrist_neur_decoder_vaf(repeatnum,foldnum,:) = 1 - square_sum(wrist2_pred_neur-wrist2_true)./SS_total;
            wrist_decoder_vaf(repeatnum,foldnum,:) = 1 - square_sum(wrist2_pred_hand-wrist2_true)./SS_total;
            % display
            fprintf('\tEvaluated fold number %d at time %f\n',foldnum,toc(fold_tic))
        end
    end
    
    neur_decoder_vaf = reshape(neur_decoder_vaf,num_repeats*num_folds,2);
    hand_decoder_vaf = reshape(hand_decoder_vaf,num_repeats*num_folds,2);
    reverse_neur_decoder_vaf = reshape(reverse_neur_decoder_vaf,num_repeats*num_folds,2);
    reverse_hand_decoder_vaf = reshape(reverse_hand_decoder_vaf,num_repeats*num_folds,2);
    other_neur_decoder_vaf = reshape(other_neur_decoder_vaf,num_repeats*num_folds,2);
    other_hand_decoder_vaf = reshape(other_hand_decoder_vaf,num_repeats*num_folds,2);
    wrist_neur_decoder_vaf = reshape(wrist_neur_decoder_vaf,num_repeats*num_folds,2);
    wrist_decoder_vaf = reshape(wrist_decoder_vaf,num_repeats*num_folds,2);
    
    vaf_diff = wrist_neur_decoder_vaf-wrist_decoder_vaf;
    mean_vaf = mean(vaf_diff);
    correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
    % upp = tinv(0.975,99);
    low = tinv(0.01,num_folds*num_repeats-1);
    % vafcihi = mean_vaf + upp * sqrt(correction*vardiff);
    vafcilo = mean_vaf + low * sqrt(correction*var(vaf_diff));
    p_val = tcdf(-mean_vaf./sqrt(correction*var(vaf_diff)),num_folds*num_repeats-1);

    % plot vafs
    barwidth = 0.4;
    figure
    bar([1 3 2 4]*0.20,mean([hand_decoder_vaf neur_decoder_vaf]),barwidth,'facecolor',[0.5 0.5 0.5],'edgecolor','none')
    hold on
    plot([1 2]*0.20,[hand_decoder_vaf(:,1) neur_decoder_vaf(:,1)]','-k','linewidth',1)
    plot([3 4]*0.20,[hand_decoder_vaf(:,2) neur_decoder_vaf(:,2)]','-k','linewidth',1)
    plot([1 3 2 4]*0.20, [hand_decoder_vaf neur_decoder_vaf]','.k','markersize',30)

    % plot reverse model
    % figure
    bar([1 3 2 4]*0.20+1,mean([reverse_hand_decoder_vaf reverse_neur_decoder_vaf]),barwidth,'facecolor',[0.5 0.5 0.5],'edgecolor','none')
    hold on
    plot([1 2]*0.20+1,[reverse_hand_decoder_vaf(:,1) reverse_neur_decoder_vaf(:,1)]','-k','linewidth',1)
    plot([3 4]*0.20+1,[reverse_hand_decoder_vaf(:,2) reverse_neur_decoder_vaf(:,2)]','-k','linewidth',1)
    plot([1 3 2 4]*0.20+1, [reverse_hand_decoder_vaf reverse_neur_decoder_vaf]','.k','markersize',30)

    % plot control
    % figure
    bar([1 3 2 4]*0.20+2,mean([wrist_decoder_vaf wrist_neur_decoder_vaf]),barwidth,'facecolor',[0.5 0.5 0.5],'edgecolor','none')
    hold on
    plot([1 2]*0.20+2,[wrist_decoder_vaf(:,1) wrist_neur_decoder_vaf(:,1)]','-k','linewidth',1)
    plot([3 4]*0.20+2,[wrist_decoder_vaf(:,2) wrist_neur_decoder_vaf(:,2)]','-k','linewidth',1)
    plot([1 3 2 4]*0.20+2, [wrist_decoder_vaf wrist_neur_decoder_vaf]','.k','markersize',30)

    % plot niceness
    xtick = ones(3,1)*(1:4)*0.20 + (0:2)';
    xtick = reshape(xtick',12,1);
    xticklabel = {'Hand-only pos','Hand+Neuron pos','Hand vel','Hand+Neuron vel',...
        'Elbow-only pos','Elbow+Neuron pos','Elbow vel','Elbow+Neuron vel',...
        'HandControl-only pos','HandControl+Neuron pos','HandControl vel','HandControl+Neuron vel'};
    set(gca,'box','off','tickdir','out','xtick',xtick,...
        'xticklabel',xticklabel,...
        'xlim',[0 3],'ylim',[0 1])
    ylabel 'Fraction VAF'
    xlabel 'Model'
    title 'Decoding performance'

