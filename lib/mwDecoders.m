function results = mwDecoders(td,params)
%% decoder analysis - predict hand movement and arm movement from spikes
    %% Setup
    % to predict elbow from both neurons and hand
    elbow_idx = 28:30;
    hand_idx = 4:6;
    neur_name = 'S1_spikes';
    add_shift = true;
    num_repeats = 10;
    num_folds = 10;
    if nargin > 1
        assignParams(who,params)
    end

    % preallocate vaf holders
    neur_decoder_vaf = zeros(num_repeats,num_folds,2);
    hand_decoder_vaf = zeros(num_repeats,num_folds,2);

    % set model parameters
    % decoder parameters
    if add_shift
        neural_signals = {neur_name,'all';sprintf('%s_shift',neur_name),'all'};
    else
        neural_signals = {neur_name,'all'};
    end
    elbow_signals = {'markers',elbow_idx;'marker_vel',elbow_idx};
    hand_signals = {'markers',hand_idx;'marker_vel',hand_idx};

    % predict elbow from either just hand or hand and neurons
    neur_decoder_params = struct('model_type','linmodel','model_name','neur_decoder',...
        'in_signals',{[neural_signals;hand_signals]},...
        'out_signals',{elbow_signals});
    hand_decoder_params = struct('model_type','linmodel','model_name','hand_decoder',...
        'in_signals',{hand_signals},...
        'out_signals',{elbow_signals});
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

            % evaluate models and save into array
            td_test = td(inds==foldnum);
            td_test = getModel(td_test,neur_decoder);
            td_test = getModel(td_test,hand_decoder);

            % get error on elbow
            elbow_true = get_vars(td_test,{'markers',elbow_idx;'marker_vel',elbow_idx});
            elbow_mean = mean(elbow_true);
            elbow_pred_neur = get_vars(td_test,check_signals(td_test,'linmodel_neur_decoder'));
            elbow_pred_hand = get_vars(td_test,check_signals(td_test,'linmodel_hand_decoder'));

            % calculate fraction of variance explained
            SS_total = square_sum(elbow_true-elbow_mean);
            neur_decoder_vaf(repeatnum,foldnum,:) = 1 - square_sum(elbow_pred_neur-elbow_true)./SS_total;
            hand_decoder_vaf(repeatnum,foldnum,:) = 1 - square_sum(elbow_pred_hand-elbow_true)./SS_total;

            % display
            fprintf('\tEvaluated fold number %d at time %f\n',foldnum,toc(fold_tic))
        end
    end
    
    neur_decoder_vaf = reshape(neur_decoder_vaf,num_repeats*num_folds,2);
    hand_decoder_vaf = reshape(hand_decoder_vaf,num_repeats*num_folds,2);

    %% Package results
    results.neur_decoder_vaf = neur_decoder_vaf;
    results.hand_decoder_vaf = hand_decoder_vaf;
