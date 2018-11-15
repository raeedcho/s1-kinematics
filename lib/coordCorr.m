function [results] = coordCorr(td,params)
% This function calculates how much variance of muscle kinematics you can explain
% by linearly mapping them to other coordinate frames.

%% Set up models
    num_folds = 10; % 5 is default number of folds, no need to pass in
    num_repeats = 10; % 20 is default number of repeats, no need to pass in
    num_musc_pcs = 5;
    model_type = 'linmodel';
    model_aliases = {'ext','handelbow','musc'};
    assignParams(who,params);

    model_names = strcat(model_type,'_',model_aliases,'_model');
    num_models = length(model_names);
    
    % set up model parameters
    musc_len_idx = find(contains(td(1).opensim_names,'_len'));
    musc_vel_idx = find(contains(td(1).opensim_names,'_muscVel'));
    lenModel_params = cell(num_models,1);
    velModel_params = cell(num_models,1);
    for modelnum = 1:num_models
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
            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals',{{'opensim_len_pca',1:num_musc_pcs}},...
                                    'out_signals',{{'opensim',musc_len_idx}});
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
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals',{{'opensim_muscVel_pca',1:num_musc_pcs}},...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        case 'ext'
            markername = 'Marker_1';
            [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td(1).marker_names);
            assert(all(point_exists),'Hand marker does not exist?')
            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals',{{'markers',marker_hand_idx}},...
                                    'out_signals',{{'opensim',musc_len_idx}});
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals',{{'marker_vel',marker_hand_idx}},...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        case 'opensim_ext'
            opensim_handPos_idx = find(contains(td(1).opensim_names,'_handPos'));
            opensim_handVel_idx = find(contains(td(1).opensim_names,'_handVel'));
            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals',{{'opensim',opensim_handPos_idx}},...
                                    'out_signals',{{'opensim',musc_len_idx}});
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals',{{'opensim',opensim_handVel_idx}},...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        case 'ego'
            % add in spherical coordinates
            td = addCoordPoint2TD(td,struct('method','markers','coord','sph','point','hand'));
            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals','markers_sph_hand_pos',...
                                    'out_signals',{{'opensim',musc_len_idx}});
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals','markers_sph_hand_vel',...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        case 'opensim_ego'
            % add in spherical coordinates
            td = addCoordPoint2TD(td,struct('method','opensim','coord','sph','point','hand'));
            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals','opensim_sph_hand_pos',...
                                    'out_signals',{{'opensim',musc_len_idx}});
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals','opensim_sph_hand_vel',...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        case 'cyl'
            % add in cylindrical coordinates
            td = addCoordPoint2TD(td,struct('method','markers','coord','cyl','point','hand'));
            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals','markers_cyl_hand_pos',...
                                    'out_signals',{{'opensim',musc_len_idx}});
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals','markers_cyl_hand_vel',...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        case 'opensim_cyl'
            % add in cylindrical coordinates
            td = addCoordPoint2TD(td,struct('method','opensim','coord','cyl','point','hand'));
            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals','opensim_cyl_hand_pos',...
                                    'out_signals',{{'opensim',musc_len_idx}});
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals','opensim_cyl_hand_vel',...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        case 'joint'
            opensim_jointAng_idx = find(contains(td(1).opensim_names,'_ang'));
            opensim_jointVel_idx = find(contains(td(1).opensim_names,'_vel'));
            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals',{{'opensim',opensim_jointAng_idx}},...
                                    'out_signals',{{'opensim',musc_len_idx}});
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals',{{'opensim',opensim_jointVel_idx}},...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        case 'handelbow'
            % indices for cartesian hand coordinates
            markername = 'Marker_1';
            [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td(1).marker_names);
            assert(all(point_exists),'Hand marker does not exist?')

            markername = 'Pronation_Pt1';
            [point_exists,marker_elbow_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td(1).marker_names);
            assert(all(point_exists),'Elbow marker does not exist?')

            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals',{{'markers',[marker_hand_idx marker_elbow_idx]}},...
                                    'out_signals',{{'opensim',musc_len_idx}});
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals',{{'marker_vel',[marker_hand_idx marker_elbow_idx]}},...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        case 'markers_pca'
            % Get PCA for marker space
            td = dimReduce(td,struct('signals','markers'));
            td = dimReduce(td,struct('signals','marker_vel'));
            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals',{{'markers_pca',1:num_musc_pcs}},...
                                    'out_signals',{{'opensim',musc_len_idx}});
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals',{{'marker_vel_pca',1:num_musc_pcs}},...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        case 'opensim_handelbow'
            opensim_handPos_idx = find(contains(td(1).opensim_names,'_handPos'));
            opensim_handVel_idx = find(contains(td(1).opensim_names,'_handVel'));
            opensim_elbowPos_idx = find(contains(td(1).opensim_names,'_elbowPos'));
            opensim_elbowVel_idx = find(contains(td(1).opensim_names,'_elbowVel'));
            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals',{{'opensim',[opensim_handPos_idx opensim_elbowPos_idx]}},...
                                    'out_signals',{{'opensim',musc_len_idx}});
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals',{{'opensim',[opensim_handVel_idx opensim_elbowVel_idx]}},...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        case 'ego_handelbow'
            % transform into new coord system
            td = addCoordPoint2TD(td,struct('method','markers','coord','sph','point','hand'));
            td = addCoordPoint2TD(td,struct('method','markers','coord','sph','point','elbow'));
            lenModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_lenModel'],...
                                    'in_signals',{{'markers_sph_hand_pos';'markers_sph_elbow_pos'}},...
                                    'out_signals',{{'opensim',musc_len_idx}});
            velModel_params{modelnum} = struct('model_type',model_type,...
                                    'model_name',[model_aliases{modelnum} '_velModel'],...
                                    'in_signals',{{'markers_sph_hand_vel';'markers_sph_elbow_vel'}},...
                                    'out_signals',{{'opensim',musc_vel_idx}});
        otherwise
            error('Unrecognized model_alias')
        end
    end

%% cross-validation loop for model predictions
    td_test = cell(num_repeats,num_folds);
    lenPredVAF = zeros(num_repeats,num_folds,num_models);
    velPredVAF = zeros(num_repeats,num_folds,num_models);
    repeat_tic = tic;
    fprintf('Starting %dx%d cross-validation at %f s\n',num_repeats,num_folds,toc(repeat_tic))
    for repeatnum = 1:num_repeats
        indices = crossvalind('Kfold',length(td),num_folds);
        fold_tic = tic;
        fprintf('\tStarting %d fold cross-validation at %f s\n',num_folds,toc(fold_tic))
        for foldnum = 1:num_folds
            % split into testing and training
            test_idx = (indices==foldnum);
            train_idx = ~test_idx;
            td_train = td(train_idx);
            td_test{repeatnum,foldnum} = td(test_idx);

            % Fit models on training data
            % model_tic = tic;
            for modelnum = 1:num_models
                [lenPredVAF(repeatnum,foldnum,modelnum),td_test{repeatnum,foldnum}] = evalFold(td_train,td_test{repeatnum,foldnum},lenModel_params{modelnum});
                [velPredVAF(repeatnum,foldnum,modelnum),td_test{repeatnum,foldnum}] = evalFold(td_train,td_test{repeatnum,foldnum},velModel_params{modelnum});
                % fprintf('\t\tFinished model %d of %d at %f s\n',modelnum,num_models,toc(model_tic))
            end
            fprintf('\tFinished fold %d of %d at %f s\n',foldnum,num_folds,toc(fold_tic))
        end
        fprintf('Finished repeat %d of %d at %f s\n',repeatnum,num_repeats,toc(repeat_tic))
    end

%% package output
    results.params.num_folds = num_folds;
    results.params.num_repeats = num_repeats;
    results.params.num_musc_pcs = num_musc_pcs;
    results.params.model_type = model_type;
    results.params.model_aliases = model_aliases;
    results.model_names = model_names;
    results.lenModel_params = lenModel_params;
    results.velModel_params = velModel_params;
    results.td_test = td_test;
    results.lenPredVAF = lenPredVAF;
    results.velPredVAF = velPredVAF;

end

function [predVAF,td_test] = evalFold(td_train,td_test,model_params)
% evaluates predictive capability of model given model params
    [~,model_info] = getModel(td_train,model_params);
    
    % predict firing rates for td_test
    td_test = getModel(td_test,model_info);
    
    % evaluate the models
    Y = getSig(td_test,model_params.out_signals);
    Y_hat = getSig(td_test,strcat(model_params.model_type,'_',model_params.model_name));
    Y_bar = repmat(mean(Y),size(Y,1),1);

    SSE = sum(sum((Y-Y_hat).^2,1),2);
    SST = sum(sum((Y-Y_bar).^2,1),2);

    predVAF = 1 - SSE/SST;
end
