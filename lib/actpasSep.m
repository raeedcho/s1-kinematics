function results = actpasSep(td,params)

    %% Prep the data
        % split into active and passive
        [~,td_act] = getTDidx(td,'ctrHoldBump',false);
        % clean nans out...?
        nanners = isnan(cat(1,td_act.tgtDir));
        td_act = td_act(~nanners);
        [~,td_pas] = getTDidx(td,'ctrHoldBump',true);

        % find the movmement onsets
        td_act = getNorm(td_act,struct('signals','vel','norm_name','speed'));
        td_act = getMoveOnsetAndPeak(td_act,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
        td_pas = getNorm(td_pas,struct('signals','vel','norm_name','speed'));
        td_pas = getMoveOnsetAndPeak(td_pas,struct('start_idx','idx_bumpTime','end_idx','idx_goCueTime','method','peak','min_ds',1));

        % trim to just movements
        td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',14});
        td_pas = trimTD(td_pas,{'idx_movement_on',0},{'idx_movement_on',14});

        % find average over the movement
        td_act = binTD(td_act,'average');
        td_pas = binTD(td_pas,'average');

        % even out sizes
        minsize = min(length(td_act),length(td_pas));
        td_act = td_act(1:minsize);
        td_pas = td_pas(1:minsize);
        td = cat(2,td_act,td_pas);

    %% set up model variables
        num_folds = 5; % 5 is default number of folds, no need to pass in
        num_repeats = 20; % 20 is default number of repeats, no need to pass in
        num_musc_pcs = 5;
        model_type = 'glm';
        model_aliases = {'ext','extforce','joint','musc','handelbow'};
        arrayname = 'S1';
        assignParams(who,params);
        neural_signals = [arrayname '_FR'];

        model_names = [strcat(model_type,'_',model_aliases,'_model') {neural_signals}];
        num_models = length(model_names);

        % set up glm parameters
        glm_params = cell(num_models-1,1);
        for modelnum = 1:num_models-1
            switch model_aliases{modelnum}
            case 'musc'
                % Do PCA on muscle space
                % do PCA on muscles, training on only the training set
                % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
                % PCAparams = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_len') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
                %                     'do_plot',true);
                PCAparams = struct('signals','muscle_len', 'do_plot',false);
                [td,~] = dimReduce(td,PCAparams);
                % get velocity PCA
                % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
                % PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
                %                     'do_plot',true);
                PCAparams_vel = struct('signals','muscle_vel', 'do_plot',false);
                [td,~] = dimReduce(td,PCAparams_vel);
                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'muscle_len_pca',1:num_musc_pcs;'muscle_vel_pca',1:num_musc_pcs}},...
                                        'out_signals',neural_signals);
            case 'ext'
                markername = 'Marker_1';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td(1).marker_names);
                assert(all(point_exists),'Hand marker does not exist?')
                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'markers',marker_hand_idx;'marker_vel',marker_hand_idx}},...
                    'out_signals',neural_signals);
            case 'handle_ext'
                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'pos';'vel'}},...
                    'out_signals',neural_signals);
            case 'opensim_ext'
                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals','opensim_hand_pos',...
                                        'out_signals',neural_signals);
            case 'extforce'
                markername = 'Marker_1';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td(1).marker_names);
                assert(all(point_exists),'Hand marker does not exist?')
                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'markers',marker_hand_idx;'marker_vel',marker_hand_idx;'force',1:3}},...
                    'out_signals',neural_signals);
            case 'handle_extforce'
                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'pos',1:2;'vel',1:2;'force',1:2}},...
                    'out_signals',neural_signals);
            case 'ego'
                % add in spherical coordinates
                td = addCoordPoint2TD(td,struct('method','markers','coord','sph','point','hand'));
                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'markers_sph_hand_pos';'markers_sph_hand_vel'}},...
                                        'out_signals',neural_signals);
            case 'opensim_ego'
                % add in spherical coordinates
                error('Opensim_ego model not compatible with new trial_data structure yet')
                % td = addCoordPoint2TD(td,struct('method','opensim','coord','sph','point','hand'));
                % glm_params{modelnum} = struct('model_type',model_type,...
                %                         'model_name',[model_aliases{modelnum} '_model'],...
                %                         'in_signals',{{'opensim_sph_hand_pos';'opensim_sph_hand_vel'}},...
                %                         'out_signals',neural_signals);
            case 'cyl'
                % add in cylindrical coordinates
                td = addCoordPoint2TD(td,struct('method','markers','coord','cyl','point','hand'));
                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'markers_cyl_hand_pos';'markers_cyl_hand_vel'}},...
                                        'out_signals',neural_signals);
            case 'opensim_cyl'
                % add in cylindrical coordinates
                error('Opensim_cyl model not compatible with new trial_data structure yet')
                % td = addCoordPoint2TD(td,struct('method','opensim','coord','cyl','point','hand'));
                % glm_params{modelnum} = struct('model_type',model_type,...
                %                         'model_name',[model_aliases{modelnum} '_model'],...
                %                         'in_signals',{{'opensim_cyl_hand_pos';'opensim_cyl_hand_vel'}},...
                %                         'out_signals',neural_signals);
            case 'joint'
                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'joint_ang';'joint_vel'}},...
                                        'out_signals',neural_signals);
            case 'handelbow'
                % indices for cartesian hand coordinates
                markername = 'Marker_1';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td(1).marker_names);
                assert(all(point_exists),'Hand marker does not exist?')

                markername = 'Pronation_Pt1';
                [point_exists,marker_elbow_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td(1).marker_names);
                assert(all(point_exists),'Elbow marker does not exist?')

                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'markers',[marker_hand_idx marker_elbow_idx];'marker_vel',[marker_hand_idx marker_elbow_idx]}},...
                                        'out_signals',neural_signals);
            case 'markers_pca'
                % Get PCA for marker space
                td = dimReduce(td,struct('signals','markers'));
                td = dimReduce(td,struct('signals','marker_vel'));
                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'markers_pca',1:num_musc_pcs;'marker_vel_pca',1:num_musc_pcs}},...
                                        'out_signals',neural_signals);
            case 'opensim_handelbow'
                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'opensim_hand_pos';'opensim_hand_vel';'opensim_elbow_pos';'opensim_elbow_vel'}},...
                                        'out_signals',neural_signals);
            case 'ego_handelbow'
                % transform into new coord system
                td = addCoordPoint2TD(td,struct('method','markers','coord','sph','point','hand'));
                td = addCoordPoint2TD(td,struct('method','markers','coord','sph','point','elbow'));
                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'markers_sph_hand_pos';'markers_sph_hand_vel';'markers_sph_elbow_pos';'markers_sph_elbow_vel'}},...
                                        'out_signals',neural_signals);
            otherwise
                error('Unrecognized model_alias')
            end
        end

    %% cross-validate the separabilities
        [~,td_act] = getTDidx(td,'ctrHoldBump',false);
        [~,td_pas] = getTDidx(td,'ctrHoldBump',true);

        repeat_tic = tic;
        fprintf('Starting %dx%d cross-validation at time: %f\n',num_repeats,num_folds,toc(repeat_tic))
        meta_table = cell2table({td(1).monkey,td(1).date_time,td(1).task},...
            'VariableNames',{'monkey','date_time','task'});
        meta_table.Properties.VariableDescriptions = repmat({'meta'},1,3);
        [trial_table_cell,lda_table_cell] = deal(cell(num_repeats,num_folds));
        for repeatnum = 1:num_repeats
            foldidx = crossvalind('kfold',minsize,num_folds);
            fold_tic = tic;
            for foldnum = 1:num_folds
                % get crossval table entry
                crossval_table = table(uint16([repeatnum foldnum]),'VariableNames',{'crossvalID'});
                crossval_table.Properties.VariableDescriptions = {'meta'};

                % split into training and testing
                train_idx = (foldidx~=foldnum);
                td_train = cat(2,td_act(train_idx),td_pas(train_idx));
                test_idx = (foldidx==foldnum);
                td_test = cat(2,td_act(test_idx),td_pas(test_idx));

                % train and test models
                glm_info = cell(1,length(model_names)-1);
                [lda_coeff,seps,model_fr] = deal(cell(1,length(model_names)));
                train_class = cat(1,td_train.ctrHoldBump);
                test_class = cat(1,td_test.ctrHoldBump);
                for modelnum = 1:length(model_names)
                    if modelnum~=length(model_names)
                        [td_train,glm_info{modelnum}] = getModel(td_train,glm_params{modelnum});

                        % predict firing rates
                        td_test = getModel(td_test,glm_info{modelnum});
                    end
        
                    % get LDA models
                    train_fr = cat(1,td_train.(model_names{modelnum}));
                    lda_mdl = fitcdiscr(train_fr,train_class);
                    lda_coeff{modelnum} = [lda_mdl.Coeffs(1,2).Const;lda_mdl.Coeffs(1,2).Linear]';

                    % get separabilities from LDA models
                    model_fr{modelnum} = cat(1,td_test.(model_names{modelnum}));
                    seps{modelnum} = sum(predict(lda_mdl,model_fr{modelnum}) == test_class)/length(test_class);
                end

                % compile trial table
                % get direction of trials for trial table
                bump_dir = cat(1,td_test.bumpDir);
                tgt_dir = cat(1,td_test.tgtDir);
                trial_dir = zeros(size(test_class));
                trial_dir(test_class) = bump_dir(test_class)*pi/180;
                trial_dir(~test_class) = tgt_dir(~test_class)*pi/180;
                trial_dir_table = table(trial_dir,'VariableNames',{'trialDir'});
                trial_dir_table.Properties.VariableDescriptions = {'circular'};
                % get class for trial table
                test_class_table = table(test_class,'VariableNames',{'isPassive'});
                test_class_table.Properties.VariableDescriptions = {'meta'};
                % get model fr for trial table
                model_fr_table = table(model_fr{:},'VariableNames',[strcat(model_aliases,'_predFR') {neural_signals}]);
                model_fr_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_names));
                trial_table_cell{repeatnum,foldnum} = horzcat(...
                    repmat(meta_table,length(trial_dir),1),...
                    repmat(crossval_table,length(trial_dir),1),...
                    test_class_table,...
                    trial_dir_table,...
                    model_fr_table);

                % construct LDA table entry
                lda_coeff_table = cell2table(lda_coeff,'VariableNames',strcat([model_aliases {neural_signals}],'_lda_coeff'));
                lda_coeff_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_names));
                sep_table = cell2table(seps,'VariableNames',strcat([model_aliases {neural_signals}],'_sep'));
                sep_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_names));
                lda_table_cell{repeatnum,foldnum} = horzcat(...
                    meta_table,...
                    crossval_table,...
                    lda_coeff_table,...
                    sep_table);
        
                fprintf('\tEvaluated fold %d at time: %f\n',foldnum,toc(fold_tic))
            end
            fprintf('Evaluated repeat %d at time: %f\n',repeatnum,toc(repeat_tic))
        end

    %% Package results
        % get one separability table
        lda_table = vertcat(lda_table_cell{:});
        trial_table = vertcat(trial_table_cell{:});

        results = struct('lda_table',lda_table,'trial_table',trial_table);
