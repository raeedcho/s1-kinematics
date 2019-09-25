function results = sepTime(td_bin,params)

    %% set up model variables
        num_folds = 5; % 5 is default number of folds, no need to pass in
        num_repeats = 20; % 20 is default number of repeats, no need to pass in
        crossval_lookup = []; % table with columns for crossvalID and trialID (to replicate a crossval split from a previous run)
        num_musc_pcs = 5;
        num_pcs = 5; % number of PCs to train LDA on
        model_type = 'glm';
        model_aliases = {'ext','extforce','handelbow'};
        neural_signals = 'S1_FR';
        which_units = 'all'; % replace with a list of indices for which units to use for separability
        timeshifts = [];
        assignParams(who,params);

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
                [td_bin,~] = dimReduce(td_bin,PCAparams);
                % get velocity PCA
                % need to drop a muscle: for some reason, PCA says rank of muscle kinematics matrix is 38, not 39.
                % PCAparams_vel = struct('signals',{{'opensim',find(contains(td(1).opensim_names,'_muscVel') & ~contains(td(1).opensim_names,'tricep_lat'))}},...
                %                     'do_plot',true);
                PCAparams_vel = struct('signals','muscle_vel', 'do_plot',false);
                [td_bin,~] = dimReduce(td_bin,PCAparams_vel);
                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'muscle_len_pca',1:num_musc_pcs;'muscle_vel_pca',1:num_musc_pcs}},...
                                        'out_signals',neural_signals);
            case 'ext'
                markername = 'Marker_3';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_bin(1).marker_names);
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
                markername = 'Marker_3';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_bin(1).marker_names);
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
                td_bin = addCoordPoint2TD(td_bin,struct('method','markers','coord','sph','point','hand'));
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
                td_bin = addCoordPoint2TD(td_bin,struct('method','markers','coord','cyl','point','hand'));
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
                markername = 'Marker_3';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_bin(1).marker_names);
                assert(all(point_exists),'Hand marker does not exist?')

                markername = 'Pronation_Pt1';
                [point_exists,marker_elbow_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_bin(1).marker_names);
                assert(all(point_exists),'Elbow marker does not exist?')

                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'markers',[marker_hand_idx marker_elbow_idx];'marker_vel',[marker_hand_idx marker_elbow_idx]}},...
                                        'out_signals',neural_signals);
            case 'markers_pca'
                % Get PCA for marker space
                td_bin = dimReduce(td_bin,struct('signals','markers'));
                td_bin = dimReduce(td_bin,struct('signals','marker_vel'));
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
                td_bin = addCoordPoint2TD(td_bin,struct('method','markers','coord','sph','point','hand'));
                td_bin = addCoordPoint2TD(td_bin,struct('method','markers','coord','sph','point','elbow'));
                glm_params{modelnum} = struct('model_type',model_type,...
                                        'model_name',[model_aliases{modelnum} '_model'],...
                                        'in_signals',{{'markers_sph_hand_pos';'markers_sph_hand_vel';'markers_sph_elbow_pos';'markers_sph_elbow_vel'}},...
                                        'out_signals',neural_signals);
            case 'ext_actpasbaseline' % Note: only works when each trial has only one bin (avg firing rate per trial)
                % Use hand position and whether or not the trial was active or passive
                % effectively allows active and passive trials to have a different baseline FR in the model
                warning('ext_actpasbaseline model will currently only work when each trial has one bin (average firing rate for trial)')
                markername = 'Marker_3';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_bin(1).marker_names);
                assert(all(point_exists),'Hand marker does not exist?')
                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'markers',marker_hand_idx;'marker_vel',marker_hand_idx;'ctrHoldBump',1}},...
                    'out_signals',neural_signals);
            otherwise
                error('Unrecognized model_alias')
            end
        end

    %% cross-validate the separabilities
        [~,td_act] = getTDidx(td_bin,'ctrHoldBump',false);
        [~,td_pas] = getTDidx(td_bin,'ctrHoldBump',true);

        repeat_tic = tic;
        fprintf('Starting %dx%d cross-validation at time: %f\n',num_repeats,num_folds,toc(repeat_tic))
        meta_table = cell2table({td_bin(1).monkey,td_bin(1).date_time,td_bin(1).task},...
            'VariableNames',{'monkey','date_time','task'});
        meta_table.Properties.VariableDescriptions = repmat({'meta'},1,3);
        [trial_table_cell,lda_table_cell] = deal(cell(num_repeats,num_folds));
        for repeatnum = 1:num_repeats
            foldidx = crossvalind('kfold',length(td_act),num_folds);
            fold_tic = tic;
            for foldnum = 1:num_folds
                % get crossval table entry
                crossval_table = table(uint16([repeatnum foldnum]),'VariableNames',{'crossvalID'});
                crossval_table.Properties.VariableDescriptions = {'meta'};

                % split into training and testing
                if isempty(crossval_lookup)
                    train_idx = (foldidx~=foldnum);
                    td_train = cat(2,td_act(train_idx),td_pas(train_idx));
                    test_idx = (foldidx==foldnum);
                    td_test = cat(2,td_act(test_idx),td_pas(test_idx));
                else
                    % read the crossval_lookup
                    [~,current_crossval] = getNTidx(crossval_lookup,'crossvalID',[repeatnum foldnum]);
                    td_whole = cat(2,td_act,td_pas);
                    whole_ids = cat(1,td_whole.trialID);
                    whole_isPassive = cat(1,false(length(td_act),1),true(length(td_pas),1));
                    test_idx_whole = ismember([whole_ids whole_isPassive],current_crossval{:,{'trialID' 'isPassive'}},'rows');
                    train_idx_whole = ~test_idx_whole;
                    
                    td_test = td_whole(test_idx_whole);
                    td_train = td_whole(train_idx_whole);
                end

                % set up meta for trial table
                trial_id = table(cat(1,td_test.trialID),'VariableNames',{'trialID'});
                trial_id.Properties.VariableDescriptions = {'meta'};
                bump_dir = cat(1,td_test.bumpDir);
                tgt_dir = cat(1,td_test.tgtDir);
                % get class for trial table
                train_class = cat(1,td_train.ctrHoldBump);
                test_class = cat(1,td_test.ctrHoldBump);
                test_class_table = table(test_class,'VariableNames',{'isPassive'});
                test_class_table.Properties.VariableDescriptions = {'meta'};
                % get direction of trials for trial table
                trial_dir = zeros(size(test_class));
                trial_dir(test_class) = bump_dir(test_class)*pi/180;
                trial_dir(~test_class) = tgt_dir(~test_class)*pi/180;
                trial_dir_table = table(trial_dir,'VariableNames',{'trialDir'});
                trial_dir_table.Properties.VariableDescriptions = {'circular'};

                % train and test models
                [lda_coeff,pca_coeff,pca_mu] = deal(cell(1,length(timeshifts)));
                % [neural_seps,neural_margin] = deal(zeros(1,length(timeshifts)));
                neural_seps = zeros(1,length(timeshifts));
                neural_margin = zeros(length(td_test),length(timeshifts));

                [input_lda_coeff] = deal(cell(1,length(model_aliases)));
                [input_seps,input_margin] = deal(cell(1,length(model_aliases)));
                for modelnum = 1:length(model_aliases)
                    input_lda_coeff{modelnum} = cell(length(timeshifts),1);
                    input_seps{modelnum} = zeros(1,length(timeshifts));
                end

                for timeshiftnum = 1:length(timeshifts)
                    % trim TD for this shift
                    shift = timeshifts(timeshiftnum);
                    td_train_shift = trimTD(td_train,{'idx_movement_on',shift},{'idx_movement_on',shift-1+floor(0.05/td_train(1).bin_size)});
                    td_train_shift = binTD(td_train_shift,'average');
                    td_test_shift = trimTD(td_test,{'idx_movement_on',shift},{'idx_movement_on',shift-1+floor(0.05/td_test(1).bin_size)});
                    td_test_shift = binTD(td_test_shift,'average');

                    for modelnum = 1:length(model_aliases)
                        % get input separability
                        if strcmpi(model_aliases{modelnum},'ext_actpasbaseline')
                            % actpas input has zero variance and perfectly separates the classes
                            % So just put down the coefficients known to separate
                            input_signals = getSig(td_train_shift,glm_params{modelnum}.in_signals);

                            % set offset coeff to -0.5
                            input_lda_coeff{modelnum}{timeshiftnum} = zeros(1,size(input_signals,2)+1);
                            input_lda_coeff{modelnum}{timeshiftnum}(1) = -0.5;

                            % figure out which column is ctrHoldBump
                            colnums = cumsum(cellfun(@length,glm_params{modelnum}.in_signals(:,2)));
                            actpas_col = colnums(strcmpi(glm_params{modelnum}.in_signals(:,1),'ctrHoldBump'));
                            input_lda_coeff{modelnum}{timeshiftnum}(actpas_col+1) = 1;

                            % actpas input should have 100% separability
                            input_seps{modelnum}(timeshiftnum) = 1;
                         else
                            input_signals = getSig(td_train_shift,glm_params{modelnum}.in_signals);
                            input_lda_mdl = fitcdiscr(input_signals,train_class);
                            input_lda_coeff{modelnum}{timeshiftnum} = [input_lda_mdl.Coeffs(2,1).Const;input_lda_mdl.Coeffs(2,1).Linear]';

                            input_signals = getSig(td_test_shift,glm_params{modelnum}.in_signals);
                            input_seps{modelnum}(:,timeshiftnum) = sum(predict(input_lda_mdl,input_signals) == test_class)/length(test_class);
                        end

                        % calculate input margins
                        input_signals = getSig(td_test_shift,glm_params{modelnum}.in_signals);
                        input_margin{modelnum}(:,timeshiftnum) = (input_signals*input_lda_coeff{modelnum}{timeshiftnum}(2:end)' + input_lda_coeff{modelnum}{timeshiftnum}(1))/norm(input_lda_coeff{modelnum}{timeshiftnum}(2:end));
                    end

                    % get PCA for neural signals
                    [~,pca_info] = dimReduce(td_train_shift,struct('signals',{{neural_signals,which_units}}));
                    td_test_shift = dimReduce(td_test_shift,pca_info);
                    pca_coeff{timeshiftnum} = pca_info.w(:,1:num_pcs);
                    pca_mu{timeshiftnum} = pca_info.mu;

                    % get LDA models
                    train_fr = getSig(td_train_shift,{neural_signals,which_units});
                    lda_mdl = fitcdiscr((train_fr-pca_mu{timeshiftnum})*pca_coeff{timeshiftnum},train_class);
                    lda_coeff{timeshiftnum} = [lda_mdl.Coeffs(2,1).Const;lda_mdl.Coeffs(2,1).Linear]';

                    % get separabilities from LDA models
                    test_fr = getSig(td_test_shift,{neural_signals,which_units});
                    neural_seps(timeshiftnum) = sum(predict(lda_mdl,(test_fr-pca_mu{timeshiftnum})*pca_coeff{timeshiftnum}) == test_class)/length(test_class);

                    % calculate margin from lda boundary
                    neural_margin(:,timeshiftnum) = ((test_fr-pca_mu{timeshiftnum})*pca_coeff{timeshiftnum}*lda_coeff{timeshiftnum}(2:end)' + lda_coeff{timeshiftnum}(1))/norm(lda_coeff{timeshiftnum}(2:end));
                end
                        
                % get margins for trial table
                input_margin_table = table(input_margin{:},'VariableNames',strcat(model_aliases,'_input_margin'));
                input_margin_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_aliases));
                neural_margin_table = table(neural_margin,'VariableNames',strcat({neural_signals},'_self_margin'));
                neural_margin_table.Properties.VariableDescriptions = {'linear'};

                % assemble trial table
                trial_table_cell{repeatnum,foldnum} = horzcat(...
                    repmat(meta_table,length(trial_dir),1),...
                    repmat(crossval_table,length(trial_dir),1),...
                    trial_id,...
                    test_class_table,...
                    trial_dir_table,...
                    input_margin_table,...
                    neural_margin_table);

                % construct LDA table entry
                input_lda_coeff_table = cell2table(input_lda_coeff,'VariableNames',strcat(model_aliases,'_input_lda_coeff'));
                input_lda_coeff_table.Properties.VariableDescriptions = repmat({'cell'},1,length(model_aliases));
                lda_coeff_table = table(lda_coeff,'VariableNames',strcat({neural_signals},'_lda_coeff'));
                lda_coeff_table.Properties.VariableDescriptions = {'cell'};
                pca_coeff_table = table(pca_coeff,'VariableNames',strcat({neural_signals},'_pca_coeff'));
                pca_coeff_table.Properties.VariableDescriptions = {'cell'};
                pca_mu_table = table(pca_mu,'VariableNames',strcat({neural_signals},'_pca_mu'));
                pca_mu_table.Properties.VariableDescriptions = {'cell'};
                input_sep_table = cell2table(input_seps,'VariableNames',strcat(model_aliases,'_input_sep'));
                input_sep_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_aliases));
                neural_sep_table = table(neural_seps,'VariableNames',strcat({neural_signals},'_self_sep'));
                neural_sep_table.Properties.VariableDescriptions = {'linear'};
                lda_table_cell{repeatnum,foldnum} = horzcat(...
                    meta_table,...
                    crossval_table,...
                    input_lda_coeff_table,...
                    lda_coeff_table,...
                    pca_coeff_table,...
                    pca_mu_table,...
                    input_sep_table,...
                    neural_sep_table);
        
                fprintf('\tEvaluated fold %d at time: %f\n',foldnum,toc(fold_tic))
            end
            fprintf('Evaluated repeat %d at time: %f\n',repeatnum,toc(repeat_tic))
        end

    %% Package results
        % get one separability table
        lda_table = vertcat(lda_table_cell{:});
        trial_table = vertcat(trial_table_cell{:});

        results = struct('lda_table',lda_table,'trial_table',trial_table);
