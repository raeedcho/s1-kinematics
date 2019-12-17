function results = actpasSep(td_bin,params)

    %% set up model variables
        num_folds = 5; % 5 is default number of folds, no need to pass in
        num_repeats = 20; % 20 is default number of repeats, no need to pass in
        crossval_lookup = []; % table with columns for crossvalID and trialID (to replicate a crossval split from a previous run)
        num_musc_pcs = 5;
        num_pcs = 5; % number of PCs to train LDA on
        model_type = 'glm';
        model_aliases = {'ext','handelbow'};
        neural_signals = 'S1_FR';
        which_units = 'all'; % replace with a list of indices for which units to use for separability
        unit_guide = td_bin.S1_unit_guide;
        get_margins = false;
        do_population_stuff = false;
        assignParams(who,params);

        if ~strcmpi(which_units,'all')
            warning('which_units applies to only LDA calculations; code will calculate individual separabilities for all neurons')
        end

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
                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
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
                    'in_signals',{{'markers',marker_hand_idx;'marker_vel',marker_hand_idx;'force',1:6}},...
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
                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
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
            case 'handelbow_actpasbaseline' % Note: only works when each trial has only one bin (avg firing rate per trial)
                % Use hand position and whether or not the trial was active or passive
                % effectively allows active and passive trials to have a different baseline FR in the model
                warning('handelbow_actpasbaseline model will currently only work when each trial has one bin (average firing rate for trial)')
                markername = 'Marker_3';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_bin(1).marker_names);
                assert(all(point_exists),'Hand marker does not exist?')
                markername = 'Pronation_Pt1';
                [point_exists,marker_elbow_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_bin(1).marker_names);
                assert(all(point_exists),'Elbow marker does not exist?')
                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'markers',[marker_hand_idx marker_elbow_idx];'marker_vel',[marker_hand_idx marker_elbow_idx];'ctrHoldBump',1}},...
                    'out_signals',neural_signals);
            case 'ext_actpasbaseline'
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
            case 'handelbow_surprise'
                % get all the hand/elbow stuff ready
                markername = 'Marker_3';
                [point_exists,marker_hand_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_bin(1).marker_names);
                assert(all(point_exists),'Hand marker does not exist?')
                markername = 'Pronation_Pt1';
                [point_exists,marker_elbow_idx] = ismember(strcat(markername,'_',{'x','y','z'}),td_bin(1).marker_names);
                assert(all(point_exists),'Elbow marker does not exist?')
                
                % set up the td for a one-hot encoding of target direction only in the passive case
                target_dirs = unique(horzcat(td_bin.bumpDir));
                for trialnum = 1:length(td_bin)
                    td_bin(trialnum).target_surprise = ismember(target_dirs,td_bin(trialnum).bumpDir) & td_bin(trialnum).ctrHoldBump;
                end
                glm_params{modelnum} = struct(...
                    'model_type',model_type,...
                    'model_name',[model_aliases{modelnum} '_model'],...
                    'in_signals',{{'markers',[marker_hand_idx marker_elbow_idx];'marker_vel',[marker_hand_idx marker_elbow_idx];'target_surprise',1:length(target_dirs)}},...
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
        [trial_table_cell,lda_table_cell,foldEval] = deal(cell(num_repeats,num_folds));
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
                [lda_coeff,self_seps,true_seps,model_fr,lda_mdl,pca_coeff,pca_mu,self_margin,true_margin] = deal(cell(1,length(model_names)));
                [input_lda_coeff,input_seps,input_lda_mdl,input_margin] = deal(cell(1,length(model_aliases)));
                model_eval = cell(1,length(model_aliases));
                indiv_seps = zeros(length(unit_guide),length(model_names));
                for modelnum = 1:length(model_names)
                    if modelnum~=length(model_names)
                        % get eval metrics for neurons trained and tested on same condition
                        cond = {'act','pas'};
                        train_cond_eval = cell(1,length(cond));
                        for condnum = 1:2
                            [~,td_train_cond] = getTDidx(td_train,'ctrHoldBump',strcmpi(cond{condnum},'pas'));
                            [~,td_test_cond] = getTDidx(td_test,'ctrHoldBump',strcmpi(cond{condnum},'pas'));
                            
                            [~,glm_info] = getModel(td_train_cond,glm_params{modelnum});
                            td_test_cond = getModel(td_test_cond,glm_info);
                            
                            eval_params = glm_info;
                            eval_params.eval_metric = 'pr2';
                            eval_params.num_boots = 1;
                            eval_params.block_trials = true;
                            eval_params.trial_idx = 1:length(td_test_cond);
                            train_cond_eval{condnum} = evalModel(td_test_cond,eval_params);
                        end
                        model_eval_train_cond = table(train_cond_eval{:},...
                            'VariableNames',strcat(model_names{modelnum},'_train_',cond,'_eval'));
                        model_eval_train_cond.Properties.VariableDescriptions = repmat({'linear'},1,width(model_eval_train_cond));
                        
                        % train on half of full data...
                        [~,td_train_half] = getTDidx(td_train,'rand',floor(length(td_train)/2));
                        [~,td_test_half] = getTDidx(td_test,'rand',floor(length(td_test)/2));
                        [~,glm_info] = getModel(td_train_half,glm_params{modelnum});
                        td_test_half = getModel(td_test_half,glm_info);
                        eval_params = glm_info;
                        eval_params.eval_metric = 'pr2';
                        eval_params.num_boots = 1;
                        eval_params.block_trials = true;
                        act_eval_params = eval_params;
                        act_eval_params.trial_idx = getTDidx(td_test_half,'ctrHoldBump',false);
                        pas_eval_params = eval_params;
                        pas_eval_params.trial_idx = getTDidx(td_test_half,'ctrHoldBump',true);
                        model_eval_half_full_train = table(...
                            evalModel(td_test_half,act_eval_params),...
                            evalModel(td_test_half,pas_eval_params),...
                            'VariableNames',strcat(model_names{modelnum},{'_half_full_train_act_eval','_half_full_train_pas_eval'}));
                        model_eval_half_full_train.Properties.VariableDescriptions = repmat({'linear'},1,2);
                        
                        % Train and test over both conditions
                        [td_train,glm_info] = getModel(td_train,glm_params{modelnum});

                        % predict firing rates
                        td_test = getModel(td_test,glm_info);

                        % get individual neural evaluation metrics
                        eval_params = glm_info;
                        eval_params.eval_metric = 'pr2';
                        eval_params.num_boots = 1;
                        eval_params.block_trials = true;
                        eval_params.trial_idx = 1:length(td_test);
                        act_eval_params = eval_params;
                        act_eval_params.trial_idx = getTDidx(td_test,'ctrHoldBump',false);
                        pas_eval_params = eval_params;
                        pas_eval_params.trial_idx = getTDidx(td_test,'ctrHoldBump',true);
                        model_eval{modelnum} = table(...
                            evalModel(td_test,eval_params),...
                            evalModel(td_test,act_eval_params),...
                            evalModel(td_test,pas_eval_params),...
                            'VariableNames',strcat(model_names{modelnum},{'_eval','_act_eval','_pas_eval'}));
                        model_eval{modelnum}.Properties.VariableDescriptions = repmat({'linear'},1,3);
                        model_eval{modelnum} = horzcat(model_eval{modelnum},model_eval_train_cond,model_eval_half_full_train);

                        % get input LDA models
                        if endsWith(model_aliases{modelnum},'_actpasbaseline') || endsWith(model_aliases{modelnum},'_surprise')
                            % actpas input has zero variance and perfectly separates the classes
                            % So just put down the coefficients known to separate
                            input_signals = getSig(td_train,glm_params{modelnum}.in_signals);

                            % set offset coeff to -0.5
                            input_lda_coeff{modelnum} = zeros(1,size(input_signals,2)+1);
                            input_lda_coeff{modelnum}(1) = -0.5;

                            % figure out which column is ctrHoldBump
                            colnums = cumsum(cellfun(@length,glm_params{modelnum}.in_signals(:,2)));
                            actpas_col = colnums(strcmpi(glm_params{modelnum}.in_signals(:,1),'ctrHoldBump'));
                            input_lda_coeff{modelnum}(actpas_col+1) = 1;
                        else
                            input_signals = getSig(td_train,glm_params{modelnum}.in_signals);
                            input_lda_mdl{modelnum} = fitcdiscr(input_signals,train_class);
                            input_lda_coeff{modelnum} = [input_lda_mdl{modelnum}.Coeffs(2,1).Const;input_lda_mdl{modelnum}.Coeffs(2,1).Linear]';
                        end
                    end
                    
                    if do_population_stuff
                        % try sqrt transform (doesn't mess with model fitting because neural signals are last
                        % td_train = sqrtTransform(td_train,struct('signals',model_names{modelnum}));
                        % td_test = sqrtTransform(td_test,struct('signals',model_names{modelnum}));

                        % get PCA
                        [~,pca_info] = dimReduce(td_train,struct('signals',{{model_names{modelnum},which_units}}));
                        td_test = dimReduce(td_test,pca_info);
                        pca_coeff{modelnum} = pca_info.w(:,1:num_pcs);
                        pca_mu{modelnum} = pca_info.mu;

                        % get LDA models
                        % train_fr = cat(1,td_train.(model_names{modelnum}));
                        train_fr = getSig(td_train,{model_names{modelnum},which_units});
                        % lda_mdl{modelnum} = fitcdiscr(train_fr,train_class);
                        lda_mdl{modelnum} = fitcdiscr((train_fr-pca_mu{modelnum})*pca_coeff{modelnum},train_class);
                        lda_coeff{modelnum} = [lda_mdl{modelnum}.Coeffs(2,1).Const;lda_mdl{modelnum}.Coeffs(2,1).Linear]';
                    end

                    % get individual neuron separabilities
                    if modelnum == length(model_names)
                        train_fr = getSig(td_train,model_names{modelnum});
                        test_fr = getSig(td_test,model_names{modelnum});
                        for neuronnum = 1:size(train_fr,2)
                            indiv_lda_mdl = fitcdiscr(train_fr(:,neuronnum),train_class);
                            indiv_seps(neuronnum,modelnum) = sum(predict(indiv_lda_mdl,test_fr(:,neuronnum)) == test_class)/length(test_class);
                        end
                    end
                end

                
                for modelnum = 1:length(model_names)
                    if modelnum ~= length(model_names)
                        % get input separability
                        input_signals = getSig(td_test,glm_params{modelnum}.in_signals);
                        if endsWith(model_aliases{modelnum},'_actpasbaseline') || endsWith(model_aliases{modelnum},'_surprise')
                            % actpas input should have 100% separability
                            input_seps{modelnum} = 1;
                        else
                            input_seps{modelnum} = sum(predict(input_lda_mdl{modelnum},input_signals) == test_class)/length(test_class);
                        end

                        % calculate input margins
                        if get_margins
                            input_margin{modelnum} = (input_signals*input_lda_coeff{modelnum}(2:end)' + input_lda_coeff{modelnum}(1))/norm(input_lda_coeff{modelnum}(2:end));
                        end
                    end

                    % get separabilities from LDA models
                    % model_fr{modelnum} = cat(1,td_test.(model_names{modelnum}));
                    model_fr{modelnum} = getSig(td_test,{model_names{modelnum},which_units});

                    if do_population_stuff
                        % self_seps{modelnum} = sum(predict(lda_mdl{modelnum},model_fr{modelnum}) == test_class)/length(test_class);
                        % true_seps{modelnum} = sum(predict(lda_mdl{end},model_fr{modelnum}) == test_class)/length(test_class);
                        self_seps{modelnum} = sum(predict(lda_mdl{modelnum},(model_fr{modelnum}-pca_mu{modelnum})*pca_coeff{modelnum}) == test_class)/length(test_class);
                        true_seps{modelnum} = sum(predict(lda_mdl{end},(model_fr{modelnum}-pca_mu{end})*pca_coeff{end}) == test_class)/length(test_class);

                        % calculate margin from lda boundary
                        if get_margins
                            self_margin{modelnum} = ((model_fr{modelnum}-pca_mu{modelnum})*pca_coeff{modelnum}*lda_coeff{modelnum}(2:end)' + lda_coeff{modelnum}(1))/norm(lda_coeff{modelnum}(2:end));
                            true_margin{modelnum} = ((model_fr{modelnum}-pca_mu{end})*pca_coeff{end}*lda_coeff{end}(2:end)' + lda_coeff{end}(1))/norm(lda_coeff{end}(2:end));
                        end
                    end
                end

                % put together individual model evaluations
                foldEval{repeatnum,foldnum} = horzcat(...
                    makeNeuronTableStarter(...
                        td_train,...
                        struct('out_signal_names',unit_guide,'meta',struct('crossvalID',crossval_table.crossvalID))),...
                    model_eval{:},...
                    array2table(...
                        indiv_seps,...
                        'VariableNames',strcat(model_names,'_indiv_sep')));
                foldEval{repeatnum,foldnum}.Properties.VariableDescriptions((end-length(model_names)+1):end) = repmat({'linear'},1,length(model_names));

                % get model fr for trial table
                model_fr_table = table(model_fr{:},'VariableNames',[strcat(model_aliases,'_predFR') {neural_signals}]);
                model_fr_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_names));

                % assemble trial table
                trial_table_cell{repeatnum,foldnum} = horzcat(...
                    repmat(meta_table,length(trial_dir),1),...
                    repmat(crossval_table,length(trial_dir),1),...
                    trial_id,...
                    test_class_table,...
                    trial_dir_table,...
                    model_fr_table);

                if get_margins
                    % get margins for trial table
                    input_margin_table = table(input_margin{:},'VariableNames',strcat(model_aliases,'_input_margin'));
                    input_margin_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_aliases));
                    if do_population_stuff
                        trial_table_cell{repeatnum,foldnum} = horzcat(...
                            trial_table_cell{repeatnum,foldnum},...
                            input_margin_table,...
                            self_margin_table,...
                            true_margin_table);
                    else
                        self_margin_table = table(self_margin{:},'VariableNames',strcat([model_aliases {neural_signals}],'_self_margin'));
                        self_margin_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_names));
                        true_margin_table = table(true_margin{:},'VariableNames',strcat([model_aliases {neural_signals}],'_true_margin'));
                        true_margin_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_names));
                        trial_table_cell{repeatnum,foldnum} = horzcat(...
                            trial_table_cell{repeatnum,foldnum},...
                            input_margin_table);
                    end
                end

                % construct LDA table entry
                input_lda_coeff_table = cell2table(input_lda_coeff,'VariableNames',strcat(model_aliases,'_input_lda_coeff'));
                input_lda_coeff_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_aliases));
                lda_coeff_table = cell2table(lda_coeff,'VariableNames',strcat([model_aliases {neural_signals}],'_lda_coeff'));
                lda_coeff_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_names));
                pca_coeff_table = cell2table(pca_coeff,'VariableNames',strcat([model_aliases {neural_signals}],'_pca_coeff'));
                pca_coeff_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_names));
                pca_mu_table = cell2table(pca_mu,'VariableNames',strcat([model_aliases {neural_signals}],'_pca_mu'));
                pca_mu_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_names));
                input_sep_table = cell2table(input_seps,'VariableNames',strcat(model_aliases,'_input_sep'));
                input_sep_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_aliases));
                self_sep_table = cell2table(self_seps,'VariableNames',strcat([model_aliases {neural_signals}],'_self_sep'));
                self_sep_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_names));
                true_sep_table = cell2table(true_seps,'VariableNames',strcat([model_aliases {neural_signals}],'_true_sep'));
                true_sep_table.Properties.VariableDescriptions = repmat({'linear'},1,length(model_names));
                lda_table_cell{repeatnum,foldnum} = horzcat(...
                    meta_table,...
                    crossval_table,...
                    input_lda_coeff_table,...
                    lda_coeff_table,...
                    pca_coeff_table,...
                    pca_mu_table,...
                    input_sep_table,...
                    self_sep_table,...
                    true_sep_table);
        
                fprintf('\tEvaluated fold %d at time: %f\n',foldnum,toc(fold_tic))
            end
            fprintf('Evaluated repeat %d at time: %f\n',repeatnum,toc(repeat_tic))
        end

    %% Package results
        % get one separability table
        lda_table = vertcat(lda_table_cell{:});
        trial_table = vertcat(trial_table_cell{:});
        neuron_eval_table = vertcat(foldEval{:});

        results = struct('lda_table',lda_table,'trial_table',trial_table,'neuron_eval_table',neuron_eval_table);
