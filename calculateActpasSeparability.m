%% Get act vs pas
    [~,td] = getTDidx(trial_data,'result','R');

    % Remove unsorted channels
    keepers = (td(1).S1_unit_guide(:,2)~=0);
    for trial = 1:length(td)
        td(trial).S1_unit_guide = td(trial).S1_unit_guide(keepers,:);
        td(trial).S1_spikes = td(trial).S1_spikes(:,keepers);
    end

    % remove low firing neurons
    td = removeBadNeurons(td,struct('min_fr',0.1));
    
    td = getMoveOnsetAndPeak(td,struct('start_idx','idx_goCueTime','end_idx','idx_endTime','method','peak','min_ds',1));
    td = smoothSignals(td,struct('signals','markers'));
    td = getDifferential(td,struct('signal','markers','alias','marker_vel'));
    % add firing rates rather than spike counts
    td = addFiringRates(td,struct('array','S1'));
    
    [~,td_act] = getTDidx(td,'ctrHoldBump',false);
    % clean nans out...?
    nanners = isnan(cat(1,td_act.target_direction));
    td_act = td_act(~nanners);
    td_act = trimTD(td_act,{'idx_movement_on',0},{'idx_movement_on',14});
    td_act = binTD(td_act,'average');
    [~,td_pas] = getTDidx(td,'ctrHoldBump',true);
    td_pas = trimTD(td_pas,{'idx_bumpTime',0},{'idx_bumpTime',14});
    td_pas = binTD(td_pas,'average');

    % even out sizes
    minsize = min(length(td_act),length(td_pas));
    td_act = td_act(1:minsize);
    td_pas = td_pas(1:minsize);
    td = cat(2,td_act,td_pas);

%% get splits on array
    if exist('ant_chans','var')
        % array_signals = find(ismember(td(1).S1_unit_guide(:,1),ant_chans));
        array_signals = find(~ismember(td(1).S1_unit_guide(:,1),ant_chans));
        % array_signals = sort(randperm(length(td(1).S1_unit_guide(:,1)),floor(length(td(1).S1_unit_guide(:,1))/2))');
        % array_signals = 1:length(td(1).S1_unit_guide);
    else
        array_signals = 1:length(td(1).S1_unit_guide);
    end

%% suppress getTDfields warning
    getTDfields(td,'time');
    onetime_warn = warning('query','last'); 
    warning('off',onetime_warn.identifier)

%% crossvalidate separabilities
    num_folds = 5;
    num_repeats = 20;
    num_pcs = 5;
    hand_idx = 1:3;
    elbow_idx = 28:30;
    
    sep_true = zeros(num_repeats,num_folds);
    sep_ext = zeros(num_repeats,num_folds);
    sep_forcekin = zeros(num_repeats,num_folds);
    sep_handelbow = zeros(num_repeats,num_folds);
    
    repeat_tic = tic;
    fprintf('Starting %dx%d cross-validation at time: %f\n',num_repeats,num_folds,toc(repeat_tic))
    for repeatnum = 1:num_repeats
        foldidx = crossvalind('kfold',minsize,num_folds);
        fold_tic = tic;
        for foldnum = 1:num_folds
            % do the foldy thing
            train_idx = (foldidx~=foldnum);
            td_train = cat(2,td_act(train_idx),td_pas(train_idx));
            % % Hand kinematics model
            % [~,ext_model_info] = getModel(td_train,struct('model_type','glm',...
            %     'model_name','ext','in_signals',{{'pos','vel'}},...
            %     'out_signals',{{'S1_FR',array_signals}}));
            % % Hand forcekin model
            % [~,forcekin_model_info] = getModel(td_train,struct('model_type','glm',...
            %     'model_name','forcekin','in_signals',{{'pos','all';'vel','all';'force','all'}},...
            %     'out_signals',{{'S1_FR',array_signals}}));
            % Hand kinematics model
            [~,ext_model_info] = getModel(td_train,struct('model_type','glm',...
                'model_name','ext','in_signals',{{'markers',hand_idx;'marker_vel',hand_idx}},...
                'out_signals',{{'S1_FR',array_signals}}));
            % Hand forcekin model
            [~,forcekin_model_info] = getModel(td_train,struct('model_type','glm',...
                'model_name','forcekin','in_signals',{{'markers',hand_idx;'marker_vel',hand_idx;'force','all'}},...
                'out_signals',{{'S1_FR',array_signals}}));
            % Hand/Elbow model
            [~,handelbow_model_info] = getModel(td_train,struct('model_type','glm',...
                'model_name','handelbow','in_signals',{{'markers',[hand_idx elbow_idx];'marker_vel',[hand_idx elbow_idx]}},...
                'out_signals',{{'S1_FR',array_signals}}));
    
            % apply to test set
            test_idx = (foldidx==foldnum);
            td_test = cat(2,td_act(test_idx),td_pas(test_idx));
            td_test = getModel(td_test,ext_model_info);
            td_test = getModel(td_test,forcekin_model_info);
            td_test = getModel(td_test,handelbow_model_info);
    
            % get pcas
            td_train = sqrtTransform(td_train,'S1_FR');
            [td_train,pca_info] = dimReduce(td_train,struct('signals',{{'S1_FR',array_signals}}));
    
            td_test = sqrtTransform(td_test,'S1_FR');
            td_test = sqrtTransform(td_test,'glm_ext');
            td_test = sqrtTransform(td_test,'glm_forcekin');
            td_test = sqrtTransform(td_test,'glm_handelbow');
    
            pca_info_ext = pca_info;
            pca_info_ext.signals = {'glm_ext',array_signals};
            pca_info_ext.recenter_for_proj = true;
            pca_info_forcekin = pca_info;
            pca_info_forcekin.signals = {'glm_forcekin',array_signals};
            pca_info_forcekin.recenter_for_proj = true;
            pca_info_handelbow = pca_info;
            pca_info_handelbow.signals = {'glm_handelbow',array_signals};
            pca_info_handelbow.recenter_for_proj = true;
    
            td_test = dimReduce(td_test,pca_info);
            td_test = dimReduce(td_test,pca_info_ext);
            td_test = dimReduce(td_test,pca_info_forcekin);
            td_test = dimReduce(td_test,pca_info_handelbow);
    
            % get LDA models
            [~,actual_lda] = test_sep(td_train,struct('signals',{{'S1_FR_pca',1:num_pcs}}));
            % [~,ext_lda] = test_sep(td_train,struct('signals',{{'glm_ext_pca',1:5}}));
            % [~,forcekin_lda] = test_sep(td_train,struct('signals',{{'glm_forcekin_pca',1:5}}));
            % [~,handelbow_lda] = test_sep(td_train,struct('signals',{{'glm_handelbow_pca',1:5}}));
    
            % check LDA models
            if repeatnum==num_repeats
                do_plot = true;
                if foldnum==1
                    h1 = figure(1);
                    h2 = figure(2);
                    h3 = figure(3);
                    h4 = figure(4);
                end
            else
                do_plot = false;
                h1 = [];
                h2 = [];
                h3 = [];
                h4 = [];
            end
            sep_true(repeatnum,foldnum) = test_sep(td_test,struct('signals',{{'S1_FR_pca',1:num_pcs}},'mdl',actual_lda,'do_plot',do_plot,'fig_handle',h1));
            sep_ext(repeatnum,foldnum) = test_sep(td_test,struct('signals',{{'glm_ext_pca',1:num_pcs}},'mdl',actual_lda,'do_plot',do_plot,'fig_handle',h2));
            sep_forcekin(repeatnum,foldnum) = test_sep(td_test,struct('signals',{{'glm_forcekin_pca',1:num_pcs}},'mdl',actual_lda,'do_plot',do_plot,'fig_handle',h3));
            sep_handelbow(repeatnum,foldnum) = test_sep(td_test,struct('signals',{{'glm_handelbow_pca',1:num_pcs}},'mdl',actual_lda,'do_plot',do_plot,'fig_handle',h4));
    
            fprintf('\tEvaluated fold %d at time: %f\n',foldnum,toc(fold_tic))
        end
        fprintf('Evaluated repeat %d at time: %f\n',repeatnum,toc(repeat_tic))
    end

%% package results
    % get one separability table
    sep_true = sep_true(:);
    sep_ext = sep_ext(:);
    sep_forcekin = sep_forcekin(:);
    sep_handelbow = sep_handelbow(:);

    seps = table(sep_true,sep_ext,sep_forcekin,sep_handelbow,...
        'VariableNames',{'true','ext','forcekin','handelbow'});

%% bar plot for separabilities of models
    % bar colors
    bar_colors(1,:) = [0.5 0.5 0.5];
    bar_colors(2,:) = [247, 148, 30]/255;
    bar_colors(3,:) = [247, 192, 30]/255;
    bar_colors(4,:) = [193, 25, 47]/255;

    % make plot
    datadir = '~/Projects/limblab/data-td/ForceKin/Results/separability';
    filenames = {'Han_20170203_actpasSeps_run20180818.mat','Chips_20170913_actpasSeps_run20180818.mat'};
    num_monks = length(filenames);
    monk_x = (2:3:((num_monks-1)*3+2))/10;
    template_x = linspace(-0.5,0.5,4)/10;
    model_spacing = mode(diff(template_x));
    figure('defaultaxesfontsize',18)
    for monkeynum = 1:num_monks
        load(fullfile(datadir,filenames{monkeynum}))

        % calculate stats
        mean_seps = mean(seps{:,:});
        var_seps = var(seps{:,:});
        correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
        std_err_seps = sqrt(var_seps*correction);

        for modelnum = 1:length(mean_seps)
            xval = monk_x(monkeynum) + template_x(modelnum);
            bar(xval,mean_seps(modelnum),model_spacing,'facecolor',bar_colors(modelnum,:),'edgecolor','none')
            hold on
            plot([xval xval],[mean_seps(modelnum)-std_err_seps(modelnum) mean_seps(modelnum)+std_err_seps(modelnum)],'k','linewidth',3)
        end
    end
    plot([0 monk_x(end)+0.2],[0.5 0.5],'--k','linewidth',2)
    plot([0 monk_x(end)+0.2],[1 1],'--k','linewidth',2)
    set(gca,'box','off','tickdir','out',...
        'xtick',monk_x,'xticklabel',filenames,'xlim',[0 monk_x(end)+0.2],...
        'ytick',[0 0.5 1],'yticklabel',{'','50%','100%'},...
        'ticklabelinterpreter','none')

%% Turn warnings back on
    warning('on',onetime_warn.identifier)
