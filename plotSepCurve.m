%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes plots from results saved by calculateSepTimeCurve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up plotting variables
    datadir = '/data/raeed/project-data/limblab/s1-kinematics/Results/Separability';
    % filename = {'Han_20171101_TRT_encodingResults_run20180809.mat','Chips_20170915_TRT_encodingResults_run20180809.mat','Lando_20170802_encodingResults_run20180809.mat'};
    % files = dir(fullfile(datadir,'*separationResults_run20190228_rerun20190806.mat'));
    files = dir(fullfile(datadir,'*sepCurveResults_10pcs_run20190228_rerun20190903.mat'));
    filename = horzcat({files.name});
    
    % for figure saving
    figdir = '/home/raeed/Wiki/Projects/limblab/s1-kinematics/figures/Separability';
    run_date = char(datetime('today','format','yyyyMMdd'));

    monkey_names = {'Chips','Han'};
    neural_signals = 'S1_FR';
    models_to_plot = {neural_signals,'ext','extforce','handelbow','ext_actpasbaseline'};
    fr_names = {neural_signals,'ext_predFR','extforce_predFR','handelbow_predFR','ext_actpasbaseline_predFR'};
    model_titles = {'Actual Firing','Extrinsic','Extrinsic + Force','Hand/Elbow','Hand+Baseline'};
    num_pcs = 3;
    timeshifts = -20:5:19;

    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;

%% Compile information over all files
    [sep_table_cell, margin_corr_table_cell] = deal(cell(length(filename),1));
    fileclock = tic;
    fprintf('Started loading files...\n')
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        % get info from this file
        num_repeats = max(sepCurveResults.trial_table.crossvalID(:,1));
        num_folds = max(sepCurveResults.trial_table.crossvalID(:,2));

        % % extract data
        % % get separability
        % meta_cols = strcmpi(sepResults.lda_table.Properties.VariableDescriptions,'meta');
        % sep_cols = endsWith(sepResults.lda_table.Properties.VariableNames,'true_sep');
        % sep_table_cell{filenum} = sepResults.lda_table(:,meta_cols | sep_cols);

        % % extract margin correlations
        % meta_cols = strcmpi(sepResults.trial_table.Properties.VariableDescriptions,'meta');
        % margin_corr_entries = cell(num_repeats,num_folds);
        % for repeatnum = 1:num_repeats
        %     for foldnum = 1:num_folds
        %         [~,fold_trials] = getNTidx(sepResults.trial_table,'crossvalID',[repeatnum foldnum]);
        %         isActive = ~fold_trials.isPassive;

        %         % get S1 margin
        %         S1_margin = fold_trials.(strcat(neural_signals,'_self_margin'));

        %         % initialize
        %         [model_margin_corr, model_margin_corr_act, model_margin_corr_pas,model_input_margin_corr, model_input_margin_corr_act, model_input_margin_corr_pas] = deal(zeros(1,length(models_to_plot)));
        %         % compile model margins
        %         for modelnum = 1:length(models_to_plot)
        %             % get model margin
        %             model_margin = fold_trials.(sprintf('%s_self_margin',models_to_plot{modelnum}));

        %             model_margin_corr(modelnum) = corr(model_margin,S1_margin);
        %             model_margin_corr_act(modelnum) = corr(model_margin(isActive),S1_margin(isActive));
        %             model_margin_corr_pas(modelnum) = corr(model_margin(~isActive),S1_margin(~isActive));

        %             % get model input margin
        %             if strcmpi(models_to_plot{modelnum},neural_signals)
        %                 % Set it to S1 margin I guess...
        %                 model_input_margin = S1_margin;
        %             else
        %                 model_input_margin = fold_trials.(sprintf('%s_input_margin',models_to_plot{modelnum}));
        %             end

        %             model_input_margin_corr(modelnum) = corr(model_input_margin,S1_margin);
        %             model_input_margin_corr_act(modelnum) = corr(model_input_margin(isActive),S1_margin(isActive));
        %             model_input_margin_corr_pas(modelnum) = corr(model_input_margin(~isActive),S1_margin(~isActive));
        %         end
        %         model_margin_corr_table = array2table(model_margin_corr,...
        %             'VariableNames',strcat(models_to_plot,'_margin_corr'));
        %         model_margin_corr_act_table = array2table(model_margin_corr_act,...
        %             'VariableNames',strcat(models_to_plot,'_margin_corr_act'));
        %         model_margin_corr_pas_table = array2table(model_margin_corr_pas,...
        %             'VariableNames',strcat(models_to_plot,'_margin_corr_pas'));
        %         model_input_margin_corr_table = array2table(model_input_margin_corr,...
        %             'VariableNames',strcat(models_to_plot,'_input_margin_corr'));
        %         model_input_margin_corr_act_table = array2table(model_input_margin_corr_act,...
        %             'VariableNames',strcat(models_to_plot,'_input_margin_corr_act'));
        %         model_input_margin_corr_pas_table = array2table(model_input_margin_corr_pas,...
        %             'VariableNames',strcat(models_to_plot,'_input_margin_corr_pas'));
        %         margin_corr_entries{repeatnum,foldnum} = horzcat(fold_trials(1,'crossvalID'),...
        %             model_margin_corr_table,...
        %             model_margin_corr_act_table,...
        %             model_margin_corr_pas_table,...
        %             model_input_margin_corr_table,...
        %             model_input_margin_corr_act_table,...
        %             model_input_margin_corr_pas_table);
        %     end
        % end
        % margin_corr_table_cell{filenum} = join(sepResults.lda_table(:,meta_cols),vertcat(margin_corr_entries{:}));

        % % compose trial table for one crossval run
        % repeatnum = 3;
        % trial_table_cell = cell(num_folds,1);
        % for foldnum = 1:num_folds
        %     [~,trial_table_cell{foldnum}] = getNTidx(sepResults.trial_table,'crossvalID',[repeatnum foldnum]);
        % end
        % trial_table = vertcat(trial_table_cell{:});

        % % plot out neural population scatter
        % isActive = ~trial_table.isPassive;
        % [dirs,~,dir_idx] = unique(trial_table.trialDir);
        % dir_colors = linspecer(length(dirs));
        % figure('defaultaxesfontsize',18)
        % for modelnum = 1:length(models_to_plot)
        %     subplot(1,length(models_to_plot),modelnum)

        %     model_fr = trial_table.(fr_names{modelnum});
        %     [~,model_pca] = pca(model_fr);
        %     model_pca = model_pca(:,1:num_pcs);
        %     lda_mdl = fitcdiscr(model_pca,isActive);
        %     lda_vec = lda_mdl.Coeffs(2,1).Linear;
        %     lda_vec = lda_vec/sqrt(sum(lda_vec.^2)); % make into unit vector
        %     null_basis = null(lda_vec');

        %     % make 3d plot
        %     % scatter3(...
        %     %     model_pca(isActive,:)*lda_vec,...
        %     %     model_pca(isActive,:)*null_basis(:,1),...
        %     %     model_pca(isActive,:)*null_basis(:,2),...
        %     %     [],dir_colors(dir_idx(isActive),:),'filled')
        %     % hold on
        %     % scatter3(...
        %     %     model_pca(~isActive,:)*lda_vec,...
        %     %     model_pca(~isActive,:)*null_basis(:,1),...
        %     %     model_pca(~isActive,:)*null_basis(:,2),...
        %     %     [],dir_colors(dir_idx(~isActive),:))
        %     % % plot lines
        %     % plot3([0 0],ylim,[0 0],'--k','linewidth',2)
        %     % plot3([0 0],[0 0],zlim,'--k','linewidth',2)
        %     % axis equal

        %     % 2D plot
        %     % scatter(...
        %     %     model_pca(isActive,:)*lda_vec,...
        %     %     model_pca(isActive,:)*null_basis(:,1),...
        %     %     [],dir_colors(dir_idx(isActive),:),'filled')
        %     % hold on
        %     % scatter(...
        %     %     model_pca(~isActive,:)*lda_vec,...
        %     %     model_pca(~isActive,:)*null_basis(:,1),...
        %     %     [],dir_colors(dir_idx(~isActive),:))
        %     % plot([0 0],ylim,'--k','linewidth',2)
        %     % axis equal

        %     % 2D plot with no dir color
        %     scatter(...
        %         model_pca(isActive,:)*lda_vec,...
        %         model_pca(isActive,:)*null_basis(:,1),...
        %         [],'k','filled')
        %     hold on
        %     scatter(...
        %         model_pca(~isActive,:)*lda_vec,...
        %         model_pca(~isActive,:)*null_basis(:,1),...
        %         [],'r','filled')
        %     plot([0 0],ylim,'--k','linewidth',2)
        %     axis equal
        %     axis off

        %     title(model_titles{modelnum})
        % end
        % suptitle(sprintf('%s-%s',trial_table.monkey{1},trial_table.date_time{1}))

        % % plot out margin comparisons between models
        % S1_margin = trial_table.(strcat(neural_signals,'_self_margin'));
        % figure('defaultaxesfontsize',18)
        % for modelnum = 1:length(models_to_plot)
        %     subplot(1,length(models_to_plot),modelnum)

        %     % get model margin
        %     model_margin = trial_table.(sprintf('%s_self_margin',models_to_plot{modelnum}));
        %     
        %     % plot it out
        %     scatter(S1_margin(isActive),model_margin(isActive),[],'k','filled')
        %     % scatter(S1_margin(~isPassive),handelbow_input_margin(~isPassive),[],'k','filled')
        %     hold on
        %     scatter(S1_margin(~isActive),model_margin(~isActive),[],'r','filled')
        %     % scatter(S1_margin(isPassive),handelbow_input_margin(isPassive),[],'r','filled')
        %     plot(xlim,[0 0],'-k')
        %     plot([0 0],ylim,'-k')
        %     plot(xlim,xlim,'--k')
        %     set(gca,'box','off','tickdir','out','xlim',xlim,'ylim',xlim)
        %     axis equal
        %     title(sprintf('%s vs. S1 Margin comparison',getModelTitles(models_to_plot{modelnum})))
        %     ylabel(sprintf('%s model margin',getModelTitles(models_to_plot{modelnum})))
        %     xlabel('S1 FR margin')
        % end
        % suptitle(sprintf('%s-%s',trial_table.monkey{1},trial_table.date_time{1}))

        % % plot margin as a function of trial number (to look for trends)
        % figure('defaultaxesfontsize',18)
        % scatter(trial_table.trialID(isActive),S1_margin(isActive),[],'k','filled')
        % hold on
        % scatter(trial_table.trialID(~isActive),S1_margin(~isActive),[],'r','filled')
        % plot(xlim,[0 0],'-k','linewidth',0.5)
        % set(gca,'box','off','tickdir','out')
        % ylabel('S1 Margin')
        % xlabel('Trial number')
        % title(sprintf('%s-%s',trial_table.monkey{1},trial_table.date_time{1}))

        % Plot out separability as a function of time shift
        figure('defaultaxesfontsize',18)
        plot(10*timeshifts([1 end]),[0.5 0.5],'--k','linewidth',2)
        hold on
        plot(repmat(timeshifts*10,num_repeats*num_folds,1)',...
            sepCurveResults.lda_table.(strcat(neural_signals,'_self_sep'))',...
            'linewidth',0.5,'color','k');
        for modelnum = 2:length(models_to_plot)
            plot(repmat(timeshifts*10,num_repeats*num_folds,1)',...
                sepCurveResults.lda_table.(strcat(models_to_plot{modelnum},'_input_sep'))',...
                'linewidth',0.5,'color',getModelColors(models_to_plot{modelnum}));
        end
        set(gca,'box','off','tickdir','out','ylim',[0 1])
        xlabel('Time from movement onset')
        ylabel('Separability')
        title(sprintf('%s-%s',sepCurveResults.trial_table.monkey{1},sepCurveResults.trial_table.date_time{1}))

        figure('defaultaxesfontsize',18)
        avgTrial = neuronAverage(sepCurveResults.trial_table,struct(...
            'keycols',{{'monkey','date_time','task','trialID','isPassive'}}));
        [~,act_trials] = getNTidx(avgTrial,'isPassive',false);
        [~,pas_trials] = getNTidx(avgTrial,'isPassive',true);
        plot(10*timeshifts([1 end]),zeros(1,2),'--k','linewidth',2)
        hold on
        h = plot(repmat(timeshifts*10,height(act_trials),1)',...
            act_trials.(strcat(neural_signals,'_self_margin'))',...
            'linewidth',0.5,'color','k');
        [h.Color] = deal([0 0 0 0.5]);
        h = plot(repmat(timeshifts*10,height(pas_trials),1)',...
            pas_trials.(strcat(neural_signals,'_self_margin'))',...
            'linewidth',0.5,'color','r');
        [h.Color] = deal([1 0 0 0.5]);
        set(gca,'box','off','tickdir','out')
        xlabel('Time from movement onset')
        ylabel('Margin')
        title(sprintf('%s-%s',sepCurveResults.trial_table.monkey{1},sepCurveResults.trial_table.date_time{1}))

        % output a counter
        fprintf('Processed file %d of %d at time %f\n',filenum,length(filename),toc(fileclock))
    end
    sep_table = vertcat(sep_table_cell{:});
    margin_corr_table = vertcat(margin_corr_table_cell{:});

%% Make summary plots
    % plot session average connected by lines...
        figure('defaultaxesfontsize',18)
        alpha = 0.05;
        model_x = (2:3:((length(models_to_plot)-1)*3+2))/10;
        for monkeynum = 1:length(monkey_names)
            subplot(length(monkey_names),1,monkeynum)
            
            % figure out what sessions we have for this monkey
            [~,monkey_seps] = getNTidx(sep_table,'monkey',monkey_names{monkeynum});
            session_datetimes = unique(monkey_seps.date_time);

            for sessionnum = 1:length(session_datetimes)
                [~,session_seps] = getNTidx(monkey_seps,'date_time',session_datetimes{sessionnum});

                % estimate error bars
                [~,cols] = ismember(strcat(models_to_plot,'_true_sep'),session_seps.Properties.VariableNames);
                num_repeats = double(max(session_seps.crossvalID(:,1)));
                num_folds = double(max(session_seps.crossvalID(:,2)));
                crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
                yvals = mean(session_seps{:,cols});
                var_seps = var(session_seps{:,cols});
                upp = tinv(1-alpha/2,num_folds*num_repeats-1);
                low = tinv(alpha/2,num_folds*num_repeats-1);
                CI_lo = yvals + low * sqrt(crossval_correction*var_seps);
                CI_hi = yvals + upp * sqrt(crossval_correction*var_seps);
                
                % get value of true separability
                [~,S1_col] = ismember(strcat(neural_signals,'_true_sep'),session_seps.Properties.VariableNames);
                S1_val = mean(session_seps{:,S1_col});

                % plot dots and lines
                % plot(model_x',yvals','-','linewidth',0.5,'color',ones(1,3)*0.5)
                % plot(model_x',yvals','-','linewidth',0.5,'color',ones(1,3)*0.5)
                plot([min(model_x)-0.2 max(model_x)+0.2],[S1_val S1_val],'--','linewidth',2,'color',session_colors(sessionnum,:))
                hold on
                plot(repmat(model_x,2,1),[CI_lo;CI_hi],'-','linewidth',2,'color',session_colors(sessionnum,:))
                scatter(model_x(:),yvals(:),50,session_colors(sessionnum,:),'filled')
            end
            plot([min(model_x)-0.2 max(model_x)+0.2],[0.5 0.5],'--k','linewidth',2)
            ylabel('Separability (%)')
            set(gca,'box','off','tickdir','out',...
                'xlim',[min(model_x)-0.2 max(model_x)+0.2],...
                'xtick',model_x,'xticklabel',model_titles,...
                'ylim',[0.5 1],'ytick',[0 0.5 1])
        end
        % saveas(gcf,fullfile(figdir,sprintf('actpasSeparability_run%s.pdf',run_date)))

    % plot margin corr connected by lines...
        figure('defaultaxesfontsize',18)
        alpha = 0.05;
        model_x = (2:3:((length(models_to_plot)-1)*3+2))/10;
        for monkeynum = 1:length(monkey_names)
            % figure out what sessions we have for this monkey
            [~,monkey_margin_corr] = getNTidx(margin_corr_table,'monkey',monkey_names{monkeynum});
            session_datetimes = unique(monkey_margin_corr.date_time);

            margin_varnames_extra = {'','_act','_pas'};
            for marginnum = 1:length(margin_varnames_extra)
                % set subplot
                subplot(length(monkey_names),length(margin_varnames_extra),...
                    (monkeynum-1)*length(margin_varnames_extra)+marginnum)

                for sessionnum = 1:length(session_datetimes)
                    [~,session_margin_corr] = getNTidx(monkey_margin_corr,'date_time',session_datetimes{sessionnum});

                    % estimate error bars
                    [~,cols] = ismember(...
                        strcat(models_to_plot,'_margin_corr',margin_varnames_extra{marginnum}),...
                        session_margin_corr.Properties.VariableNames);
                    num_repeats = double(max(session_margin_corr.crossvalID(:,1)));
                    num_folds = double(max(session_margin_corr.crossvalID(:,2)));
                    crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
                    yvals = mean(session_margin_corr{:,cols});
                    var_margin_corr = var(session_margin_corr{:,cols});
                    upp = tinv(1-alpha/2,num_folds*num_repeats-1);
                    low = tinv(alpha/2,num_folds*num_repeats-1);
                    CI_lo = yvals + low * sqrt(crossval_correction*var_margin_corr);
                    CI_hi = yvals + upp * sqrt(crossval_correction*var_margin_corr);
                    
                    % plot dots and lines
                    % plot(model_x',yvals','-','linewidth',0.5,'color',ones(1,3)*0.5)
                    plot(repmat(model_x,2,1),[CI_lo;CI_hi],'-','linewidth',2,'color',session_colors(sessionnum,:))
                    hold on
                    scatter(model_x(:),yvals(:),50,session_colors(sessionnum,:),'filled')
                end
                plot([min(model_x)-0.2 max(model_x)+0.2],[0 0],'--k','linewidth',2)
                set(gca,'box','off','tickdir','out',...
                    'xlim',[min(model_x)-0.2 max(model_x)+0.2],...
                    'xtick',model_x,'xticklabel',model_titles,'xticklabelrotation',45,...
                    'ylim',[-1 1],'ytick',[-1 0 1])
            end
        end
        % add labels
        for monkeynum = 1:length(monkey_names)
            % set subplot
            subplot(length(monkey_names),length(margin_varnames_extra),...
                (monkeynum-1)*length(margin_varnames_extra)+1)

            ylabel({monkey_names{monkeynum};'Margin Correlation'})
        end
        for marginnum = 1:length(margin_varnames_extra)
            % set subplot
            subplot(length(monkey_names),length(margin_varnames_extra),...
                marginnum)

            title(sprintf('MarginCorrelation%s',margin_varnames_extra{marginnum}),'interpreter','none')
        end
        % saveas(gcf,fullfile(figdir,sprintf('actpasSeparability_run%s.pdf',run_date)))

