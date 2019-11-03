%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes plots from results saved by calculateActPasSeparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up plotting variables
    datadir = '/data/raeed/project-data/limblab/s1-kinematics/Results/Separability';
    % filename = {'Han_20171101_TRT_encodingResults_run20180809.mat','Chips_20170915_TRT_encodingResults_run20180809.mat','Lando_20170802_encodingResults_run20180809.mat'};
    % files = dir(fullfile(datadir,'*separationResults_run20190228_rerun20190806.mat'));
    % files = dir(fullfile(datadir,'*separationResults_run20190228.mat'));
    % files = dir(fullfile(datadir,'*separationResults_50msLag_run20190228_rerun20190809.mat'));
    % files = dir(fullfile(datadir,'*separationResults_run20190228_rerun20191010.mat'));
    files = dir(fullfile(datadir,'*separationResults_run20190228_rerun20191016.mat'));
    filename = horzcat({files.name});
    
    % for figure saving
    figdir = '/home/raeed/Wiki/Projects/limblab/s1-kinematics/figures/Separability';
    run_date = char(datetime('today','format','yyyyMMdd'));

    monkey_names = {'Chips','Han'};
    neural_signals = 'S1_FR';
    % models_to_plot = {neural_signals,'ext','extforce','handelbow','ext_actpasbaseline'};
    models_to_plot = {neural_signals,'ext','extforce','handelbow'};
    fr_names = {neural_signals,'ext_predFR','extforce_predFR','handelbow_predFR','ext_actpasbaseline_predFR'};
    model_titles = {'Actual Firing','Extrinsic','Extrinsic + Force','Hand/Elbow','Hand+Baseline'};
    num_pcs = 3;

    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;

%% Single neuron figures/analyses
    fileclock = tic;
    neuron_eval_cell = cell(length(filename),1);
    fprintf('Started loading files...\n')
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        % replace infs with nans
        numeric_cols = strcmpi(sepResults.neuron_eval_table.Properties.VariableDescriptions,'linear');
        numeric_vals = sepResults.neuron_eval_table(:,numeric_cols).Variables;
        infidx = isinf(numeric_vals);
        numeric_vals(infidx) = NaN;
        sepResults.neuron_eval_table(:,numeric_cols).Variables = numeric_vals;

        % start with the tables we need...
        average_trials = neuronAverage(...
            sepResults.trial_table,...
            struct('keycols',{{'monkey','date_time','task','trialID','isPassive'}},'do_ci',false));
        avg_neuron_eval = neuronAverage(sepResults.neuron_eval_table,struct(...
            'keycols',{{'monkey','date','task','signalID'}},...
            'do_ci',false,...
            'do_nanmean',true));

        % try a histogram version of neuron v firing rate
        figure('defaultaxesfontsize',18)
        ax = zeros(size(average_trials.S1_FR,2),1);
        % here we know that there are only a finite number of possibilities
        possible_FR = unique(average_trials.S1_FR);
        % split into active and passive
        [~,act_trials] = getNTidx(average_trials,'isPassive',false);
        [~,pas_trials] = getNTidx(average_trials,'isPassive',true);
        for neuronnum = 1:size(average_trials.S1_FR,2)
            ax(neuronnum) = subplot(1,size(average_trials.S1_FR,2),neuronnum);

            % get counts of fr in the unique bins
            act_counts = histcounts(act_trials.S1_FR(:,neuronnum),[possible_FR;Inf]);
            pas_counts = histcounts(pas_trials.S1_FR(:,neuronnum),[possible_FR;Inf]);

            % plot bars for each
            % plot([0 0],[possible_FR(1) possible_FR(end)])
            barh(possible_FR',act_counts,1,'FaceColor','k','EdgeColor','none','FaceAlpha',0.5)
            hold on
            barh(possible_FR',pas_counts,1,'FaceColor','r','EdgeColor','none','FaceAlpha',0.5)

            % print out separabilities...
            title(sprintf('%2.0f',avg_neuron_eval.S1_FR_indiv_sep(neuronnum,:)*100))

            % set(gca,'box','off','tickdir','out')
            axis off
        end
        subplot(1,size(average_trials.S1_FR,2),1)
        axis on
        set(gca,'box','off','tickdir','out','xtick',[])
        ylabel('Firing rate (Hz)')
        suptitle(sprintf('%s %s',sepResults.neuron_eval_table.monkey{1},sepResults.neuron_eval_table.date{1}))
        linkaxes(ax,'y')
        saveas(gcf,fullfile(figdir,sprintf(...
            '%s_%s_actpasFR_run%s.pdf',...
            sepResults.neuron_eval_table.monkey{1},...
            strrep(sepResults.neuron_eval_table.date{1},'/',''),...
            run_date)))

        % plot separability of each neuron
        figure('defaultaxesfontsize',18)
        signalIDs = avg_neuron_eval.signalID;
        plot([0 size(signalIDs,1)+1],[0.5 0.5],'--k','linewidth',2)
        hold on
        for neuronnum = 1:size(signalIDs,1)
            [~,single_neuron_eval] = getNTidx(sepResults.neuron_eval_table,'signalID',signalIDs(neuronnum,:));
            scatter(repmat(neuronnum,1,height(single_neuron_eval)),single_neuron_eval.S1_FR_indiv_sep,25,'k','filled','markerfacealpha',0.2)
            scatter(neuronnum,avg_neuron_eval.S1_FR_indiv_sep(neuronnum,:),100,'k','filled')
        end
        set(gca,'box','off','tickdir','out','ylim',[0 1],'xlim',[0 size(signalIDs,1)+1])
        saveas(gcf,fullfile(figdir,sprintf(...
            '%s_%s_indivneuron_separability_run%s.pdf',...
            sepResults.neuron_eval_table.monkey{1},...
            strrep(sepResults.neuron_eval_table.date{1},'/',''),...
            run_date)))

        % compile neuron eval table together
        neuron_eval_cell{filenum} = sepResults.neuron_eval_table;
        
        % output a counter
        fprintf('Processed file %d of %d at time %f\n',filenum,length(filename),toc(fileclock))
    end

    % extract things
    neuron_eval = vertcat(neuron_eval_cell{:});
    neuron_eval.Properties.VariableNames = strrep(neuron_eval.Properties.VariableNames,'glm_','');
    neuron_eval.Properties.VariableNames = strrep(neuron_eval.Properties.VariableNames,'model_','');
    avg_neuron_eval = neuronAverage(neuron_eval,struct(...
        'keycols',{{'monkey','date','task','signalID'}},...
        'do_ci',false,...
        'do_nanmean',true));

    % compare pR2 of handelbow vs ext
    figure('defaultaxesfontsize',18)
    model_pairs = {'ext','handelbow';'extforce','handelbow'};
    for pairnum = 1:size(model_pairs,1)
        for monkeynum = 1:length(monkey_names)
            % set subplot...
            subplot(length(monkey_names),size(model_pairs,1),...
                (monkeynum-1)*size(model_pairs,1)+pairnum)
            plot([-1 1],[-1 1],'k--','linewidth',0.5)
            hold on
            plot([0 0],[-1 1],'k-','linewidth',0.5)
            plot([-1 1],[0 0],'k-','linewidth',0.5)

            % get sessions
            [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
            session_dates = unique(monkey_evals.date);

            for sessionnum = 1:length(session_dates)
                [~,session_evals] = getNTidx(monkey_evals,'date',session_dates{sessionnum});
                pr2_winners = compareEncoderMetrics(session_evals,struct(...
                    'bonferroni_correction',6,...
                    'models',{models_to_plot},...
                    'model_pairs',{model_pairs},...
                    'postfix','_eval'));

                [~,avg_pR2] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});
                % scatter filled circles if there's a winner, empty circles if not
                no_winner =  cellfun(@isempty,pr2_winners(pairnum,:));
                scatter(...
                    avg_pR2.(strcat(model_pairs{pairnum,1},'_eval'))(no_winner),...
                    avg_pR2.(strcat(model_pairs{pairnum,2},'_eval'))(no_winner),...
                    [],session_colors(sessionnum,:))
                scatter(...
                    avg_pR2.(strcat(model_pairs{pairnum,1},'_eval'))(~no_winner),...
                    avg_pR2.(strcat(model_pairs{pairnum,2},'_eval'))(~no_winner),...
                    [],session_colors(sessionnum,:),'filled')
            end

            % make axes pretty
            % set(gca,'box','off','tickdir','out',...
            %     'xlim',[-0.3 0.6],'ylim',[-0.3 0.6])
            axis image
            if monkeynum ~= 1 || pairnum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,1})))
            ylabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,2})))
        end
    end
    saveas(gcf,fullfile(figdir,sprintf('actpas_indivneuron_modelpR2compare_run%s.pdf',run_date)))

    % Plot within condition vs across condition pR2 for each neuron in all sessions
    conds = {'act','pas'};
    for modelnum = 2:length(models_to_plot)
        figure('defaultaxesfontsize',18)
        for monkeynum = 1:length(monkey_names)
            for condnum = 1:2
                % set subplot
                subplot(length(monkey_names),2,(monkeynum-1)*2+condnum)
                plot([-1 1],[-1 1],'k--','linewidth',0.5)
                hold on
                plot([0 0],[-1 1],'k-','linewidth',0.5)
                plot([-1 1],[0 0],'k-','linewidth',0.5)

                % get sessions
                [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
                session_dates = unique(monkey_evals.date);

                % plot out each session
                for sessionnum = 1:length(session_dates)
                    [~,avg_pR2] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});
                    scatter(...
                        avg_pR2.(sprintf('%s_%s_eval',models_to_plot{modelnum},conds{condnum})),...
                        avg_pR2.(sprintf('%s_eval',models_to_plot{modelnum})),...
                        [],session_colors(sessionnum,:),'filled')
                end
                % make axes pretty
                % set(gca,'box','off','tickdir','out',...
                %     'xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
                axis image
                if monkeynum ~= 1 || condnum ~= 1
                    set(gca,'box','off','tickdir','out',...
                        'xtick',[],'ytick',[])
                end
                xlabel(sprintf('%s %s pR2',getModelTitles(models_to_plot{modelnum}),conds{condnum}))
                ylabel(sprintf('%s full pR2',getModelTitles(models_to_plot{modelnum})))
            end
        end
        suptitle('Full pR^2 vs within condition pR^2')
        saveas(gcf,fullfile(figdir,sprintf('actpas_%s_pr2_full_v_within_run%s.pdf',models_to_plot{modelnum},run_date)))
    end

    % get correlation values for each crossval
    keycols = {'monkey','date','task','crossvalID'};
    keyTable = unique(neuron_eval(:,keycols));
    corr_cell = cell(height(keyTable),1);
    for key_idx = 1:height(keyTable)
        key = keyTable(key_idx,:);
        cond_idx = ismember(neuron_eval(:,keycols),key);
        neuron_eval_select = neuron_eval(cond_idx,:);

        % get correlations
        model_corr = cell(1,length(models_to_plot)-1);
        for modelnum = 2:length(models_to_plot)
            vaf_modelsep_neuronsep = ...
                1 - ...
                sum((neuron_eval_select.S1_FR_indiv_sep-neuron_eval_select.(sprintf('%s_indiv_sep',models_to_plot{modelnum}))).^2)/...
                sum((neuron_eval_select.S1_FR_indiv_sep-mean(neuron_eval_select.S1_FR_indiv_sep)).^2);
            corr_modelsep_neuronsep = corr(neuron_eval_select.S1_FR_indiv_sep,neuron_eval_select.(sprintf('%s_indiv_sep',models_to_plot{modelnum})),'rows','complete');
            corr_pr2_neuronsep = corr(neuron_eval_select.S1_FR_indiv_sep,neuron_eval_select.(sprintf('%s_eval',models_to_plot{modelnum})),'rows','complete');
            corr_actpr2_neuronsep = corr(neuron_eval_select.S1_FR_indiv_sep,neuron_eval_select.(sprintf('%s_act_eval',models_to_plot{modelnum})),'rows','complete');
            corr_paspr2_neuronsep = corr(neuron_eval_select.S1_FR_indiv_sep,neuron_eval_select.(sprintf('%s_pas_eval',models_to_plot{modelnum})),'rows','complete');

            model_corr{modelnum-1} = table(...
                vaf_modelsep_neuronsep,...
                corr_modelsep_neuronsep,...
                corr_pr2_neuronsep,...
                corr_actpr2_neuronsep,...
                corr_paspr2_neuronsep,...
                'VariableNames',strcat(models_to_plot{modelnum},{'_vaf_modelsep_neuronsep','_corr_modelsep_neuronsep','_corr_pr2_neuronsep','_corr_actpr2_neuronsep','_corr_paspr2_neuronsep'}));
        end

        % put together in table
        corr_cell{key_idx} = horzcat(model_corr{:});
    end
    neuron_corr_table = horzcat(keyTable,vertcat(corr_cell{:}));

    % plot separability against full and within condition pR2
    conds = {'','act_','pas_'};
    for modelnum = 2:length(models_to_plot)
        figure('defaultaxesfontsize',18)
        for monkeynum = 1:length(monkey_names)
            for condnum = 1:length(conds)
                % set subplot
                subplot(length(monkey_names),length(conds),(monkeynum-1)*length(conds)+condnum)
                plot([0 0],[0 1],'k-','linewidth',0.5)
                hold on
                plot([-1 1],[0 0],'k-','linewidth',0.5)
                plot([-1 1],[0.5 0.5],'k--','linewidth',0.5)

                % get sessions
                [~,monkey_evals] = getNTidx(neuron_eval,'monkey',monkey_names{monkeynum});
                session_dates = unique(monkey_evals.date);

                % plot out each session
                for sessionnum = 1:length(session_dates)
                    [~,avg_pR2] = getNTidx(avg_neuron_eval,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});
                    scatter(...
                        avg_pR2.(sprintf('%s_%seval',models_to_plot{modelnum},conds{condnum})),...
                        avg_pR2.(sprintf('S1_FR_indiv_sep',models_to_plot{modelnum})),...
                        [],session_colors(sessionnum,:),'filled')

                    % fit quick linear model to plot fit line
                    lm = fitlm(...
                        avg_pR2.(sprintf('%s_%seval',models_to_plot{modelnum},conds{condnum})),...
                        avg_pR2.S1_FR_indiv_sep);
                    plot([-1;1],lm.predict([-1;1]),...
                        '--','color',session_colors(sessionnum,:),'linewidth',1)
                end
                % make axes pretty
                set(gca,'box','off','tickdir','out',...
                    'xlim',[-1 1],'ylim',[0.4 1])
                if monkeynum ~= 1 || condnum ~= 1
                    set(gca,'box','off','tickdir','out',...
                        'xtick',[],'ytick',[])
                end
                xlabel(sprintf('%s %s pR2',getModelTitles(models_to_plot{modelnum}),conds{condnum}))
                ylabel(sprintf('Neural Separability',getModelTitles(models_to_plot{modelnum})))
            end
        end
        suptitle('Neural separability vs pR^2')
        saveas(gcf,fullfile(figdir,sprintf('actpas_%s_sepvpr2_run%s.pdf',models_to_plot{modelnum},run_date)))
    end

    % make figure for VAF of model separability
    figure('defaultaxesfontsize',18)
    alpha = 0.05;
    model_x = (2:3:((length(models_to_plot)-2)*3+2))/10;
    for monkeynum = 1:length(monkey_names)
        subplot(1,length(monkey_names),monkeynum)
        plot([min(model_x)-0.2 max(model_x)+0.2],[0 0],'-k','linewidth',2)
        hold on
        
        % figure out what sessions we have for this monkey
        [~,monkey_corrs] = getNTidx(neuron_corr_table,'monkey',monkey_names{monkeynum});
        session_dates = unique(monkey_corrs.date);

        for sessionnum = 1:length(session_dates)
            [~,session_corrs] = getNTidx(monkey_corrs,'date',session_dates{sessionnum});

            % estimate error bars
            [~,cols] = ismember(strcat(models_to_plot(2:end),'_vaf_modelsep_neuronsep'),session_corrs.Properties.VariableNames);
            num_repeats = double(max(session_corrs.crossvalID(:,1)));
            num_folds = double(max(session_corrs.crossvalID(:,2)));
            crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
            yvals = mean(session_corrs{:,cols});
            var_corrs = var(session_corrs{:,cols});
            upp = tinv(1-alpha/2,num_folds*num_repeats-1);
            low = tinv(alpha/2,num_folds*num_repeats-1);
            CI_lo = yvals + low * sqrt(crossval_correction*var_corrs);
            CI_hi = yvals + upp * sqrt(crossval_correction*var_corrs);

            % plot dots and lines
            plot(repmat(model_x,2,1),[CI_lo;CI_hi],'-','linewidth',2,'color',session_colors(sessionnum,:))
            scatter(model_x(:),yvals(:),50,session_colors(sessionnum,:),'filled')
        end
        ylabel('VAF of model predicted separability')
        title(monkey_names{monkeynum})
        set(gca,'box','off','tickdir','out',...
            'xlim',[min(model_x)-0.2 max(model_x)+0.2],...
            'xtick',model_x,'xticklabel',getModelTitles(models_to_plot(2:end)),...
            'ylim',[-1 1],'ytick',[-1 0 0.5 1])
    end
    saveas(gcf,fullfile(figdir,sprintf('actpas_indivneuron_modelseparabilityVAF_run%s.pdf',run_date)))

    % make figure for correlations of model separability with actual separability
    figure('defaultaxesfontsize',18)
    alpha = 0.05;
    model_x = (2:3:((length(models_to_plot)-2)*3+2))/10;
    for monkeynum = 1:length(monkey_names)
        subplot(1,length(monkey_names),monkeynum)
        plot([min(model_x)-0.2 max(model_x)+0.2],[0 0],'-k','linewidth',2)
        hold on
        
        % figure out what sessions we have for this monkey
        [~,monkey_corrs] = getNTidx(neuron_corr_table,'monkey',monkey_names{monkeynum});
        session_dates = unique(monkey_corrs.date);

        for sessionnum = 1:length(session_dates)
            [~,session_corrs] = getNTidx(monkey_corrs,'date',session_dates{sessionnum});

            % estimate error bars
            [~,cols] = ismember(strcat(models_to_plot(2:end),'_corr_modelsep_neuronsep'),session_corrs.Properties.VariableNames);
            num_repeats = double(max(session_corrs.crossvalID(:,1)));
            num_folds = double(max(session_corrs.crossvalID(:,2)));
            crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
            yvals = mean(session_corrs{:,cols});
            var_corrs = var(session_corrs{:,cols});
            upp = tinv(1-alpha/2,num_folds*num_repeats-1);
            low = tinv(alpha/2,num_folds*num_repeats-1);
            CI_lo = yvals + low * sqrt(crossval_correction*var_corrs);
            CI_hi = yvals + upp * sqrt(crossval_correction*var_corrs);

            % plot dots and lines
            plot(repmat(model_x,2,1),[CI_lo;CI_hi],'-','linewidth',2,'color',session_colors(sessionnum,:))
            scatter(model_x(:),yvals(:),50,session_colors(sessionnum,:),'filled')
        end
        ylabel('Correlation between model predicted separability and actual neural separability')
        title(monkey_names{monkeynum})
        set(gca,'box','off','tickdir','out',...
            'xlim',[min(model_x)-0.2 max(model_x)+0.2],...
            'xtick',model_x,'xticklabel',getModelTitles(models_to_plot(2:end)),...
            'ylim',[-1 1],'ytick',[-1 0 0.5 1])
    end
    saveas(gcf,fullfile(figdir,sprintf('actpas_indivneuron_modelseparabilityCorr_run%s.pdf',run_date)))

    % make figure for correlations of pR2 handelbow model with separability
    alpha = 0.05;
    xvals = [2 5 8]/10;
    for modelnum = 2:length(models_to_plot)
        figure('defaultaxesfontsize',18)
        for monkeynum = 1:length(monkey_names)
            subplot(1,length(monkey_names),monkeynum)
            plot([min(xvals)-0.2 max(xvals)+0.2],[0 0],'-k','linewidth',2)
            hold on
            
            % figure out what sessions we have for this monkey
            [~,monkey_corrs] = getNTidx(neuron_corr_table,'monkey',monkey_names{monkeynum});
            session_dates = unique(monkey_corrs.date);

            for sessionnum = 1:length(session_dates)
                [~,session_corrs] = getNTidx(monkey_corrs,'date',session_dates{sessionnum});

                % estimate error bars
                [~,cols] = ismember(...
                    strcat(models_to_plot{modelnum},{'_corr_pr2_neuronsep','_corr_actpr2_neuronsep','_corr_paspr2_neuronsep'}),...
                    session_corrs.Properties.VariableNames);
                num_repeats = double(max(session_corrs.crossvalID(:,1)));
                num_folds = double(max(session_corrs.crossvalID(:,2)));
                crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
                yvals = mean(session_corrs{:,cols});
                var_corrs = var(session_corrs{:,cols});
                upp = tinv(1-alpha/2,num_folds*num_repeats-1);
                low = tinv(alpha/2,num_folds*num_repeats-1);
                CI_lo = yvals + low * sqrt(crossval_correction*var_corrs);
                CI_hi = yvals + upp * sqrt(crossval_correction*var_corrs);
                
                % plot dots and lines
                plot(repmat(xvals,2,1),[CI_lo;CI_hi],'-','linewidth',2,'color',session_colors(sessionnum,:))
                scatter(xvals(:),yvals(:),50,session_colors(sessionnum,:),'filled')
            end
            ylabel('Correlation with active/passive separability')
            title(monkey_names{monkeynum})
            set(gca,'box','off','tickdir','out',...
                'xlim',[min(xvals)-0.2 max(xvals)+0.2],...
                'xtick',xvals,'xticklabel',{'Full pR^2','Active pR^2','Passive pR^2'},...
                'ylim',[-1 1],'ytick',[-1 -0.5 0 0.5 1])
        end
        suptitle(models_to_plot{modelnum})
        saveas(gcf,fullfile(figdir,sprintf('actpas_indivneuron_%s_pr2separabilityCorr_run%s.pdf',models_to_plot{modelnum},run_date)))
    end

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

    % plot individual neural separability on array map
        avg_neuron_eval = horzcat(avg_neuron_eval,table(avg_neuron_eval.signalID(:,1),'VariableNames',{'chan'}));
        figure('defaultaxesfontsize',18)
        lm_table = plotArrayMap(avg_neuron_eval,struct('map_plot','S1_FR_indiv_sep','clims',[0.5 0.9],'calc_linmodels',true));
        suptitle('Active/Passive Separability')
        for sessionnum = 1:height(lm_table)
            coefTest(lm_table{sessionnum,3}{1})
        end
        figure('defaultaxesfontsize',18)
        lm_table = plotArrayMap(avg_neuron_eval,struct('map_plot','handelbow_act_eval','clims',[-0.5 0.5],'calc_linmodels',true));
        suptitle('Active/Passive Separability')
        for sessionnum = 1:height(lm_table)
            coefTest(lm_table{sessionnum,3}{1})
        end

