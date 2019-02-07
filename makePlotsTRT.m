%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes plots from results saved by calculateEncoders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up plotting variables
    datadir = '/home/raeed/data/project-data/limblab/s1-kinematics/Results/Encoding';
    % filename = {'Han_20171101_TRT_encodingResults_run20180809.mat','Chips_20170915_TRT_encodingResults_run20180809.mat','Lando_20170802_encodingResults_run20180809.mat'};
    files = dir(fullfile(datadir,'*encodingResults_allModels_run20190206.mat'));
    filename = horzcat({files.name});

    monkey_names = {'Chips','Han','Lando'};
    models_to_plot = {'ego','ext','musc','handelbow'};

    % colors for pm, dl conditions
    cond_colors = [...
        231,138,195;...
        166,216,84]/255;

    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;

%% Compile information over all files
    [model_eval,tuning_corr,shift_err] = deal(cell(length(monkey_names),size(session_colors,1)));
    session_ctr = zeros(length(monkey_names),1);
    fileclock = tic;
    fprintf('Started loading files...\n')
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        % classify monkey and session number
        monkey_idx = find(strcmpi(encoderResults.crossEval.monkey{1},monkey_names));
        session_ctr(monkey_idx) = session_ctr(monkey_idx) + 1;

        % We already have evaluation table in crossEval... just extract the models we want
        model_eval{monkey_idx,session_ctr(monkey_idx)} = encoderResults.crossEval(:,contains(encoderResults.crossEval.Properties.VariableDescriptions,'meta'));
        model_eval_cell = cell(1,length(models_to_plot));
        for modelnum = 1:length(models_to_plot)
            model_eval_cell{modelnum} = table(encoderResults.crossEval.(sprintf('glm_%s_model_eval',models_to_plot{modelnum})),...
                'VariableNames',strcat(models_to_plot(modelnum),'_eval'));
            model_eval_cell{modelnum}.Properties.VariableDescriptions = {'linear'};
        end

        model_eval{monkey_idx,session_ctr(monkey_idx)} = horzcat(...
            model_eval{monkey_idx,session_ctr(monkey_idx)},...
            model_eval_cell{:});

        % Get tuning curve correlation table
        tuning_corr{monkey_idx,session_ctr(monkey_idx)} = calculateEncoderTuningCorr(...
            encoderResults,struct('model_aliases',{models_to_plot},'neural_signal','S1_FR'));

        % Get PD shift error table
        shift_err{monkey_idx,session_ctr(monkey_idx)} = calculateEncoderPDShiftErr(...
            encoderResults,struct('model_aliases',{models_to_plot}));

        % output a counter
        fprintf('Processed file %d of %d at time %f\n',filenum,length(filename),toc(fileclock))
    end

%% Get pR2 pairwise comparisons for all model pairs and all neurons
    % find winners of pR2
        pr2_winners = cell(length(monkey_names),size(session_colors,1));
        for monkeynum = 1:length(monkey_names)
            for sessionnum = 1:session_ctr(monkeynum)
                [pr2_winners{monkeynum,sessionnum},model_pairs] = compareEncoderMetrics(...
                        model_eval{monkeynum,sessionnum},struct(...
                            'models',{models_to_plot},...
                            'postfix','_eval'));
            end
        end

    % make the pairwise comparison scatter plot
        figure
        for monkeynum = 1:length(monkey_names)
            for pairnum = 1:size(model_pairs,1)
                % set subplot
                subplot(length(monkey_names),size(model_pairs,1),...
                    (monkeynum-1)*size(model_pairs,1)+pairnum)
                plot([-1 1],[-1 1],'k--','linewidth',0.5)
                hold on
                plot([0 0],[-1 1],'k-','linewidth',0.5)
                plot([-1 1],[0 0],'k-','linewidth',0.5)
                for sessionnum = 1:session_ctr(monkeynum)
                    avg_pR2 = neuronAverage(model_eval{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
                    % scatter filled circles if there's a winner, empty circles if not
                    no_winner =  cellfun(@isempty,pr2_winners{monkeynum,sessionnum}(pairnum,:));
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
                set(gca,'box','off','tickdir','out',...
                    'xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
                axis square
                if monkeynum ~= length(monkey_names) || pairnum ~= 1
                    set(gca,'box','off','tickdir','out',...
                        'xtick',[],'ytick',[])
                end
                xlabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,1})))
                ylabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,2})))
            end
        end

    % make the winner dot plot
        figure('defaultaxesfontsize',18)
        % y coordinate of individual models
        model_y = (2:3:(length(models_to_plot)*3+2))/10;
        for monkeynum = 1:length(monkey_names)
            % different subplots for different monkeys
            subplot(1,3,monkeynum)

            % template for within model session differences
            template_y = linspace(-0.3,0.3,session_ctr(monkeynum))/10;

            % make dotplot
            for sessionnum = 1:session_ctr(monkeynum)
                % get avg pR2
                avg_pR2 = neuronAverage(model_eval{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
                model_cols = endsWith(avg_pR2.Properties.VariableNames,'_eval');
                model_colnames = strrep(avg_pR2.Properties.VariableNames(model_cols),'_eval','');

                % reset x value
                model_xval = zeros(length(models_to_plot)+1,1);
                for neuronnum = 1:height(avg_pR2)
                    yval = model_y + template_y(sessionnum);

                    % find highest pR2
                    [~,model_sort_idx] = sort(avg_pR2{neuronnum,model_cols});
                    model_order = model_colnames(model_sort_idx);
                    best_model = model_order{end};
                    runnerup_model = model_order{end-1};

                    % check decisiveness of victory
                    best_model_wins = strcmpi(pr2_winners{monkeynum,sessionnum}(:,neuronnum),best_model);
                    runnerup_wins = strcmpi(pr2_winners{monkeynum,sessionnum}(:,neuronnum),runnerup_model);
                    if sum(best_model_wins) == length(models_to_plot)-1
                        yval_plot = yval(model_sort_idx(end));
                        model_xval(model_sort_idx(end)) = model_xval(model_sort_idx(end)) + 1;
                        xval_plot = model_xval(model_sort_idx(end));
                        scatter(repmat(xval_plot,size(yval_plot)),yval_plot,[],session_colors(sessionnum,:),'filled')
                    elseif sum(best_model_wins) == length(models_to_plot)-2
                        yval_plot = yval(model_sort_idx(end));
                        model_xval(model_sort_idx(end)) = model_xval(model_sort_idx(end)) + 1;
                        xval_plot = model_xval(model_sort_idx(end));
                        scatter(repmat(xval_plot,size(yval_plot)),yval_plot,[],session_colors(sessionnum,:),'<','filled')

                        if sum(runnerup_wins) == length(models_to_plot)-2
                            yval_plot = yval(model_sort_idx(end-1));
                            model_xval(model_sort_idx(end-1)) = model_xval(model_sort_idx(end-1)) + 1;
                            xval_plot = model_xval(model_sort_idx(end-1));
                            scatter(repmat(xval_plot,size(yval_plot)),yval_plot,[],session_colors(sessionnum,:),'<','filled')
                        end
                    else
                        yval_plot = yval(end);
                        model_xval(end) = model_xval(end) + 1;
                        xval_plot = model_xval(end);
                        scatter(repmat(xval_plot,size(yval_plot)),yval_plot,[],session_colors(sessionnum,:),'filled')
                    end

                    % scatter(repmat(model_xval,size(yval)),yval,[],avg_pR2{monkeynum,sessionnum}(neuronnum,:),'filled')
                    hold on

                    % increment x value
                end
            end
            title(monkey_names{monkeynum})
            axis ij
            if monkeynum == 1
                set(gca,'box','off','tickdir','out',...
                    'ylim',[model_y(1)+template_y(1)-0.1 model_y(end)+template_y(end)+0.1],...
                    'ytick',model_y,...
                    'xlim',[0 30],...
                    'xtick',0:10:30,...
                    'yticklabel',[getModelTitles(models_to_plot);{'None'}])
            else
                set(gca,'box','off','tickdir','out',...
                    'ylim',[model_y(1)+template_y(1)-0.1 model_y(end)+template_y(end)+0.1],...
                    'ytick',model_y,...
                    'xlim',[0 30],...
                    'xtick',0:10:30,...
                    'yticklabel',{})
            end
            xlabel('Number of neurons')
        end

    % plot by neuron
    figure('defaultaxesfontsize',18)
    plotModelMetric(model_eval,struct(...
        'monkey_names',{monkey_names},...
        'session_ctr',session_ctr,...
        'session_colors',session_colors,...
        'models_to_plot',{models_to_plot},...
        'postfix','_eval',...
        'marginal_col','crossvalID',...
        'line_sparsity',0));
    xlabel('Model Pseudo-R^2')

    % plot by crossval run
    figure('defaultaxesfontsize',18)
    plotModelMetric(model_eval,struct(...
        'monkey_names',{monkey_names},...
        'session_ctr',session_ctr,...
        'session_colors',session_colors,...
        'models_to_plot',{models_to_plot},...
        'postfix','_eval',...
        'marginal_col','signalID',...
        'line_sparsity',0));
    xlabel('Model Pseudo-R^2')

%% Tuning curve shape comparison
    % plot by neuron
    figure('defaultaxesfontsize',18)
    plotModelMetric(tuning_corr,struct(...
        'monkey_names',{monkey_names},...
        'session_ctr',session_ctr,...
        'session_colors',session_colors,...
        'models_to_plot',{models_to_plot},...
        'postfix','_tuningCorr',...
        'marginal_col','crossvalID',...
        'line_sparsity',0));
    xlabel('Modeled tuning curve correlations')

    % plot by crossval run
    figure('defaultaxesfontsize',18)
    plotModelMetric(tuning_corr,struct(...
        'monkey_names',{monkey_names},...
        'session_ctr',session_ctr,...
        'session_colors',session_colors,...
        'models_to_plot',{models_to_plot},...
        'postfix','_tuningCorr',...
        'marginal_col','signalID',...
        'line_sparsity',0));
    xlabel('Modeled tuning curve correlations')

%% PD shifts over all monkeys
    file_shifts = cell(length(filename),length(models_to_plot)); % shift tables for each model in each file
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        shift_tables = calculatePDShiftTables(encoderResults,[strcat('glm_',models_to_plot,'_model') 'S1_FR']);
        mean_shifts = cell(length(models_to_plot),1);
        for modelnum = 1:length(models_to_plot)+1
            mean_shifts{modelnum} = neuronAverage(shift_tables{modelnum},struct(...
                'keycols',{{'monkey','date','task','signalID'}}));
            [~,file_shifts{filenum,modelnum}] = getNTidx(mean_shifts{modelnum},'signalID',encoderResults.tunedNeurons);
        end
    end

    allFileShifts_real = vertcat(file_shifts{:,end});

    % Make histograms and scatters
    % hists = figure('defaultaxesfontsize',18);
    total_hists = figure('defaultaxesfontsize',18);
    scatters = figure('defaultaxesfontsize',18);
    for monkeynum = 1:length(monkey_names)
        % get monkey specific session dates
        [~,monkey_shifts_real] = getNTidx(allFileShifts_real,'monkey',monkey_names{monkeynum});
        session_dates = unique(monkey_shifts_real.date);

        % make histogram combining sessions
            % first figure out ylim
            % ylim_high = 10*floor(height(monkey_shifts_real)/20);
            ylim_high = 40;
            % actual PD shift histogram
            figure(total_hists)
            subplot(length(monkey_names),length(models_to_plot)+1,(monkeynum-1)*(length(models_to_plot)+1)+1)
            h = histogram(gca,monkey_shifts_real.velPD*180/pi,'BinWidth',10,'DisplayStyle','bar');
            set(h,'facecolor','k','edgecolor','none')
            subplot(length(monkey_names),length(models_to_plot)+1,(monkeynum-1)*(length(models_to_plot)+1)+1)
            set(gca,'box','off','tickdir','out','xlim',[-180 180],'xtick',[-180 0 180],'ylim',[0 ylim_high],'ytick',[0 ylim_high/2 ylim_high],'view',[-90 90])
            ylabel(monkey_names{monkeynum})
            if monkeynum == 1
                title('Actual PD Shift')
            end
            for modelnum = 1:length(models_to_plot)
                allFileShifts_model = vertcat(file_shifts{:,modelnum});
                [~,monkey_shifts_model] = getNTidx(allFileShifts_model,'monkey',monkey_names{monkeynum});

                % modeled PD shift histogram
                subplot(length(monkey_names),length(models_to_plot)+1,(monkeynum-1)*(length(models_to_plot)+1)+modelnum+1)
                h = histogram(gca,monkey_shifts_model.velPD*180/pi,'BinWidth',10,'DisplayStyle','bar');
                set(h,'facecolor','k','edgecolor','none')
                subplot(length(monkey_names),length(models_to_plot)+1,(monkeynum-1)*(length(models_to_plot)+1)+modelnum+1)
                set(gca,'box','off','tickdir','out','xlim',[-180 180],'xtick',[-180 0 180],'ylim',[0 ylim_high],'ytick',[0 ylim_high/2 ylim_high],'view',[-90 90])
                if monkeynum == 1
                    title(sprintf('%s modeled PD shift',getModelTitles(models_to_plot{modelnum})))
                end
            end

        % make plots separating sessions
            for sessionnum = 1:length(session_dates)
                % get real shifts for this session
                [~,session_shifts_real] = getNTidx(allFileShifts_real,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});

                % first the real PD shift histogram
                % figure(hists)
                % subplot(length(monkey_names),length(models_to_plot)+1,(monkeynum-1)*(length(models_to_plot)+1)+1)
                % h = histogram(gca,session_shifts_real.velPD*180/pi,'BinWidth',10,'DisplayStyle','bar');
                % set(h,'facecolor',session_colors(sessionnum,:),'edgecolor','none')
                % hold on

                % now the models
                for modelnum = 1:length(models_to_plot)
                    % get the modeled shifts for this session
                    allFileShifts_model = vertcat(file_shifts{:,modelnum});
                    [~,session_shifts_model] = getNTidx(allFileShifts_model,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});

                    % % modeled PD shift histogram
                    % figure(hists)
                    % subplot(length(monkey_names),length(models_to_plot)+1,(monkeynum-1)*(length(models_to_plot)+1)+modelnum+1)
                    % h = histogram(gca,session_shifts_model.velPD*180/pi,'BinWidth',10,'DisplayStyle','bar');
                    % set(h,'facecolor',session_colors(sessionnum,:),'edgecolor','none')
                    % hold on

                    % scatter plots
                    figure(scatters)
                    subplot(length(monkey_names),length(models_to_plot),(monkeynum-1)*(length(models_to_plot))+modelnum)
                    % hsh = scatterhist(180/pi*session_shifts_model.velPD,180/pi*session_shifts_real.velPD,...
                    %     'markersize',50,'group',allFileShifts_real.monkey,'location','NorthWest',...
                    %     'direction','out','plotgroup','off','color',monkey_colors,'marker','...',...
                    %     'nbins',[15 15],'style','stairs');
                    % hsh(3).Children.EdgeColor = [0 0 0];
                    % hsh(2).Children.EdgeColor = model_colors(modelnum,:);
                    scatter(180/pi*session_shifts_model.velPD,180/pi*session_shifts_real.velPD,50,session_colors(sessionnum,:),'filled')
                    hold on
                end
            end
            % make axes pretty
            % figure(hists)
            % subplot(length(monkey_names),length(models_to_plot)+1,(monkeynum-1)*(length(models_to_plot)+1)+1)
            % set(gca,'box','off','tickdir','out','xlim',[-180 180],'xtick',[-180 0 180],'ylim',[0 15],'ytick',[0 15 30],'view',[-90 90])
            % ylabel(monkey_names{monkeynum})
            % if monkeynum == 1
            %     title('Actual PD Shift')
            % end
            for modelnum = 1:length(models_to_plot)
                % histograms
                % figure(hists)
                % subplot(length(monkey_names),length(models_to_plot)+1,(monkeynum-1)*(length(models_to_plot)+1)+modelnum+1)
                % set(gca,'box','off','tickdir','out','xlim',[-180 180],'xtick',[-180 0 180],'ylim',[0 15],'ytick',[0 15 30],'view',[-90 90])
                % if monkeynum == 1
                %     title(sprintf('%s modeled PD shift',getModelTitles(models_to_plot{modelnum})))
                % end

                % scatter plots
                figure(scatters)
                subplot(length(monkey_names),length(models_to_plot),(monkeynum-1)*length(models_to_plot)+modelnum)
                plot([-180 180],[0 0],'-k','linewidth',2)
                plot([0 0],[-180 180],'-k','linewidth',2)
                plot([-180 180],[-180 180],'--k','linewidth',2)
                axis equal
                set(gca,'box','off','tickdir','out','xtick',[-180 180],'ytick',[-180 180],'xlim',[-180 180],'ylim',[-180 180])
                % labels
                if monkeynum == length(monkey_names)
                    xlabel 'Modeled PD Shift'
                end
                if modelnum == 1
                    ylabel({monkey_names{monkeynum};'Actual PD Shift'})
                    ylbl = get(gca,'ylabel');
                    set(ylbl,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right')
                end
                if monkeynum == 1
                    title(sprintf('%s model',getModelTitles(models_to_plot{modelnum})),'interpreter','none')
                end
            end
    end

%% PD shift error dotplots
    % plot by neuron
    figure('defaultaxesfontsize',18)
    plotModelMetric(shift_err,struct(...
        'monkey_names',{monkey_names},...
        'session_ctr',session_ctr,...
        'session_colors',session_colors,...
        'models_to_plot',{models_to_plot},...
        'postfix','_err',...
        'marginal_col','crossvalID',...
        'line_sparsity',0));
    xlabel('PD Shift Model Error')

    % plot by crossval run
    figure('defaultaxesfontsize',18)
    plotModelMetric(shift_err,struct(...
        'monkey_names',{monkey_names},...
        'session_ctr',session_ctr,...
        'session_colors',session_colors,...
        'models_to_plot',{models_to_plot},...
        'postfix','_err',...
        'marginal_col','signalID',...
        'line_sparsity',0));
    xlabel('PD Shift Model Error')

%% Get example tuning curves for all models
    for monkeynum = 1%:num_monks
        clear encoderResults

        % load data
        load(fullfile(datadir,filename{monkeynum}))

        %% Plot out tuning curves
            % compare PM and DL tuning for each model
            for modelnum = 1:num_models
                figure('defaultaxesfontsize',18)
                % figure
                compareTuning(encoderResults.tuning_curves(:,modelnum),encoderResults.pdTables(:,modelnum),struct('which_units',find(encoderResults.isTuned),'cond_colors',cond_colors))
                % compareTuning(encoderResults.tuning_curves(:,modelnum),encoderResults.pdTables(:,modelnum),struct('which_units',find(encoderResults.isTuned),'cond_colors',cond_colors,'maxFR',1))
                title(encoderResults.params.model_names{modelnum},'interpreter','none')
            end
    end

%% Plot out example firing rates
    % load data
    filenum = 4;
    load(fullfile(datadir,filename{filenum}))

    % set plotting params
    num_trials = 10;
    trials_to_plot = randperm(length(encoderResults.td_tuning{2}),num_trials);

    % just for fun...
    model_random = models_to_plot(randperm(length(models_to_plot)));

    for neuronnum = 1:size(encoderResults.td_tuning{2}(1).S1_spikes,2)
        h = figure('defaultaxesfontsize',18);
        ax = zeros(length(models_to_plot),2);
        for spacenum = 1:2
            ax(spacenum) = subplot(2,1,spacenum);
            plotExampleFR(encoderResults.td_tuning{spacenum},...
                struct('neuron_idx',neuronnum,'models',{models_to_plot},'trial_idx',trials_to_plot))
            title(sprintf('Neuron %d, spacenum %d',neuronnum,spacenum))
        end
        linkaxes(ax(:),'y')
        waitfor(h)
    end

%% Extra stuff/in progress...
    %% Plot handle positions
        figure('defaultaxesfontsize',18)
        pos_dl = cat(1,results.td_test{2}.pos);
        plot(pos_dl(:,1),pos_dl(:,2),'r')
        hold on
        pos_pm = cat(1,results.td_test{1}.pos);
        plot(pos_pm(:,1),pos_pm(:,2),'b')
        axis equal
        
        % clean up
        clearvars pos_*
