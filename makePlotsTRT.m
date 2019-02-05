%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes plots from results given by analyzeTRT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% Main line figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 1 - Task and classic analysis methods
    % 1a - Monkey using manipulandum (in illustrator)

    % 1b - example movements in two workspaces (with neural firing dots?)
    % First, get example trials
    [~,td_pm_ex] = getTDidx(trial_data,'spaceNum',1,'result','R','rand',1);
    [~,td_dl_ex] = getTDidx(trial_data,'spaceNum',2,'result','R','rand',1);
    % trim to just go from target start to end
    td_ex = trimTD([td_pm_ex td_dl_ex],{'idx_ctHoldTime',0},{'idx_endTime',0});
    % plot the example trials
    figure('defaultaxesfontsize',18)
    plotTRTTrials(td_ex);
    % plot neural firing?
    unit_idx = 1;
    plotSpikesOnHandle(td_ex,struct('unit_idx',unit_idx,'spikespec','b.','spikesize',10));
    % plot of same muscle movement given different Jacobians?

    % Switch to classical PDs?
    % 1c - example directional rasters and tuning curves?

%% Set up plotting variables
    datadir = '/home/raeed/data/project-data/limblab/s1-kinematics/Results/Encoding';
    % filename = {'Han_20171101_TRT_encodingResults_run20180809.mat','Chips_20170915_TRT_encodingResults_run20180809.mat','Lando_20170802_encodingResults_run20180809.mat'};
    files = dir(fullfile(datadir,'*encodingResults_allModels_run20190127.mat'));
    filename = horzcat({files.name});

    monkey_names = {'Chips','Han','Lando'};
    models_to_plot = {'ext','ego','musc','handelbow'};

    % colors for pm, dl conditions
    cond_colors = [...
        231,138,195;...
        166,216,84]/255;

    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;

%% Get pR2 pairwise comparisons for all model pairs and all neurons
    % Go by file and compile
    avg_pR2 = cell(length(monkey_names),size(session_colors,1));
    winners = cell(length(monkey_names),size(session_colors,1));
    session_ctr = zeros(length(monkey_names),1);
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        % classify monkey and session number
        monkey_idx = find(strcmpi(encoderResults.crossEval.monkey{1},monkey_names));
        session_ctr(monkey_idx) = session_ctr(monkey_idx) + 1;

        % get average pR2s by neuron
        avgEval = neuronAverage(encoderResults.crossEval,contains(encoderResults.crossEval.Properties.VariableDescriptions,'meta'));

        avg_pR2{monkey_idx,session_ctr(monkey_idx)} = zeros(height(avgEval),length(models_to_plot));
        for modelnum = 1:length(models_to_plot)
            avg_pR2{monkey_idx,session_ctr(monkey_idx)}(:,modelnum) = avgEval.(sprintf('glm_%s_model_eval',models_to_plot{modelnum}));
        end

        % get comparison
        [winners{monkey_idx,session_ctr(monkey_idx)},model_pairs] = compareEncoderPR2(encoderResults,models_to_plot);
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
            % reset x value
            model_xval = zeros(length(models_to_plot)+1,1);
            for neuronnum = 1:size(avg_pR2{monkeynum,sessionnum})
                yval = model_y + template_y(sessionnum);

                % find highest pR2
                [~,model_sort_idx] = sort(avg_pR2{monkeynum,sessionnum}(neuronnum,:));
                model_order = models_to_plot(model_sort_idx);
                best_model = model_order{end};
                runnerup_model = model_order{end-1};

                % find model pairs involving best model
                best_pairs_idx = any(strcmpi(model_pairs,best_model),2);

                % check decisiveness of victory
                best_model_wins = strcmpi(winners{monkeynum,sessionnum}(best_pairs_idx,neuronnum),best_model);
                if all(best_model_wins)
                    yval = yval(model_sort_idx(end));
                    model_xval(model_sort_idx(end)) = model_xval(model_sort_idx(end)) + 1;
                    xval = model_xval(model_sort_idx(end));
                else
                    yval = yval(end);
                    model_xval(end) = model_xval(end) + 1;
                    xval = model_xval(end);
                end

                % plot dots with darkness coding for value of pR2
                scatter(repmat(xval,size(yval)),yval,[],session_colors(sessionnum,:),'filled')
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
                % find the indices of models to plot for this pair
                model1_idx = find(strcmpi(models_to_plot,model_pairs{pairnum,1}));
                model2_idx = find(strcmpi(models_to_plot,model_pairs{pairnum,2}));

                % scatter filled circles if there's a winner, empty circles if not
                no_winner =  cellfun(@isempty,winners{monkeynum,sessionnum}(pairnum,:));
                scatter(...
                    avg_pR2{monkeynum,sessionnum}(no_winner,model1_idx),...
                    avg_pR2{monkeynum,sessionnum}(no_winner,model2_idx),...
                    [],session_colors(sessionnum,:))
                scatter(...
                    avg_pR2{monkeynum,sessionnum}(~no_winner,model1_idx),...
                    avg_pR2{monkeynum,sessionnum}(~no_winner,model2_idx),...
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

%% Plot pR2 of all monkeys bar plot
    % x coordinate of individual monkey bars
    monk_x = (2:3:((num_monks-1)*3+2))/10;
    % template for within monkey bars separation
    template_x = linspace(-0.5,0.5,length(models_to_plot))/10;
    model_spacing = mode(diff(template_x));

    % make plot
    figure('defaultaxesfontsize',18)
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        avgEval = neuronAverage(encoderResults.crossEval,contains(encoderResults.crossEval.Properties.VariableDescriptions,'meta'));

        % model_idx = find(contains(err{monkeynum}.Properties.VariableNames,models_to_plot));
        % mean_err = mean(err{monkeynum}{:,model_idx});
        % var_err = var(err{monkeynum}{:,model_idx});
        % std_err_err = sqrt(correction*var_err);

        avg_pR2 = zeros(height(avgEval),length(models_to_plot));
        for modelnum = 1:length(models_to_plot)
            xval = monk_x(monkeynum) + template_x(modelnum);
            mean_pR2 = mean(avgEval.(sprintf('glm_%s_model_eval',models_to_plot{modelnum})));
            avg_pR2(:,modelnum) = avgEval.(sprintf('glm_%s_model_eval',models_to_plot{modelnum}));
            bar(xval,mean_pR2,model_spacing,'facecolor',getModelColors(models_to_plot{modelnum}),'edgecolor','none')
            hold on
            % plot([xval xval],[mean_err(modelnum)-std_err_err(modelnum) mean_err(modelnum)+std_err_err(modelnum)],'k','linewidth',3)
        end
        xval = repmat(monk_x(monkeynum)+template_x,length(avg_pR2),1);
        scatter(xval(:),avg_pR2(:),[],'k','filled')
        plot(xval',avg_pR2','-k','linewidth',1)
    end
    set(gca,'tickdir','out','box','off','xtick',monk_x,...
        'xticklabel',filename,'ytick',[0 0.25 0.5],'ticklabelinterpreter','none')
    % axis equal
    ylim([0 0.6])
    % xlim([0 1])
    ylabel('Model pseudo-R^2')

%% Tuning curve shape comparison
    % go by file and compile
        % find correlations between modeled tuning curves and true tuning curve
        tuning_corr = cell(length(monkey_names),size(session_colors,1));
        session_ctr = zeros(length(monkey_names),1);
        for filenum = 1:length(filename)
            % load data
            load(fullfile(datadir,filename{filenum}))

            % classify monkey and session number
            monkey_idx = find(strcmpi(encoderResults.crossEval.monkey{1},monkey_names));
            session_ctr(monkey_idx) = session_ctr(monkey_idx) + 1;

            tuning_corr{monkey_idx,session_ctr(monkey_idx)} = calculateEncoderTuningCorr(...
                encoderResults,struct('model_aliases',{models_to_plot},'neural_signal','S1_FR'));
        end

    % plot by neuron
        figure('defaultaxesfontsize',18)
        % y coordinate of individual monkey bars
        monkey_y = (2:3:((length(monkey_names)-1)*3+2))/10;
        % template for within monkey bars separation
        template_y = linspace(-1,1,length(models_to_plot))/10;
        for monkeynum = 1:length(monkey_names)
            for sessionnum = 1:session_ctr(monkeynum)
                % average for each neuron
                avg_corr = neuronAverage(tuning_corr{monkeynum,sessionnum},...
                    struct('keycols',{{'monkey','date','task','signalID'}},'do_ci',false));
                yval = repmat(monkey_y(monkeynum) + template_y,height(avg_corr),1);
                % add some jitter
                yval = yval+randn(size(yval,1),1)/150;

                % sparsify the lines if needed (not for neurons...)
                doplot = true(length(yval),1);
                cols = contains(avg_corr.Properties.VariableNames,'tuningCorr');
                xvals = avg_corr{:,cols};
                plot(xvals(doplot,:)',yval(doplot,:)','-','linewidth',0.5,'color',ones(1,3)*0.5)
                hold on
                scatter(xvals(:),yval(:),50,session_colors(sessionnum,:),'filled')
            end
        end
        axis ij
        ytickmarks = monkey_y + template_y';
        set(gca,'box','off','tickdir','out',...  'xlim',[0,1],'xtick',0:0.5:1.0,...
            'ytick',ytickmarks(:),'yticklabel',repmat(getModelTitles(models_to_plot),1,length(monkey_names)))
        ylabel(vertcat(monkey_names(:)))
        ylbl = get(gca,'ylabel');
        set(ylbl,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
        title('Modeled tuning curve correlations')
        xlabel('Modeled tuning curve correlations')

    % plot by crossval run
        figure('defaultaxesfontsize',18)
        % y coordinate of individual monkey bars
        monkey_y = (2:3:((length(monkey_names)-1)*3+2))/10;
        % template for within monkey bars separation
        template_y = linspace(-1,1,length(models_to_plot))/10;
        for monkeynum = 1:length(monkey_names)
            for sessionnum = 1:session_ctr(monkeynum)
                % average for each neuron
                avg_corr = neuronAverage(tuning_corr{monkeynum,sessionnum},...
                    struct('keycols',{{'monkey','date','task','crossvalID'}},'do_ci',false));
                yval = repmat(monkey_y(monkeynum) + template_y,height(avg_corr),1);
                % add some jitter
                yval = yval+randn(size(yval,1),1)/150;

                % sparsify the lines
                doplot = rand(length(yval),1)<0.5;
                cols = contains(avg_corr.Properties.VariableNames,'tuningCorr');
                xvals = avg_corr{:,cols};
                plot(xvals(doplot,:)',yval(doplot,:)','-','linewidth',0.5,'color',session_colors(sessionnum,:))
                hold on
                scatter(xvals(:),yval(:),25,ones(1,3)*0.5,'filled')
            end
        end
        axis ij
        ytickmarks = monkey_y + template_y';
        set(gca,'box','off','tickdir','out',...  'xlim',[0,1],'xtick',0:0.5:1.0,...
            'ytick',ytickmarks(:),'yticklabel',repmat(getModelTitles(models_to_plot),1,length(monkey_names)))
        ylabel(vertcat(monkey_names(:)))
        ylbl = get(gca,'ylabel');
        set(ylbl,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
        title('Modeled tuning curve correlations')
        xlabel('Modeled tuning curve correlations')

%% PD shifts over all monkeys
    file_shifts = cell(length(filename),length(models_to_plot)); % shift tables for each model in each file
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        shift_tables = calculatePDShiftTables(encoderResults,[strcat('glm_',models_to_plot,'_model') 'S1_FR']);
        mean_shifts = cell(length(models_to_plot),1);
        for modelnum = 1:length(models_to_plot)+1
            mean_shifts{modelnum} = neuronAverage(shift_tables{modelnum},contains(shift_tables{modelnum}.Properties.VariableDescriptions,'meta'));
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

%% PD shift error dotplot
    % compile error information
    shift_err = cell(length(monkey_names),size(session_colors,1));
    session_ctr = zeros(length(monkey_names),1);
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        % classify monkey and session number
        monkey_idx = find(strcmpi(encoderResults.crossEval.monkey{1},monkey_names));
        session_ctr(monkey_idx) = session_ctr(monkey_idx) + 1;

        shift_err{monkey_idx,session_ctr(monkey_idx)} = calculateEncoderPDShiftErr(encoderResults,struct('model_aliases',{models_to_plot}));
    end

    % plot by neuron
        figure('defaultaxesfontsize',18)
        monkey_y = (2:3:((length(monkey_names)-1)*3+2))/10;
        template_y = linspace(-1,1,length(models_to_plot))/10;
        for monkeynum = 1:length(monkey_names)
            for sessionnum = 1:session_ctr(monkeynum)
                % average for each TUNED neuron
                avg_err = neuronAverage(shift_err{monkeynum,sessionnum},...
                    struct('keycols',{{'monkey','date','task','signalID'}},'do_ci',false));
                yval = repmat(monkey_y(monkeynum) + template_y,height(avg_err),1);
                % add some jitter
                yval = yval+randn(size(yval,1),1)/150;
    
                % sparsify the lines
                doplot = true(length(yval),1);
                cols = contains(avg_err.Properties.VariableNames,'err');
                xvals = avg_err{:,cols};
                plot(xvals(doplot,:)',yval(doplot,:)','-','linewidth',0.5,'color',ones(1,3)*0.5)
                hold on
                scatter(xvals(:),yval(:),50,session_colors(sessionnum,:),'filled')
            end
        end
        axis ij
        ytickmarks = monkey_y + template_y';
        set(gca,'box','off','tickdir','out',...  'xlim',[0,1.5],'xtick',0:0.5:1.5,...
            'ytick',ytickmarks(:),'yticklabel',repmat(getModelTitles(models_to_plot),1,length(monkey_names)))
        ylabel(vertcat(monkey_names(:)))
        ylbl = get(gca,'ylabel');
        set(ylbl,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
        title('PD Shift Model Error')
        xlabel('PD Shift Model Error')

    % plot by crossval run
        figure('defaultaxesfontsize',18)
        monkey_y = (2:3:((length(monkey_names)-1)*3+2))/10;
        template_y = linspace(-1,1,length(models_to_plot))/10;
        for monkeynum = 1:length(monkey_names)
            for sessionnum = 1:session_ctr(monkeynum)
                % average for each TUNED neuron
                avg_err = neuronAverage(shift_err{monkeynum,sessionnum},...
                    struct('keycols',{{'monkey','date','task','crossvalID'}},'do_ci',false));
                yval = repmat(monkey_y(monkeynum) + template_y,height(avg_err),1);
                % add some jitter
                yval = yval+randn(size(yval,1),1)/150;
    
                % sparsify the lines
                doplot = rand(length(yval),1)<0.5;
                cols = contains(avg_err.Properties.VariableNames,'err');
                xvals = avg_err{:,cols};
                plot(xvals(doplot,:)',yval(doplot,:)','-','linewidth',0.5,'color',session_colors(sessionnum,:))
                hold on
                scatter(xvals(:),yval(:),25,ones(1,3)*0.5,'filled')
            end
        end
        axis ij
        ytickmarks = monkey_y + template_y';
        set(gca,'box','off','tickdir','out',...
            'xlim',[0,1.5],'xtick',0:0.5:1.5,...
            'ytick',ytickmarks(:),'yticklabel',repmat(getModelTitles(models_to_plot),1,length(monkey_names)))
        ylabel(vertcat(monkey_names(:)))
        ylbl = get(gca,'ylabel');
        set(ylbl,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
        title('PD Shift Model Error')
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
