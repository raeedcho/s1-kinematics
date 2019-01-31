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

    num_monks = length(filename);
    hyp = cell(num_monks,1);
    p_val = cell(num_monks,1);

    % colors for pm, dl conditions
    cond_colors = [...
        231,138,195;...
        166,216,84]/255;

    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;

%% Get pR2 pairwise comparisons for all model pairs and all neurons
    models_to_plot = {'ext','ego','musc','handelbow'};

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

    % make the plot
    figure('defaultaxesfontsize',18)
    % y coordinate of individual monkey
    monkey_y = (2:3:((length(monkey_names)-1)*3+2))/10;
    % template for within monkey bars separation
    template_y = linspace(-1,1,length(models_to_plot))/10;
    for monkeynum = 1:length(monkey_names)
        % reset x value
        xval = 0;
        for sessionnum = 1:session_ctr(monkeynum)
            for neuronnum = 1:size(avg_pR2{monkeynum,sessionnum})
                yval = monkey_y(monkeynum) + template_y;

                % find highest pR2
                [~,model_sort_idx] = sort(avg_pR2{monkeynum,sessionnum}(neuronnum,:));
                model_order = models_to_plot(model_sort_idx);
                best_model = model_order{end};
                runnerup_model = model_order{end-1};

                % find model pairs involving best model
                best_pairs_idx = any(strcmpi(model_pairs,best_model),2);

                % check decisiveness of victory
                best_model_wins = strcmpi(winners{monkeynum,sessionnum}(best_pairs_idx,neuronnum),best_model);
                if all(best_model_wins) == 0
                    yval = yval(model_sort_idx(end));
                else
                    yval = [];
                end

                % plot dots with darkness coding for value of pR2
                scatter(repmat(xval,size(yval)),yval,[],session_colors(sessionnum,:),'filled')
                % scatter(repmat(xval,size(yval)),yval,[],avg_pR2{monkeynum,sessionnum}(neuronnum,:),'filled')
                hold on

                % increment x value
                xval = xval + 0.1;
            end
        end
    end
    axis ij
    ytickmarks = monkey_y + template_y';
    set(gca,'box','off','tickdir','out',...
        'ytick',ytickmarks(:),'yticklabel',repmat(getModelTitles(models_to_plot),1,length(monkey_names)))
    ylabel(vertcat(monkey_names(:)))
    ylbl = get(gca,'ylabel');
    set(ylbl,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
    title('Pseudo-R^2')
    xlabel('Neuron')

%% Plot pR2 of all monkeys bar plot
    models_to_plot = {'ego','ext','musc','handelbow'};
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
    models_to_plot = {'ego','ext','musc','handelbow'};

    % find correlations between modeled tuning curves and true tuning curve
    tuning_corr = cell(length(monkey_names),size(session_colors,1));
    session_ctr = zeros(length(monkey_names),1);
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        % classify monkey and session number
        monkey_idx = find(strcmpi(encoderResults.crossEval.monkey{1},monkey_names));
        session_ctr(monkey_idx) = session_ctr(monkey_idx) + 1;

        % setup...
        num_bins = encoderResults.params.num_tuning_bins;
        num_neurons = height(encoderResults.tuning_curves{1,1});

        true_tuning_idx = contains(encoderResults.params.model_names,'S1');
        tuning_corr{monkey_idx,session_ctr(monkey_idx)} = zeros(height(encoderResults.tuning_curves{1,1}),length(models_to_plot));

        % arrange tuning curves for each neuron
        for neuron_idx = 1:num_neurons
            tuning_curve_mat = zeros(num_bins*2,length(models_to_plot)+1);
            for spacenum = 1:2
                tuning_curve_mat(num_bins*(spacenum-1)+(1:num_bins),end) = encoderResults.tuning_curves{spacenum,true_tuning_idx}(neuron_idx,:).velCurve';
                for modelnum = 1:length(models_to_plot)
                    tuning_idx = strcmp(encoderResults.params.model_aliases,models_to_plot{modelnum});
                    tuning_curve_mat(num_bins*(spacenum-1)+(1:num_bins),modelnum) = encoderResults.tuning_curves{spacenum,tuning_idx}(neuron_idx,:).velCurve';
                end
            end
            covar_mat = nancov(tuning_curve_mat);
            tuning_corr{monkey_idx,session_ctr(monkey_idx)}(neuron_idx,:) = covar_mat(end,1:end-1)./sqrt(diag(covar_mat(1:end-1,1:end-1))'*covar_mat(end,end));
            % tuning_corr(neuron_idx,:) = covar_mat(end,1:end-1);
        end

        % % bootstrap correlation values for actual tuning curves
        % boot_tuning_corr = zeros(num_neurons,num_boots);
        % boot_tic = tic;
        % for bootnum = 1:num_boots
        %     [~,boot_idx1] = datasample(encoderResults.td_tuning{1},length(encoderResults.td_tuning{1}));
        %     [~,boot_idx2] = datasample(encoderResults.td_tuning{1},length(encoderResults.td_tuning{1}));

        %     % arrange tuning curves for each neuron
        %     for neuron_idx = 1:height(encoderResults.tuning_curves{1,1})
        %         % split into two estimates of tuning curve for each workspace
        %         % then concatenate tuning curves of two workspaces for two overall tuning curve
        %         temp_curves = zeros(num_bins*2,2);
        %         for spacenum = 1:2
        %             tuning_params = struct('out_signals',{{'S1_FR',neuron_idx}},'out_signal_names',1,...
        %                 'num_bins',num_bins,'meta',struct('spaceNum',spacenum));
        %             temp_table1 = getTuningCurves(encoderResults.td_tuning{spacenum}(boot_idx1),tuning_params);
        %             temp_table2 = getTuningCurves(encoderResults.td_tuning{spacenum}(boot_idx2),tuning_params);
        %             temp_curves(num_bins*(spacenum-1)+(1:num_bins),1) = temp_table1.velCurve';
        %             temp_curves(num_bins*(spacenum-1)+(1:num_bins),2) = temp_table2.velCurve';
        %         end

        %         temp_covar = nancov(temp_curves);
        %         boot_tuning_corr(neuron_idx,bootnum) = temp_covar(1,2)/sqrt(prod(nanvar(temp_curves)));
        %     end
        %     % tuning_corr(neuron_idx,:) = covar_mat(end,1:end-1);
        %     fprintf('Bootstrap %d done at time %f\n',bootnum,toc(boot_tic))
        % end
    end

    figure('defaultaxesfontsize',18)
    % y coordinate of individual monkey bars
    monkey_y = (2:3:((length(monkey_names)-1)*3+2))/10;
    % template for within monkey bars separation
    template_y = linspace(-1,1,length(models_to_plot))/10;
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            yval = repmat(monkey_y(monkeynum) + template_y,length(tuning_corr{monkeynum,sessionnum}),1);
            % add some jitter
            yval = yval+randn(size(yval))/150;

            % sparsify the lines
            % doplot = rand(length(yval),1)<0.5;
            doplot = true(length(yval),1);
            plot(tuning_corr{monkeynum,sessionnum}(doplot,:)',yval(doplot,:)','-','linewidth',0.5,'color',ones(1,3)*0.5)
            hold on
            scatter(tuning_corr{monkeynum,sessionnum}(:),yval(:),50,session_colors(sessionnum,:),'filled')
        end
    end
    axis ij
    ytickmarks = monkey_y + template_y';
    set(gca,'box','off','tickdir','out',...
        'xlim',[0,1],'xtick',0:0.5:1.0,...
        'ytick',ytickmarks(:),'yticklabel',repmat(getModelTitles(models_to_plot),1,length(monkey_names)))
    ylabel(vertcat(monkey_names(:)))
    ylbl = get(gca,'ylabel');
    set(ylbl,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
    title('Modeled tuning curve correlations')
    xlabel('Modeled tuning curve correlations')

%% PD shifts over all monkeys
    models_to_plot = {'ext','ego','musc','handelbow'};
    model_titles = getModelTitles(models_to_plot);
    % model_colors = getModelColors(models_to_plot);
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
                    title(sprintf('%s modeled PD shift',model_titles{modelnum}))
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
                %     title(sprintf('%s modeled PD shift',model_titles{modelnum}))
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
                    title(sprintf('%s model',model_titles{modelnum}),'interpreter','none')
                end
            end
    end

%% PD shift error dotplot
    models_to_plot = {'ext','ego','musc','handelbow'};
    monkey_y = (2:3:((length(monkey_names)-1)*3+2))/10;
    template_y = linspace(-1,1,length(models_to_plot))/10;

    % compile error information
    err = cell(length(monkey_names),size(session_colors,1));
    session_ctr = zeros(length(monkey_names),1);
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        % classify monkey and session number
        monkey_idx = find(strcmpi(encoderResults.crossEval.monkey{1},monkey_names));
        session_ctr(monkey_idx) = session_ctr(monkey_idx) + 1;

        err{monkey_idx,session_ctr(monkey_idx)} = calculateEncoderPDShiftErr(encoderResults,struct('model_aliases',{models_to_plot}));
    end

    figure('defaultaxesfontsize',18)
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            yval = repmat(monkey_y(monkeynum) + template_y,length(err{monkeynum,sessionnum}{:,:}),1);
            % add some jitter
            yval = yval+randn(size(yval))/150;

            % sparsify the lines
            doplot = rand(length(yval),1)<0.2;
            plot(err{monkeynum,sessionnum}{:,:}(doplot,:)',yval(doplot,:)','-','linewidth',0.5,'color',ones(1,3)*0.5)
            hold on
            scatter(err{monkeynum,sessionnum}{:,:}(:),yval(:),50,session_colors(sessionnum,:),'filled')
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

%% Plot PD shift error on all monkeys
    num_monks = length(filename);
    correction = 1/100 + 1/4;
    models_to_plot = {'ext','ego','musc','handelbow'};
    % models_to_plot = encoderResults.params.model_aliases;
    model_colors = getModelColors(models_to_plot);
    % models_to_plot = model_aliases;
    % x coordinate of individual monkey bars
    monk_x = (2:3:((num_monks-1)*3+2))/10;
    % template for within monkey bars separation
    template_x = linspace(-1,1,length(models_to_plot))/10;
    model_spacing = mode(diff(template_x));

    % make plot
    figure('defaultaxesfontsize',18)
    for monkeynum = 1:num_monks
        % load data
        load(fullfile(datadir,filename{monkeynum}))
        err{monkeynum} = calculateEncoderPDShiftErr(encoderResults,struct('model_aliases',{models_to_plot}));

        for modelnum = 1:length(models_to_plot)
            mean_err = mean(err{monkeynum}.(models_to_plot{modelnum}));
            var_err = var(err{monkeynum}.(models_to_plot{modelnum}));
            std_err_err = sqrt(correction*var_err);

            xval = monk_x(monkeynum) + template_x(modelnum);
            bar(xval,mean_err,model_spacing,'facecolor',model_colors(modelnum,:),'edgecolor','none')
            hold on
            plot([xval xval],[mean_err-std_err_err mean_err+std_err_err],'k','linewidth',3)
        end
        % xval = repmat(monk_x(monkeynum)+template_x,length(err{monkeynum}{:,:}),1);
        % scatter(xval(:),err{monkeynum}{:,:}(:),[],'k','filled')
        % plot(xval',err{monkeynum}{:,:}','-k','linewidth',1)
    end
    set(gca,'tickdir','out','box','off','xtick',monk_x,...
        'xticklabel',filename,'ytick',[0 0.5],'ticklabelinterpreter','none')
    % axis equal
    ylim([0 0.7])
    % xlim([0 1])
    ylabel('Error of model')

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

%% Decoder stuff
    datadir = '/home/raeed/data/limblab/data-td/FullWS/Results/Decoding';
    filename = {'Han_20160325_RWhold_decodingResults_run20180813.mat','Chips_20151211_RW_decodingResults_run20181029.mat'};
    % datadir = '/home/raeed/data/limblab/data-td/ActPas/Results/Decoding';
    % filename = {'Han_20170203_COactpas_decodingResults_run20181029.mat','Chips_20170913_COactpas_decodingResults_run20181029.mat'};
    num_monks = length(filename);
    barwidth = 0.4;
    markerstyle = 'oo';
    barfig = figure('defaultaxesfontsize',18);
    % scatfig = figure('defaultaxesfontsize',18);
    monk_colors = linspecer(2);
    decoder_colors = [...
        158 33 106;...
        0 159 145]/255;
    for monkeynum = 1:num_monks
        %% load data
        load(fullfile(datadir,filename{monkeynum}))

        %% make decoder performance scatter plots
        x_offset = monkeynum-1;
        % first monkey
        figure(barfig)
        bar([1 3]*0.20+x_offset,mean(decoderResults.hand_decoder_vaf),barwidth,'facecolor',decoder_colors(1,:),'edgecolor','none')
        hold on
        bar([2 4]*0.20+x_offset,mean(decoderResults.neur_decoder_vaf),barwidth,'facecolor',decoder_colors(2,:),'edgecolor','none')
        plot([1 2]*0.20+x_offset,[decoderResults.hand_decoder_vaf(:,1) decoderResults.neur_decoder_vaf(:,1)]','-k','linewidth',1)
        plot([3 4]*0.20+x_offset,[decoderResults.hand_decoder_vaf(:,2) decoderResults.neur_decoder_vaf(:,2)]','-k','linewidth',1)
        plot([1 3 2 4]*0.20+x_offset, [decoderResults.hand_decoder_vaf decoderResults.neur_decoder_vaf]','.k','markersize',30)

        % figure(scatfig)
        % % subplot(1,2,1)
        % scatter(...
        %     decoderResults.hand_decoder_vaf(:,1),...
        %     decoderResults.neur_decoder_vaf(:,1),...
        %     75,...
        %     monk_colors(monkeynum,:),...
        %     markerstyle(monkeynum),...
        %     'filled')
        % hold on
        % plot([0 1],[0 1],'--k','linewidth',2)

        % % subplot(1,2,2)
        % scatter(...
        %     decoderResults.hand_decoder_vaf(:,2),...
        %     decoderResults.neur_decoder_vaf(:,2),...
        %     75,...
        %     monk_colors(monkeynum,:),...
        %     markerstyle(monkeynum),...
        %     'linewidth',2)
        % hold on
        % plot([0 1],[0 1],'--k','linewidth',2)
    end

    % make pretty
    figure(barfig)
    xtick = [(1:4)*0.20 1+(1:4)*0.20];
    xticklabel = repmat({'Hand-only pos','Hand+Neuron pos','Hand vel','Hand+Neuron vel'},1,2);
    set(gca,'box','off','tickdir','out','xtick',xtick,...
        'xticklabel',xticklabel,...
        'xlim',[0 num_monks],'ylim',[0 1])
    ylabel 'Fraction VAF'
    xlabel 'Model'
    title([{'Decoding performance'};filename'],'interpreter','none')

    % figure(scatfig)
    % % subplot(1,2,1)
    % % axis equal
    % % subplot(1,2,2)
    % axis equal
    % set(gca,'box', 'off', 'tickdir', 'out', 'xtick',[0 1],'ytick',[0 1])

    % plot predictions
    % load data
    monkeynum = 2;
    load(fullfile(datadir,filename{monkeynum}))
    figure('defaultaxesfontsize',18)
    last_trial_end = 0;
    elbow_idx = 28:30;
    num_trials = 7;
    for trialnum = 1:num_trials
        true_pos = 100*getSig(decoderResults.td_test(trialnum),{'markers',elbow_idx});
        true_vel = 100*getSig(decoderResults.td_test(trialnum),{'marker_vel',elbow_idx});
        hand_pred_pos = 100*getSig(decoderResults.td_test(trialnum),{'linmodel_hand_decoder',1:3});
        hand_pred_vel = 100*getSig(decoderResults.td_test(trialnum),{'linmodel_hand_decoder',4:6});
        neur_pred_pos = 100*getSig(decoderResults.td_test(trialnum),{'linmodel_neur_decoder',1:3});
        neur_pred_vel = 100*getSig(decoderResults.td_test(trialnum),{'linmodel_neur_decoder',4:6});
        binvec = (1:size(true_pos,1)) + last_trial_end;
        last_trial_end = binvec(end);
        timevec = binvec*decoderResults.td_test(1).bin_size;

        ax(1) = subplot(3,2,1);
        plot(timevec,true_pos(:,3),'-k','linewidth',2)
        hold on
        plot(timevec,hand_pred_pos(:,3),'-','linewidth',2,'color',decoder_colors(1,:))
        plot(timevec,neur_pred_pos(:,3),'-','linewidth',2,'color',decoder_colors(2,:))
        title('Elbow Position')
        ylabel('X-coordinate (cm)');
        ax(2) = subplot(3,2,2);
        plot(timevec,true_vel(:,3),'-k','linewidth',2)
        hold on
        plot(timevec,hand_pred_vel(:,3),'-','linewidth',2,'color',decoder_colors(1,:))
        plot(timevec,neur_pred_vel(:,3),'-','linewidth',2,'color',decoder_colors(2,:))
        title('Elbow Velocity')
        
        ax(3) = subplot(3,2,3);
        plot(timevec,true_pos(:,1),'-k','linewidth',2)
        hold on
        plot(timevec,hand_pred_pos(:,1),'-','linewidth',2,'color',decoder_colors(1,:))
        plot(timevec,neur_pred_pos(:,1),'-','linewidth',2,'color',decoder_colors(2,:))
        ylabel('Y-coordinate (cm)')
        ax(4) = subplot(3,2,4);
        plot(timevec,true_vel(:,1),'-k','linewidth',2)
        hold on
        plot(timevec,hand_pred_vel(:,1),'-','linewidth',2,'color',decoder_colors(1,:))
        plot(timevec,neur_pred_vel(:,1),'-','linewidth',2,'color',decoder_colors(2,:))

        ax(5) = subplot(3,2,5);
        plot(timevec,true_pos(:,2),'-k','linewidth',2)
        hold on
        plot(timevec,hand_pred_pos(:,2),'-','linewidth',2,'color',decoder_colors(1,:))
        plot(timevec,neur_pred_pos(:,2),'-','linewidth',2,'color',decoder_colors(2,:))
        ylabel('Z-coordinate (cm)')
        xlabel('Time (s)')
        ax(6) = subplot(3,2,6);
        plot(timevec,true_vel(:,2),'-k','linewidth',2)
        hold on
        plot(timevec,hand_pred_vel(:,2),'-','linewidth',2,'color',decoder_colors(1,:))
        plot(timevec,neur_pred_vel(:,2),'-','linewidth',2,'color',decoder_colors(2,:))
        xlabel('Time (s)')
    end
    linkaxes(ax,'x')

    trial_ends = cumsum(cat(1,decoderResults.td_test.idx_endTime))';
    goCues = cat(1,decoderResults.td_test.idx_goCueTime)' + [0 trial_ends(1:end-1)];
    trial_ends_time = trial_ends*decoderResults.td_test(1).bin_size;
    goCues_time = goCues*decoderResults.td_test(1).bin_size;
    for plotnum = 1:length(ax)
        subplot(3,2,plotnum)
        ylims = get(gca,'ylim')';
        targ_y = [0.9 0.1] * ylims;
        plot(repmat(trial_ends_time,2,1),repmat(ylims,1,numel(trial_ends_time)),'--k','linewidth',3)
        scatter(goCues_time(:)',repmat(targ_y,1,numel(goCues_time)),100,'r','s','filled')
        set(gca,'box','off','tickdir','out')
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
