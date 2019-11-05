%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes plots from results saved by calculateEncoders
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up plotting variables
    if ispc
        dataroot = '';
    else
        dataroot = '/data/raeed';
    end

    datadir = fullfile(dataroot,'project-data','limblab','s1-kinematics','Results','Encoding');
    % filename = {'Han_20171101_TRT_encodingResults_run20180809.mat','Chips_20170915_TRT_encodingResults_run20180809.mat','Lando_20170802_encodingResults_run20180809.mat'};
    % files = dir(fullfile(datadir,'*encodingResults_allModels_run20190206.mat'));
    % files = dir(fullfile(datadir,'*encodingResults_allModels_run20190808.mat'));
    % files = dir(fullfile(datadir,'*encodingResults_opensimControls_run20190311.mat'));
    % files = dir(fullfile(datadir,'*encodingResults_*_run20191001.mat'));
    files = dir(fullfile(datadir,'*encodingResults_*_run20191031.mat'));
    filename = horzcat({files.name});
    
    % for figure saving
    figdir = '/home/raeed/Wiki/Projects/limblab/s1-kinematics/figures/Encoding';
    run_date = char(datetime('today','format','yyyyMMdd'));

    monkey_names = {'Chips','Han','Lando'};
    models_to_plot = {'ext','handelbow','extforce'};
    hand_models = {'ext','ego'};
    arm_models = {'joint','musc','handelbow'};
    included_models = union(union(models_to_plot,hand_models),arm_models);
    not_plot_models = setdiff(included_models,models_to_plot);
    % models_to_plot = {'elbow','ext','musc','handelbow'};
    % models_to_plot = {'elbow','ext','extforce','handelbow'};
    % models_to_plot = {'ext','extforce','handelbow','handelbowforce'};
    % models_to_plot = {'ego','ext','joint','musc','handelbow'};
    % models_to_plot = {'joint','musc','opensim_handelbow','handelbow'};

    % colors for pm, dl conditions
    cond_colors = [...
        231,138,195;...
        166,216,84]/255;

    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;

%% Compile information over all files
    [model_eval,model_tuning,tuning_corr,shift_vaf,tuned_neurons] = deal(cell(length(monkey_names),size(session_colors,1)));
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
        model_eval_cell = cell(1,length(included_models));
        space_eval_cell = cell(2,length(included_models));
        for modelnum = 1:length(included_models)
            model_eval_cell{modelnum} = table(encoderResults.crossEval.(sprintf('glm_%s_model_eval',included_models{modelnum})),...
                'VariableNames',strcat(included_models(modelnum),'_eval'));
            model_eval_cell{modelnum}.Properties.VariableDescriptions = {'linear'};
            for spacenum = 1:2
                space_eval_cell{spacenum,modelnum} = table(encoderResults.crossEval.(sprintf('glm_%s_model_space%d_eval',included_models{modelnum},spacenum)),...
                    'VariableNames',{sprintf('%s_space%d_eval',included_models{modelnum},spacenum)});
                space_eval_cell{spacenum,modelnum}.Properties.VariableDescriptions = {'linear'};
            end
        end
        model_eval{monkey_idx,session_ctr(monkey_idx)} = horzcat(...
            model_eval{monkey_idx,session_ctr(monkey_idx)},...
            model_eval_cell{:},...
            space_eval_cell{:});

        % We already have tuning table in crossTuning... just extract the models we want
        model_tuning{monkey_idx,session_ctr(monkey_idx)} = encoderResults.crossTuning(:,...
            contains(encoderResults.crossTuning.Properties.VariableDescriptions,'meta') |...
            strcmpi(encoderResults.crossTuning.Properties.VariableNames,'bins'));
        model_tuning_cell = cell(1,length(included_models)+1);
        for modelnum = 1:length(included_models)
            model_tuning_cell{modelnum} = table(...
                encoderResults.crossTuning.(sprintf('glm_%s_model_velCurve',included_models{modelnum})),...
                encoderResults.crossTuning.(sprintf('glm_%s_model_velPD',included_models{modelnum})),...
                'VariableNames',strcat(included_models(modelnum),{'_velCurve','_velPD'}));
            model_tuning_cell{modelnum}.Properties.VariableDescriptions = {'linear','circular'};
        end
        model_tuning_cell{end} = table(...
            encoderResults.crossTuning.('S1_FR_velCurve'),...
            encoderResults.crossTuning.('S1_FR_velPD'),...
            'VariableNames',strcat('S1_FR',{'_velCurve','_velPD'}));
        model_tuning_cell{end}.Properties.VariableDescriptions = {'linear','circular'};
        % put it together
        model_tuning{monkey_idx,session_ctr(monkey_idx)} = horzcat(...
            model_tuning{monkey_idx,session_ctr(monkey_idx)},...
            model_tuning_cell{:});

        % Get tuning curve correlation table
        tuning_corr{monkey_idx,session_ctr(monkey_idx)} = calculateEncoderTuningCorr(...
            encoderResults,struct('model_aliases',{included_models},'neural_signal','S1_FR'));

        % Get PD shift error table
        shift_vaf{monkey_idx,session_ctr(monkey_idx)} = calculateEncoderPDShiftVAF(...
            encoderResults,struct('model_aliases',{included_models}));

        % get tuned neurons
        tuned_neurons{monkey_idx,session_ctr(monkey_idx)} = encoderResults.tunedNeurons;

        % output a counter
        fprintf('Processed file %d of %d at time %f\n',filenum,length(filename),toc(fileclock))
    end

%% Get pR2 pairwise comparisons for model pairs and all neurons
    % find winners of pR2
    pr2_winners = cell(length(monkey_names),size(session_colors,1));
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            [pr2_winners{monkeynum,sessionnum},model_pairs] = compareEncoderMetrics(...
                    model_eval{monkeynum,sessionnum},struct(...
                        'bonferroni_correction',6,...
                        'models',{models_to_plot},...
                        'model_pairs',{{'ext','handelbow';'extforce','handelbow'}},...
                        'postfix','_eval'));
        end
    end

    % figure out how many neurons the hand models could beat either of the whole-arm models
    % hand_neuron_counter = 0;
    % rowchecks = contains(model_pairs,hand_models);
    % rowchecks = xor(rowchecks(:,1),rowchecks(:,2));
    % for monkeynum = 1:length(monkey_names)
    %     for sessionnum = 1:session_ctr(monkeynum)
    %         hand_neuron_counter = hand_neuron_counter + ...
    %             sum(any(contains(pr2_winners{monkeynum,sessionnum}(rowchecks,:),{'ext','ego'})),2)
    %     end
    % end
    all_pr2_winners = horzcat(pr2_winners{:});
    ext_winners = sum(strcmpi(all_pr2_winners,'ext'),2);
    handelbow_winners = sum(strcmpi(all_pr2_winners,'handelbow'),2);
    extforce_winners = sum(strcmpi(all_pr2_winners,'extforce'),2);

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
            if monkeynum ~= 1 || pairnum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,1})))
            ylabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,2})))
        end
    end
    suptitle('Pseudo-R^2 pairwise comparisons')
    % saveas(gcf,fullfile(figdir,sprintf('pr2_pairwise_run%s.pdf',run_date)))

    % show scatter plot for hand/elbow pR2 within condition vs against condition
    for modelnum = 1:length(models_to_plot)
        figure
        for monkeynum = 1:length(monkey_names)
            for spacenum = 1:2
                % set subplot
                subplot(length(monkey_names),2,(monkeynum-1)*2+spacenum)

                % plot lines
                plot([-1 1],[-1 1],'k--','linewidth',0.5)
                hold on
                plot([0 0],[-1 1],'k-','linewidth',0.5)
                plot([-1 1],[0 0],'k-','linewidth',0.5)

                % plot out each session
                for sessionnum = 1:session_ctr(monkeynum)
                    avg_pR2 = neuronAverage(model_eval{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
                    scatter(...
                        avg_pR2.(sprintf('%s_space%d_eval',models_to_plot{modelnum},spacenum)),...
                        avg_pR2.(sprintf('%s_eval',models_to_plot{modelnum})),...
                        [],session_colors(sessionnum,:),'filled')
                end
                % make axes pretty
                set(gca,'box','off','tickdir','out',...
                    'xlim',[-0.1 0.6],'ylim',[-0.1 0.6])
                axis square
                if monkeynum ~= 1 || spacenum ~= 1
                    set(gca,'box','off','tickdir','out',...
                        'xtick',[],'ytick',[])
                end
                xlabel(sprintf('%s workspace %d pR2',getModelTitles(models_to_plot{modelnum}),spacenum))
                ylabel(sprintf('%s full pR2',getModelTitles(models_to_plot{modelnum})))
            end
        end
        suptitle('Full pR^2 vs within condition pR^2')
        saveas(gcf,fullfile(figdir,sprintf('%s_pr2_full_v_within_run%s.pdf',models_to_plot{modelnum},run_date)))
    end

    % make the winner dot plot
        % figure('defaultaxesfontsize',18)
        % plotAnalysisWinners(model_eval,pr2_winners,struct(...
        %     'monkey_names',{monkey_names},...
        %     'session_ctr',session_ctr,...
        %     'session_colors',session_colors,...
        %     'models_to_plot',{models_to_plot},...
        %     'postfix','_eval'))
        % suptitle('Pseudo-R^2 Winners')
        % saveas(gcf,fullfile(figdir,sprintf('pr2_winners_run%s.pdf',run_date)))

    % plot array map of ext model eval
        % plot all evals together
        model_eval_table = vertcat(model_eval{:});
        model_eval_table = horzcat(model_eval_table,...
            table(model_eval_table.handelbow_eval-model_eval_table.ext_eval,...
                'VariableNames',{'he_ext_eval_diff'}));
        model_eval_table = horzcat(model_eval_table,...
            table(model_eval_table.ext_eval-model_eval_table.elbow_eval,...
                'VariableNames',{'ext_elbow_eval_diff'}));
        % average over all crossvals
        model_eval_table = neuronAverage(model_eval_table,struct('keycols',{{'monkey','date','task','signalID'}},'do_ci',false));
        % extract and add chan
        model_eval_table = horzcat(model_eval_table,...
            table(model_eval_table.signalID(:,1),'VariableNames',{'chan'}));
        % extract and add color
        % for modelnum = 1:length(models_to_plot)
        %     model_eval_table = horzcat(model_eval_table,...
        %         table(model_eval_table.(sprintf('%s_eval',models_to_plot{modelnum})),...
        %             'VariableNames',{sprintf('%s_eval_color',models_to_plot{modelnum})}));
        % end
        % for modelnum = 1:length(models_to_plot)
        %     plotArrayMap(model_eval_table,struct('map_plot',sprintf('%s_eval',models_to_plot{modelnum})))
        %     suptitle(sprintf('%s pR2',getModelTitles(models_to_plot{modelnum})))
        % end
        coarseness = 1;
        figure('defaultaxesfontsize',18)
        lm_table = plotArrayMap(model_eval_table,struct(...
            'map_plot','ext_elbow_eval_diff',...
            'calc_linmodels',true,...
            'coarseness',coarseness));
        suptitle(sprintf('Extrinsic pR2 (coarseness factor: %d)',coarseness))
        % plotArrayMap(model_eval_table,struct('map_plot','ext_elbow_eval_diff','clims',[-0.15 0.15]))
        % suptitle('Hand/Elbow pR2 - Extrinsic pR2')
        for sessionnum = 1:height(lm_table)
            coefTest(lm_table{sessionnum,3}{1})
        end

%% Tuning curve shape comparison
    % find winners of tuning corr
    tuning_corr_winners = cell(length(monkey_names),size(session_colors,1));
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            [tuning_corr_winners{monkeynum,sessionnum},model_pairs] = compareEncoderMetrics(...
                    tuning_corr{monkeynum,sessionnum},struct(...
                        'bonferroni_correction',6,...
                        'models',{models_to_plot},...
                        'model_pairs',{{'ext','handelbow';'extforce','handelbow'}},...
                        'postfix','_tuningCorr'));
        end
    end

    % figure out how many neurons the hand-based models could beat either of the whole-arm models
    % hand_neuron_counter = 0;
    % rowchecks = contains(model_pairs,{'ext','ego'});
    % rowchecks = xor(rowchecks(:,1),rowchecks(:,2));
    % for monkeynum = 1:length(monkey_names)
    %     for sessionnum = 1:session_ctr(monkeynum)
    %         hand_neuron_counter = hand_neuron_counter + ...
    %             sum(any(contains(tuning_corr_winners{monkeynum,sessionnum}(rowchecks,:),{'ext','ego'})),2)
    %     end
    % end
    all_tuning_corr_winners = horzcat(tuning_corr_winners{:});
    ext_winners = sum(strcmpi(all_tuning_corr_winners,'ext'),2);
    handelbow_winners = sum(strcmpi(all_tuning_corr_winners,'handelbow'),2);
    extforce_winners = sum(strcmpi(all_tuning_corr_winners,'extforce'),2);

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
                avg_corr = neuronAverage(tuning_corr{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
                % scatter filled circles if there's a winner, empty circles if not
                no_winner =  cellfun(@isempty,tuning_corr_winners{monkeynum,sessionnum}(pairnum,:));
                scatter(...
                    avg_corr.(strcat(model_pairs{pairnum,1},'_tuningCorr'))(no_winner),...
                    avg_corr.(strcat(model_pairs{pairnum,2},'_tuningCorr'))(no_winner),...
                    [],session_colors(sessionnum,:))
                scatter(...
                    avg_corr.(strcat(model_pairs{pairnum,1},'_tuningCorr'))(~no_winner),...
                    avg_corr.(strcat(model_pairs{pairnum,2},'_tuningCorr'))(~no_winner),...
                    [],session_colors(sessionnum,:),'filled')
            end
            % make axes pretty
            set(gca,'box','off','tickdir','out',...
                'xlim',[-0.1 1],'ylim',[-0.1 1],...
                'xtick',0:0.5:1,'ytick',0:0.5:1)
            axis square
            if monkeynum ~= 1 || pairnum ~= 1
                set(gca,'box','off','tickdir','out',...
                    'xtick',[],'ytick',[])
            end
            xlabel(sprintf('%s',getModelTitles(model_pairs{pairnum,1})))
            ylabel(sprintf('%s',getModelTitles(model_pairs{pairnum,2})))
        end
    end
    suptitle('Tuning correlation pairwise comparisons')
    saveas(gcf,fullfile(figdir,sprintf('tuningCorr_pairwise_run%s.pdf',run_date)))

    % make the winner dot plot
        % figure('defaultaxesfontsize',18)
        % plotAnalysisWinners(tuning_corr,tuning_corr_winners,struct(...
        %     'monkey_names',{monkey_names},...
        %     'session_ctr',session_ctr,...
        %     'session_colors',session_colors,...
        %     'models_to_plot',{models_to_plot},...
        %     'postfix','_tuningCorr'))
        % suptitle('Tuning Correlation Winners')
        % saveas(gcf,fullfile(figdir,sprintf('tuningCorr_winners_run%s.pdf',run_date)))

    % plot array map of ext model eval
        % put together table
        tuning_corr_table = vertcat(tuning_corr{:});
        tuning_corr_table = horzcat(tuning_corr_table,...
            table(tuning_corr_table.handelbow_tuningCorr-tuning_corr_table.ext_tuningCorr,...
                'VariableNames',{'he_ext_tuningCorr_diff'}));
        tuning_corr_table.Properties.VariableDescriptions{end} = 'linear';
        tuning_corr_table = horzcat(tuning_corr_table,...
            table(tuning_corr_table.ext_tuningCorr-tuning_corr_table.elbow_tuningCorr,...
                'VariableNames',{'ext_elbow_tuningCorr_diff'}));
        tuning_corr_table.Properties.VariableDescriptions{end} = 'linear';
        % average over all crossvals
        tuning_corr_table = neuronAverage(tuning_corr_table,struct('keycols',{{'monkey','date','task','signalID'}},'do_ci',false));
        % extract and add chan
        tuning_corr_table = horzcat(tuning_corr_table,...
            table(tuning_corr_table.signalID(:,1),'VariableNames',{'chan'}));
        % plot
        coarseness = 1;
        figure('defaultaxesfontsize',18)
        lm_table = plotArrayMap(tuning_corr_table,struct(...
            'map_plot','ext_tuningCorr',...
            'calc_linmodels',true,...
            'clims',[0 1],...
            'coarseness',coarseness));
        suptitle(sprintf('Extrinsic tuning correlation (coarseness factor: %d)',coarseness))
        % plotArrayMap(tuning_corr_table,struct('map_plot','ext_elbow_eval_diff','clims',[-0.15 0.15]))
        % suptitle('Hand/Elbow pR2 - Extrinsic pR2')
        for sessionnum = 1:height(lm_table)
            coefTest(lm_table{sessionnum,3}{1})
        end

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
            h = histogram(gca,monkey_shifts_real.velPD*180/pi,'BinWidth',10,'DisplayStyle','stair');
            set(h,'facecolor','none','edgecolor',ones(1,3)*0.5)
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
                h = histogram(gca,monkey_shifts_model.velPD*180/pi,'BinWidth',10,'DisplayStyle','stair');
                set(h,'facecolor','none','edgecolor',ones(1,3)*0.5)
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
    figure(total_hists)
    suptitle('PD shift histograms')
    saveas(gcf,fullfile(figdir,sprintf('PDShiftHists_run%s.pdf',run_date)))
    figure(scatters)
    suptitle('PD shift scatter plots')
    saveas(gcf,fullfile(figdir,sprintf('PDShiftScatters_run%s.pdf',run_date)))

    % Make array maps for dPD
    dPD_map = allFileShifts_real;
    % extract and append channel 
    dPD_map = horzcat(dPD_map,table(dPD_map.signalID(:,1),'VariableNames',{'chan'}));
    % set up color table for dPD
    % dPD_map = horzcat(dPD_map,...
    %     table(interp1(linspace(0,pi,100),viridis(100),abs(dPD_map.velPD)),'VariableNames',{'velPD_color'}));
    dPD_map = horzcat(dPD_map,...
        table(abs(dPD_map.velPD),'VariableNames',{'velPD_color'}));
    % plot array maps
    figure('defaultaxesfontsize',18)
    lm_table = plotArrayMap(dPD_map,struct('map_plot','velPD_color','clims',[0 pi/2],'calc_linmodels',true));
    suptitle('dPD array maps')
    for sessionnum = 1:height(lm_table)
        coefTest(lm_table{sessionnum,3}{1})
    end

%% PD shift VAF dotplots
    % find winners of PD shift
    shift_vaf_winners = cell(length(monkey_names),size(session_colors,1));
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            [shift_vaf_winners{monkeynum,sessionnum},model_pairs] = compareEncoderMetrics(...
                    shift_vaf{monkeynum,sessionnum},struct(...
                        'bonferroni_correction',6,...
                        'models',{models_to_plot},...
                        'model_pairs',{{'ext','handelbow';'extforce','handelbow'}},...
                        'postfix','_vaf'));
        end
    end

    % Find session winners of PD shift
    shift_vaf_session_winners = cell(length(monkey_names),size(session_colors,1));
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            shift_vaf_session = neuronAverage(shift_vaf{monkeynum,sessionnum},struct(...
                'keycols',{{'monkey','date','task','crossvalID'}},...
                'do_ci',false));
            temp_tab = table(ones(height(shift_vaf_session),1),'VariableNames',{'signalID'});
            shift_vaf_session = horzcat(temp_tab,shift_vaf_session);
            [shift_vaf_session_winners{monkeynum,sessionnum},model_pairs] = compareEncoderMetrics(...
                    shift_vaf_session,struct(...
                        'bonferroni_correction',6,...
                        'models',{models_to_plot},...
                        'model_pairs',{{'ext','handelbow';'extforce','handelbow'}},...
                        'postfix','_vaf'));
        end
    end

    % find winners
    all_shift_vaf_winners = horzcat(shift_vaf_winners{:});
    ext_winners = sum(strcmpi(all_shift_vaf_winners,'ext'),2);
    handelbow_winners = sum(strcmpi(all_shift_vaf_winners,'handelbow'),2);
    extforce_winners = sum(strcmpi(all_shift_vaf_winners,'extforce'),2);

    % session winners
    all_shift_vaf_session_winners = horzcat(shift_vaf_session_winners{:});
    ext_winners = sum(strcmpi(all_shift_vaf_session_winners,'ext'),2);
    handelbow_winners = sum(strcmpi(all_shift_vaf_session_winners,'handelbow'),2);
    extforce_winners = sum(strcmpi(all_shift_vaf_session_winners,'extforce'),2);

    % get average shift vaf over all neurons
    shift_vaf_all = vertcat(shift_vaf{:});
    avg_shift_vaf_all = neuronAverage(shift_vaf_all,struct(...
        'keycols',{{'task'}},...
        'do_ci',false));

    % make the winner dot plot
    figure('defaultaxesfontsize',18)
    plotAnalysisWinners(shift_vaf,shift_vaf_winners,struct(...
        'monkey_names',{monkey_names},...
        'session_ctr',session_ctr,...
        'session_colors',session_colors,...
        'models_to_plot',{models_to_plot},...
        'postfix','_vaf'))
    suptitle('PD Shift VAF Winners')
    saveas(gcf,fullfile(figdir,sprintf('PDShiftVAF_winners_run%s.pdf',run_date)))

    % plot session averages with CI bars
    figure('defaultaxesfontsize',18)
    alpha = 0.05;
    model_x = (2:3:((length(models_to_plot)-1)*3+2))/10;
    for monkeynum = 1:length(monkey_names)
        subplot(length(monkey_names),1,monkeynum)
        for sessionnum = 1:session_ctr(monkeynum)
            % plot session average connected by lines...
            avg_shift_vaf = neuronAverage(shift_vaf{monkeynum,sessionnum},...
                struct('keycols',{{'monkey','date','task','crossvalID'}},'do_ci',false));

            % estimate error bars
            [~,cols] = ismember(strcat(models_to_plot,'_vaf'),avg_shift_vaf.Properties.VariableNames);
            num_repeats = double(max(shift_vaf{monkeynum,sessionnum}.crossvalID(:,1)));
            num_folds = double(max(shift_vaf{monkeynum,sessionnum}.crossvalID(:,2)));
            crossval_correction = 1/(num_folds*num_repeats) + 1/(num_folds-1);
            yvals = mean(avg_shift_vaf{:,cols});
            var_vaf = var(avg_shift_vaf{:,cols});
            upp = tinv(1-alpha/2,num_folds*num_repeats-1);
            low = tinv(alpha/2,num_folds*num_repeats-1);
            CI_lo = yvals + low * sqrt(crossval_correction*var_vaf);
            CI_hi = yvals + upp * sqrt(crossval_correction*var_vaf);
            
            % plot dots and lines
            plot(repmat(model_x,2,1),[CI_lo;CI_hi],'-','linewidth',2,'color',session_colors(sessionnum,:))
            hold on
            scatter(model_x(:),yvals(:),50,session_colors(sessionnum,:),'filled')
        end
        ylabel('PD Shift Circular VAF')
        set(gca,'box','off','tickdir','out',...
            'xlim',[0 1.3],'xtick',model_x,'xticklabel',getModelTitles(models_to_plot),...
            'ylim',[0 1],'ytick',[0 1])
    end
    saveas(gcf,fullfile(figdir,sprintf('PDShiftVAF_sessionAvg_run%s.pdf',run_date)))

    % plot array map of ext model shift vaf
        % put together table
        shift_vaf_table = vertcat(shift_vaf{:});
        shift_vaf_table = horzcat(shift_vaf_table,...
            table(shift_vaf_table.handelbow_vaf-shift_vaf_table.ext_vaf,...
                'VariableNames',{'he_ext_vaf_diff'}));
        shift_vaf_table.Properties.VariableDescriptions{end} = 'linear';
        shift_vaf_table = horzcat(shift_vaf_table,...
            table(shift_vaf_table.ext_vaf-shift_vaf_table.elbow_vaf,...
                'VariableNames',{'ext_elbow_vaf_diff'}));
        shift_vaf_table.Properties.VariableDescriptions{end} = 'linear';
        % average over all crossvals
        shift_vaf_table = neuronAverage(shift_vaf_table,struct('keycols',{{'monkey','date','task','signalID'}},'do_ci',false));
        % extract and add chan
        shift_vaf_table = horzcat(shift_vaf_table,...
            table(shift_vaf_table.signalID(:,1),'VariableNames',{'chan'}));
        % plot
        coarseness = 2;
        figure('defaultaxesfontsize',18)
        lm_table = plotArrayMap(shift_vaf_table,struct(...
            'map_plot','ext_elbow_vaf_diff',...
            'calc_linmodels',true,...
            'coarseness',coarseness));
        suptitle(sprintf('Extrinsic shift VAF (coarseness factor: %d)',coarseness))
        % plotArrayMap(tuning_corr_table,struct('map_plot','ext_elbow_eval_diff','clims',[-0.15 0.15]))
        % suptitle('Hand/Elbow pR2 - Extrinsic pR2')
        for sessionnum = 1:height(lm_table)
            coefTest(lm_table{sessionnum,3}{1})
        end

%% Get example tuning curves for all models
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            %% Plot out tuning curves for tuned neurons
            f = figure('defaultaxesfontsize',18);
            % plot a max of 7 neurons
            num_neurons = min(7,size(tuned_neurons{monkeynum,sessionnum},1));
            neurons_to_plot = tuned_neurons{monkeynum,sessionnum}(...
                randperm(size(tuned_neurons{monkeynum,sessionnum},1),num_neurons),...
                :);
            for neuronnum = 1:num_neurons
                % figure out maxFR over both workspaces
                [~,temp_table] = getNTidx(model_tuning{monkeynum,sessionnum},...
                    'signalID',neurons_to_plot(neuronnum,:));
                minFR = floor(min(min(temp_table{:,strcat([models_to_plot {'S1_FR'}],'_velCurve')})));
                maxFR = ceil(max(max(temp_table{:,strcat([models_to_plot {'S1_FR'}],'_velCurve')})));

                % go plot each workspace
                for spacenum = 1:size(cond_colors,1)
                    % get tuning table specifically for this neuron and workspace
                    [~,temp_table] = getNTidx(model_tuning{monkeynum,sessionnum},...
                        'signalID',neurons_to_plot(neuronnum,:),...
                        'spaceNum',spacenum);

                    % plot out actual tuning curves
                    subplot(length(models_to_plot)+1,num_neurons,neuronnum)
                    plotTuning(temp_table,...
                        struct('maxFR',maxFR,...
                            'minFR',minFR,...
                            'unroll',true,...
                            'color',cond_colors(spacenum,:),...
                            'curve_colname',sprintf('%s_velCurve','S1_FR'),...
                            'pd_colname',sprintf('%s_velPD','S1_FR'),...
                            'plot_ci',false))
                    title(strcat('Neuron ',num2str(neurons_to_plot(neuronnum,:))))

                    % plot out modeled tuning curves
                    for modelnum = 1:length(models_to_plot)
                        subplot(length(models_to_plot)+1,num_neurons,modelnum*num_neurons+neuronnum)
                        plotTuning(temp_table,...
                            struct('maxFR',maxFR,...
                                'minFR',minFR,...
                                'unroll',true,...
                                'color',cond_colors(spacenum,:),...
                                'curve_colname',sprintf('%s_velCurve',models_to_plot{modelnum}),...
                                'pd_colname',sprintf('%s_velPD',models_to_plot{modelnum}),...
                                'plot_ci',false))
                        title(getModelTitles(models_to_plot{modelnum}),'interpreter','none')
                    end
                end
            end
            set(gcf,'renderer','Painters')
            saveas(gcf,fullfile(figdir,sprintf('tuningCurves_%s_%s_run%s.pdf',...
                model_tuning{monkeynum,sessionnum}{1,'monkey'}{1},...
                strrep(model_tuning{monkeynum,sessionnum}{1,'date'}{1},'/',''),...
                run_date)))
            % close(f)
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

%% example raster plot
    % load data
    filenum = 4;
    load(fullfile(datadir,filename{filenum}))

    % set plotting params
    num_trials = 2;
    trials_to_plot = randperm(length(encoderResults.td_tuning{2}),num_trials);

    for neuronnum = 1:size(encoderResults.td_tuning{2}(1).S1_spikes,2)
        h = figure('defaultaxesfontsize',18);
        ax = zeros(length(models_to_plot),2);
        for spacenum = 1:2
            ax(spacenum) = subplot(2,1,spacenum);
            plotExampleRaster(encoderResults.td_tuning{spacenum},...
                struct('trial_idx',trials_to_plot))
            title(sprintf('Neuron %d, spacenum %d',neuronnum,spacenum))
        end
        linkaxes(ax(:),'y')
        waitfor(h)
    end

%% Within model class comparison
    % figure out what comparisons are allowed
    allowed_comparisons = vertcat(nchoosek(hand_models,2),nchoosek(arm_models,2));
    allowed_comparisons = sort(allowed_comparisons')';

    % pR2 comparison for all models
        all_pr2_winners = cell(length(monkey_names),size(session_colors,1));
        for monkeynum = 1:length(monkey_names)
            for sessionnum = 1:session_ctr(monkeynum)
                [all_pr2_winners{monkeynum,sessionnum},model_pairs] = compareEncoderMetrics(...
                        model_eval{monkeynum,sessionnum},struct(...
                            'bonferroni_correction',6,...
                            'models',{included_models},...
                            'model_pairs',{allowed_comparisons},...
                            'postfix','_eval'));
            end
        end
        figure
        for monkeynum = 1:length(monkey_names)
            for pairnum = 1:size(allowed_comparisons,1)
                % set subplot
                subplot(length(monkey_names),size(allowed_comparisons,1),...
                    (monkeynum-1)*size(model_pairs,1)+pairnum)
                plot([-1 1],[-1 1],'k--','linewidth',0.5)
                hold on
                plot([0 0],[-1 1],'k-','linewidth',0.5)
                plot([-1 1],[0 0],'k-','linewidth',0.5)
                for sessionnum = 1:session_ctr(monkeynum)
                    avg_pR2 = neuronAverage(model_eval{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
                    % scatter filled circles if there's a winner, empty circles if not
                    no_winner =  cellfun(@isempty,all_pr2_winners{monkeynum,sessionnum}(pairnum,:));
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
                if monkeynum ~= 1 || pairnum ~= 1
                    set(gca,'box','off','tickdir','out',...
                        'xtick',[],'ytick',[])
                end
                xlabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,1})))
                ylabel(sprintf('%s pR2',getModelTitles(model_pairs{pairnum,2})))
            end
        end
        suptitle('Pseudo-R^2 pairwise comparisons')
        saveas(gcf,fullfile(figdir,sprintf('allmodel_pr2_pairwise_run%s.pdf',run_date)))

    % tuning correlation comparison all models
        all_tuning_corr_winners = cell(length(monkey_names),size(session_colors,1));
        for monkeynum = 1:length(monkey_names)
            for sessionnum = 1:session_ctr(monkeynum)
                [all_tuning_corr_winners{monkeynum,sessionnum},model_pairs] = compareEncoderMetrics(...
                        tuning_corr{monkeynum,sessionnum},struct(...
                            'bonferroni_correction',6,...
                            'models',{included_models},...
                            'model_pairs',{allowed_comparisons},...
                            'postfix','_tuningCorr'));
            end
        end
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
                    avg_corr = neuronAverage(tuning_corr{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
                    % scatter filled circles if there's a winner, empty circles if not
                    no_winner =  cellfun(@isempty,all_tuning_corr_winners{monkeynum,sessionnum}(pairnum,:));
                    scatter(...
                        avg_corr.(strcat(model_pairs{pairnum,1},'_tuningCorr'))(no_winner),...
                        avg_corr.(strcat(model_pairs{pairnum,2},'_tuningCorr'))(no_winner),...
                        [],session_colors(sessionnum,:))
                    scatter(...
                        avg_corr.(strcat(model_pairs{pairnum,1},'_tuningCorr'))(~no_winner),...
                        avg_corr.(strcat(model_pairs{pairnum,2},'_tuningCorr'))(~no_winner),...
                        [],session_colors(sessionnum,:),'filled')
                end
                % make axes pretty
                set(gca,'box','off','tickdir','out',...
                    'xlim',[-0.1 1],'ylim',[-0.1 1],...
                    'xtick',0:0.5:1,'ytick',0:0.5:1)
                axis square
                if monkeynum ~= 1 || pairnum ~= 1
                    set(gca,'box','off','tickdir','out',...
                        'xtick',[],'ytick',[])
                end
                xlabel(sprintf('%s',getModelTitles(model_pairs{pairnum,1})))
                ylabel(sprintf('%s',getModelTitles(model_pairs{pairnum,2})))
            end
        end
        suptitle('Tuning correlation pairwise comparisons')
        saveas(gcf,fullfile(figdir,sprintf('allmodel_tuningCorr_pairwise_run%s.pdf',run_date)))

    % PD shift for non-plotted models
        file_shifts = cell(length(filename),length(not_plot_models)); % shift tables for each model in each file
        for filenum = 1:length(filename)
            % load data
            load(fullfile(datadir,filename{filenum}))

            shift_tables = calculatePDShiftTables(encoderResults,[strcat('glm_',not_plot_models,'_model') 'S1_FR']);
            mean_shifts = cell(length(not_plot_models),1);
            for modelnum = 1:length(not_plot_models)+1
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
                subplot(length(monkey_names),length(not_plot_models)+1,(monkeynum-1)*(length(not_plot_models)+1)+1)
                h = histogram(gca,monkey_shifts_real.velPD*180/pi,'BinWidth',10,'DisplayStyle','stair');
                set(h,'facecolor','none','edgecolor',ones(1,3)*0.5)
                set(gca,'box','off','tickdir','out','xlim',[-180 180],'xtick',[-180 0 180],'ylim',[0 ylim_high],'ytick',[0 ylim_high/2 ylim_high],'view',[-90 90])
                ylabel(monkey_names{monkeynum})
                if monkeynum == 1
                    title('Actual PD Shift')
                end
                for modelnum = 1:length(not_plot_models)
                    allFileShifts_model = vertcat(file_shifts{:,modelnum});
                    [~,monkey_shifts_model] = getNTidx(allFileShifts_model,'monkey',monkey_names{monkeynum});

                    % modeled PD shift histogram
                    subplot(length(monkey_names),length(not_plot_models)+1,(monkeynum-1)*(length(not_plot_models)+1)+modelnum+1)
                    h = histogram(gca,monkey_shifts_model.velPD*180/pi,'BinWidth',10,'DisplayStyle','stair');
                    set(h,'facecolor','none','edgecolor',ones(1,3)*0.5)
                    set(gca,'box','off','tickdir','out','xlim',[-180 180],'xtick',[-180 0 180],'ylim',[0 ylim_high],'ytick',[0 ylim_high/2 ylim_high],'view',[-90 90])
                    if monkeynum == 1
                        title(sprintf('%s modeled PD shift',getModelTitles(not_plot_models{modelnum})))
                    end
                end

            % make plots separating sessions
                for sessionnum = 1:length(session_dates)
                    % get real shifts for this session
                    [~,session_shifts_real] = getNTidx(allFileShifts_real,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});

                    % first the real PD shift histogram
                    % figure(hists)
                    % subplot(length(monkey_names),length(not_plot_models)+1,(monkeynum-1)*(length(not_plot_models)+1)+1)
                    % h = histogram(gca,session_shifts_real.velPD*180/pi,'BinWidth',10,'DisplayStyle','bar');
                    % set(h,'facecolor',session_colors(sessionnum,:),'edgecolor','none')
                    % hold on

                    % now the models
                    for modelnum = 1:length(not_plot_models)
                        % get the modeled shifts for this session
                        allFileShifts_model = vertcat(file_shifts{:,modelnum});
                        [~,session_shifts_model] = getNTidx(allFileShifts_model,'monkey',monkey_names{monkeynum},'date',session_dates{sessionnum});

                        % % modeled PD shift histogram
                        % figure(hists)
                        % subplot(length(monkey_names),length(not_plot_models)+1,(monkeynum-1)*(length(not_plot_models)+1)+modelnum+1)
                        % h = histogram(gca,session_shifts_model.velPD*180/pi,'BinWidth',10,'DisplayStyle','bar');
                        % set(h,'facecolor',session_colors(sessionnum,:),'edgecolor','none')
                        % hold on

                        % scatter plots
                        figure(scatters)
                        subplot(length(monkey_names),length(not_plot_models),(monkeynum-1)*(length(not_plot_models))+modelnum)
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
                % subplot(length(monkey_names),length(not_plot_models)+1,(monkeynum-1)*(length(not_plot_models)+1)+1)
                % set(gca,'box','off','tickdir','out','xlim',[-180 180],'xtick',[-180 0 180],'ylim',[0 15],'ytick',[0 15 30],'view',[-90 90])
                % ylabel(monkey_names{monkeynum})
                % if monkeynum == 1
                %     title('Actual PD Shift')
                % end
                for modelnum = 1:length(not_plot_models)
                    % histograms
                    % figure(hists)
                    % subplot(length(monkey_names),length(not_plot_models)+1,(monkeynum-1)*(length(not_plot_models)+1)+modelnum+1)
                    % set(gca,'box','off','tickdir','out','xlim',[-180 180],'xtick',[-180 0 180],'ylim',[0 15],'ytick',[0 15 30],'view',[-90 90])
                    % if monkeynum == 1
                    %     title(sprintf('%s modeled PD shift',getModelTitles(not_plot_models{modelnum})))
                    % end

                    % scatter plots
                    figure(scatters)
                    subplot(length(monkey_names),length(not_plot_models),(monkeynum-1)*length(not_plot_models)+modelnum)
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
                        title(sprintf('%s model',getModelTitles(not_plot_models{modelnum})),'interpreter','none')
                    end
                end
        end
        figure(total_hists)
        suptitle('PD shift histograms')
        saveas(gcf,fullfile(figdir,sprintf('extra_PDShiftHists_run%s.pdf',run_date)))
        figure(scatters)
        suptitle('PD shift scatter plots')
        saveas(gcf,fullfile(figdir,sprintf('extra_PDShiftScatters_run%s.pdf',run_date)))

%% Extra stuff/in progress...
    %% Collect model prediction differences and plot
    for filenum = 1%:length(filename)
        load(fullfile(datadir,filename{filenum}))

        td_tuning = horzcat(encoderResults.td_tuning{:});
        td_tuning = binTD(td_tuning,'average');

        % get difference between hand/elbow and extrinsic model predictions
        % handelbow_ext_pred_diff = getSig(td_tuning,'glm_handelbow_model')-getSig(td_tuning,'glm_ext_model');
        handelbow_ext_pred_diff = vertcat(td_tuning.glm_handelbow_model)-vertcat(td_tuning.glm_ext_model);

        % normalize by z-scoring
        handelbow_ext_pred_diff = zscore(handelbow_ext_pred_diff);

        % plot image
        figure
        imagesc((handelbow_ext_pred_diff'),[-10 10])
        colormap(viridis)
        colorbar

        figure
        plot(abs(handelbow_ext_pred_diff))
    end
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
