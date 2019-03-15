function plotAnalysisWinners(metric_tables,winners,params)
    % plotAnalysisWinners - Make a dotplot of "analysis winners"
    % Each model gets a dot if it wins all pairwise comparisons
    % with other models and gets a half dot if it wins all but
    % one comparison.
    % Inputs:
    %   metric_tables - cell array of tables of metrics to compare models on
    %       Each cell corresponds to a metric table for a monkey session
    %   winners - cell array (corresponding to metric_tables) of winners found
    %       by compareEncoderMetrics
    %   params - parameters struct
    %       .monkey_names - names of monkeys
    %       .session_ctr - array matching monkey_names with number of sessions
    %       .session_colors - colors for sessions
    %       .models_to_plot - models to plot
    %       .postfix - postfix of model name in table column names,
    %           e.g. '_eval' for column name 'ext_eval' and model name 'ext'

    % parameters
    monkey_names = {};
    session_ctr = [];
    session_colors = [];
    models_to_plot = {};
    postfix = '';
    assignParams(who,params)

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
            avg_metric = neuronAverage(metric_tables{monkeynum,sessionnum},struct('keycols','signalID','do_ci',false));
            model_cols = endsWith(avg_metric.Properties.VariableNames,postfix);
            model_colnames = strrep(avg_metric.Properties.VariableNames(model_cols),postfix,'');

            assert(all(strcmpi(model_colnames,models_to_plot)),'Models in table are in different order than specified!')

            % count all full wins and partial wins for this session
            full_win_ctr = zeros(length(models_to_plot)+1,1);
            part_win_ctr = zeros(length(models_to_plot)+1,1);
            for neuronnum = 1:height(avg_metric)
                % find highest pR2
                [~,model_sort_idx] = sort(avg_metric{neuronnum,model_cols});
                best_idx = model_sort_idx(end);
                runnerup_idx = model_sort_idx(end-1);
                best_model = model_colnames{best_idx};
                runnerup_model = model_colnames{runnerup_idx};

                % find number of wins
                best_model_wins = strcmpi(winners{monkeynum,sessionnum}(:,neuronnum),best_model);
                runnerup_wins = strcmpi(winners{monkeynum,sessionnum}(:,neuronnum),runnerup_model);
                if sum(best_model_wins) == length(models_to_plot)-1 % won all comparisons
                    full_win_ctr(best_idx) = full_win_ctr(best_idx)+1;
                elseif sum(best_model_wins) == length(models_to_plot)-2 % won all but one comparison
                    part_win_ctr(best_idx) = part_win_ctr(best_idx)+1;
                    if sum(runnerup_wins) == length(models_to_plot)-2 % there's another partial winner
                        part_win_ctr(runnerup_idx) = part_win_ctr(runnerup_idx)+1;
                    end
                else % no full or partial winners
                    full_win_ctr(end) = full_win_ctr(end)+1;
                end
                hold on
            end

            % plot dots
            % full wins first
            dot_xval = cell(1,length(full_win_ctr));
            dot_yval = cell(1,length(full_win_ctr));
            for modelnum = 1:length(full_win_ctr)
                dot_xval{modelnum} = 1:full_win_ctr(modelnum);
                dot_yval{modelnum} = repmat(model_y(modelnum) + template_y(sessionnum),1,full_win_ctr(modelnum));
            end
            dot_xval = horzcat(dot_xval{:});
            dot_yval = horzcat(dot_yval{:});
            scatter(dot_xval,dot_yval,[],session_colors(sessionnum,:),'filled')

            % then part wins
            dot_xval = cell(1,length(part_win_ctr));
            dot_yval = cell(1,length(part_win_ctr));
            for modelnum = 1:length(part_win_ctr)
                dot_xval{modelnum} = full_win_ctr(modelnum)+(1:part_win_ctr(modelnum));
                dot_yval{modelnum} = repmat(model_y(modelnum) + template_y(sessionnum),1,part_win_ctr(modelnum));
            end
            dot_xval = horzcat(dot_xval{:});
            dot_yval = horzcat(dot_yval{:});
            scatter(dot_xval,dot_yval,[],session_colors(sessionnum,:))
        end
        title(monkey_names{monkeynum})
        axis ij
        if monkeynum == 1
            set(gca,'box','off','tickdir','out',...
                'ylim',[model_y(1)+template_y(1)-0.1 model_y(end)+template_y(end)+0.1],...
                'ytick',model_y,...
                'xlim',[0 50],...
                'xtick',0:10:50,...
                'yticklabel',[getModelTitles(models_to_plot);{'None'}])
        else
            set(gca,'box','off','tickdir','out',...
                'ylim',[model_y(1)+template_y(1)-0.1 model_y(end)+template_y(end)+0.1],...
                'ytick',model_y,...
                'xlim',[0 50],...
                'xtick',0:10:50,...
                'yticklabel',{})
        end
        xlabel('Number of neurons')
    end
end
