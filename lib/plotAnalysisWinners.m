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

            % set x counter and potential y values
            model_xval = zeros(length(models_to_plot)+1,1);
            model_yval = model_y + template_y(sessionnum);
            for neuronnum = 1:height(avg_metric)

                % find highest pR2
                [~,model_sort_idx] = sort(avg_metric{neuronnum,model_cols});
                model_order = model_colnames(model_sort_idx);
                best_model = model_order{end};
                runnerup_model = model_order{end-1};

                % check decisiveness of victory
                best_model_wins = strcmpi(winners{monkeynum,sessionnum}(:,neuronnum),best_model);
                runnerup_wins = strcmpi(winners{monkeynum,sessionnum}(:,neuronnum),runnerup_model);
                if sum(best_model_wins) == length(models_to_plot)-1
                    yval_plot = model_yval(model_sort_idx(end));
                    model_xval(model_sort_idx(end)) = model_xval(model_sort_idx(end)) + 1;
                    xval_plot = model_xval(model_sort_idx(end));
                    scatter(repmat(xval_plot,size(yval_plot)),yval_plot,[],session_colors(sessionnum,:),'filled')
                elseif sum(best_model_wins) == length(models_to_plot)-2
                    yval_plot = model_yval(model_sort_idx(end));
                    model_xval(model_sort_idx(end)) = model_xval(model_sort_idx(end)) + 1;
                    xval_plot = model_xval(model_sort_idx(end));
                    scatter(repmat(xval_plot,size(yval_plot)),yval_plot,[],session_colors(sessionnum,:),'<','filled')

                    if sum(runnerup_wins) == length(models_to_plot)-2
                        yval_plot = model_yval(model_sort_idx(end-1));
                        model_xval(model_sort_idx(end-1)) = model_xval(model_sort_idx(end-1)) + 1;
                        xval_plot = model_xval(model_sort_idx(end-1));
                        scatter(repmat(xval_plot,size(yval_plot)),yval_plot,[],session_colors(sessionnum,:),'<','filled')
                    end
                else
                    yval_plot = model_yval(end);
                    model_xval(end) = model_xval(end) + 1;
                    xval_plot = model_xval(end);
                    scatter(repmat(xval_plot,size(yval_plot)),yval_plot,[],session_colors(sessionnum,:),'filled')
                end

                hold on
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
end
