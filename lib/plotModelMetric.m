function plotModelMetric(metric_tables,params)
    % plots model metric (pr2, tuning curve correlation, or PD shift error)
    % by either neuron or by crossval run
    % Inputs:
    %   metric_table - cell array of tables of metrics for models to plot
    %   params - parameters struct
    %       .monkey_names - names of monkeys
    %       .session_ctr - array matching monkey_names with numbers of sessions
    %       .session_colors - colors for sessions
    %       .models_to_plot - models to plot
    %       .postfix - postfix of model name in table column names,
    %           e.g. '_eval' for column name 'ext_eval' and model name 'ext'
    %       .marginal_col - key column(s) for which to average, e.g. 'signalID'
    %           or 'crossvalID' or {'signalID','crossvalID'}
    %       .line_sparsity - percent of lines to drop in plot

    % parameters
    monkey_names = {};
    session_ctr = [];
    session_colors = [];
    models_to_plot = {};
    postfix = '';
    marginal_col = 'signalID';
    colored_lines = false;
    line_sparsity = 0;

    assignParams(who,params)

    if ischar(marginal_col)
        marginal_col = {marginal_col};
    end
    assert(all(endsWith(marginal_col,'ID')),'Marginal col must be signalID, crossvalID or both!')

    % y coordinate of individual monkey bars
    monkey_y = (2:3:((length(monkey_names)-1)*3+2))/10;
    % template for within monkey bars separation
    template_y = linspace(-1,1,length(models_to_plot))/10;
    for monkeynum = 1:length(monkey_names)
        for sessionnum = 1:session_ctr(monkeynum)
            % average for each neuron
            avg_metric = neuronAverage(metric_tables{monkeynum,sessionnum},...
                struct('keycols',{[{'monkey','date','task'},marginal_col]},'do_ci',false));
            yval = repmat(monkey_y(monkeynum) + template_y,height(avg_metric),1);
            % add some jitter
            yval = yval+randn(size(yval,1),1)/150;

            % sparsify the lines if needed (not for neurons...)
            doplot = rand(length(yval),1)>=line_sparsity;
            [~,cols] = ismember(strcat(models_to_plot,postfix),avg_metric.Properties.VariableNames);
            xvals = avg_metric{:,cols};
            if ~colored_lines
                plot(xvals(doplot,:)',yval(doplot,:)','-','linewidth',0.5,'color',ones(1,3)*0.5)
                hold on
                scatter(xvals(:),yval(:),50,session_colors(sessionnum,:),'filled')
            else
                plot(xvals(doplot,:)',yval(doplot,:)','-','linewidth',0.5,'color',session_colors(sessionnum,:))
                hold on
                scatter(xvals(:),yval(:),25,ones(1,3)*0.5,'filled')
            end
        end
    end
    axis ij
    ytickmarks = monkey_y + template_y';
    set(gca,'box','off','tickdir','out',...  'xlim',[0,1],'xtick',0:0.5:1.0,...
        'ytick',ytickmarks(:),'yticklabel',repmat(getModelTitles(models_to_plot),1,length(monkey_names)))
    ylabel(vertcat(monkey_names(:)))
    ylbl = get(gca,'ylabel');
    set(ylbl,'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','center')
