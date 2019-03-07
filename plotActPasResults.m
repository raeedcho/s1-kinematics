%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script makes plots from results saved by calculateActPasSeparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up plotting variables
    datadir = '/home/raeed/data/project-data/limblab/s1-kinematics/Results/Separability';
    % filename = {'Han_20171101_TRT_encodingResults_run20180809.mat','Chips_20170915_TRT_encodingResults_run20180809.mat','Lando_20170802_encodingResults_run20180809.mat'};
    files = dir(fullfile(datadir,'*separationResults_run20190228.mat'));
    filename = horzcat({files.name});
    
    % for figure saving
    figdir = '/home/raeed/Wiki/Projects/limblab/s1-kinematics/figures/Encoding';
    run_date = char(datetime('today','format','yyyyMMdd'));

    monkey_names = {'Chips','Han'};
    models_to_plot = {'S1_FR','ext','extforce','handelbow'};
    fr_names = {'S1_FR','ext_predFR','extforce_predFR','handelbow_predFR'};
    model_titles = {'Actual Firing','Extrinsic','Extrinsic + Force','Hand/Elbow'};
    num_pcs = 3;

    session_colors = [...
        102,194,165;...
        252,141,98;...
        141,160,203]/255;

%% Compile information over all files
    [sep_table_cell] = deal(cell(length(filename),1));
    fileclock = tic;
    fprintf('Started loading files...\n')
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        % extract...stuff
        meta_cols = strcmpi(sepResults.lda_table.Properties.VariableDescriptions,'meta');
        sep_cols = endsWith(sepResults.lda_table.Properties.VariableNames,'true_sep');
        sep_table_cell{filenum} = sepResults.lda_table(:,meta_cols | sep_cols);

        % compose trial table for one crossval run
        repeatnum = 3;
        num_folds = max(sepResults.trial_table.crossvalID(:,2));
        trial_table_cell = cell(num_folds,1);
        for foldnum = 1:num_folds
            [~,trial_table_cell{foldnum}] = getNTidx(sepResults.trial_table,'crossvalID',[repeatnum foldnum]);
        end
        trial_table = vertcat(trial_table_cell{:});

        % plot out neural population scatter
        isActive = ~trial_table.isPassive;
        [dirs,~,dir_idx] = unique(trial_table.trialDir);
        dir_colors = linspecer(length(dirs));
        figure('defaultaxesfontsize',18)
        for modelnum = 1:length(models_to_plot)
            subplot(1,length(models_to_plot),modelnum)

            model_fr = trial_table.(fr_names{modelnum});
            [~,model_pca] = pca(model_fr);
            model_pca = model_pca(:,1:num_pcs);
            lda_mdl = fitcdiscr(model_pca,isActive);
            lda_vec = lda_mdl.Coeffs(2,1).Linear;
            lda_vec = lda_vec/sqrt(sum(lda_vec.^2)); % make into unit vector
            null_basis = null(lda_vec');

            % make 3d plot
            % scatter3(...
            %     model_pca(isActive,:)*lda_vec,...
            %     model_pca(isActive,:)*null_basis(:,1),...
            %     model_pca(isActive,:)*null_basis(:,2),...
            %     [],dir_colors(dir_idx(isActive),:),'filled')
            % hold on
            % scatter3(...
            %     model_pca(~isActive,:)*lda_vec,...
            %     model_pca(~isActive,:)*null_basis(:,1),...
            %     model_pca(~isActive,:)*null_basis(:,2),...
            %     [],dir_colors(dir_idx(~isActive),:))
            % % plot lines
            % plot3([0 0],ylim,[0 0],'--k','linewidth',2)
            % plot3([0 0],[0 0],zlim,'--k','linewidth',2)
            % axis equal

            % 2D plot
            % scatter(...
            %     model_pca(isActive,:)*lda_vec,...
            %     model_pca(isActive,:)*null_basis(:,1),...
            %     [],dir_colors(dir_idx(isActive),:),'filled')
            % hold on
            % scatter(...
            %     model_pca(~isActive,:)*lda_vec,...
            %     model_pca(~isActive,:)*null_basis(:,1),...
            %     [],dir_colors(dir_idx(~isActive),:))
            % plot([0 0],ylim,'--k','linewidth',2)
            % axis equal

            % 2D plot with no dir color
            scatter(...
                model_pca(isActive,:)*lda_vec,...
                model_pca(isActive,:)*null_basis(:,1),...
                [],'k','filled')
            hold on
            scatter(...
                model_pca(~isActive,:)*lda_vec,...
                model_pca(~isActive,:)*null_basis(:,1),...
                [],'r','filled')
            plot([0 0],ylim,'--k','linewidth',2)
            axis equal
            axis off

            title(model_titles{modelnum})
        end
        suptitle(sprintf('%s-%s',trial_table.monkey{1},trial_table.date_time{1}))

        % output a counter
        fprintf('Processed file %d of %d at time %f\n',filenum,length(filename),toc(fileclock))
    end
    sep_table = vertcat(sep_table_cell{:});

%% Make scatter
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
            
            % plot dots and lines
            plot(model_x',yvals','-','linewidth',0.5,'color',ones(1,3)*0.5)
            hold on
            plot(repmat(model_x,2,1),[CI_lo;CI_hi],'-','linewidth',2,'color',session_colors(sessionnum,:))
            scatter(model_x(:),yvals(:),50,session_colors(sessionnum,:),'filled')
        end
        plot([0 1.3],[0.5 0.5],'--k','linewidth',2)
        ylabel('Separability (%)')
        set(gca,'box','off','tickdir','out',...
            'xlim',[0 1.3],'xtick',model_x,'xticklabel',model_titles,...
            'ylim',[0 1],'ytick',[0 0.5 1])
    end
    % saveas(gcf,fullfile(figdir,sprintf('actpasSeparability_run%s.pdf',run_date)))

