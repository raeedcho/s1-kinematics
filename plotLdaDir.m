% This script plots the lda coefficients in the full neural space.
% Each line represents the set of coefficients for a particular
% cross-validation run.

%% Set up variables
    datadir = '/data/raeed/project-data/limblab/s1-kinematics/Results/Separability';
    % filename = {'Han_20171101_TRT_encodingResults_run20180809.mat','Chips_20170915_TRT_encodingResults_run20180809.mat','Lando_20170802_encodingResults_run20180809.mat'};
    files = dir(fullfile(datadir,'*separationResults_run20190228.mat'));
    filename = horzcat({files.name});
    
    % for figure saving
    figdir = '/home/raeed/Wiki/Projects/limblab/s1-kinematics/figures/Encoding';
    run_date = char(datetime('today','format','yyyyMMdd'));
    
    monkey_names = {'Chips','Han'};
    % models_to_plot = {'S1_FR','ext','extforce','handelbow'};
    % fr_names = {'S1_FR','ext_predFR','extforce_predFR','handelbow_predFR'};
    % model_titles = {'Actual Firing','Extrinsic','Extrinsic + Force','Hand/Elbow'};
    % num_pcs = 3;


%% Load and plot lda from each file
    for filenum = 1:length(filename)
        % load data
        load(fullfile(datadir,filename{filenum}))

        % extract...stuff
        true_lda_coeff = findLdaDir(sepResults.lda_table);

        % plot lda coefficiencts
        figure
        for entrynum = 1:height(true_lda_coeff)
            plot(true_lda_coeff(entrynum,:).S1_FR_lda_full_coeff,'ko')
            % bar(true_lda_coeff(entrynum,:).S1_FR_lda_full_coeff)
            hold on
        end
        plot(xlim,[0 0],'k-')
        
        set(gca,'box','off','tickdir','out');

        title(filename{filenum})
    end

